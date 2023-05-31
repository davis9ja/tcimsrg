#include <stdio.h>
#include <stdlib.h>
#include <ctime>
#include <boost/numeric/odeint.hpp>
#include <boost/format.hpp>
#include <chrono>
#include <iostream>
#include <fstream>
#include <omp.h>

#include "taco.h"
#include "dmcpp.hpp"
#include "pairinghamiltonian.hpp"
#include "occupation_factors.hpp"
#include "generator.hpp"
#include "flow_imsrg2.hpp"
#include "system.hpp"
#include "DMD.h"

#include "white.cpp"
#include "white_atan.cpp"
#include "brillouin.cpp"
#include "BACKEND_ublas.cpp"
#include "BACKEND_taco.cpp"
#include "DMD.c"

//#include "derivative.hpp"

void process_user_input(int argc, char **argv, 
                        int &numStates, 
                        int &nHoles, int &nParticles,
                        double &d, double &g, double &pb, 
                        double &t0, double &t1, double &dt,
                        int &generator_id,
                        int &reference_type,
                        int &solver, int &BACKEND_id,
                        int &dmdOn, int &dmdObs, int &dmdTrunc) {

    // SET DEFAULT PARAMS
    numStates = 8;
    nHoles = 4;
    nParticles = 4;
    d = 1.0;
    g = 0.5;
    pb = 0.0;
    t0 = 0.0;
    t1 = 5.0;
    dt = 0.1;
    generator_id = 0;
    reference_type = 0;
    solver = 0;
    BACKEND_id = 1;
    dmdOn = 0;
    dmdObs = 20;
    dmdTrunc = 4;
    
    std::cout << "Running IMSRG(2) on pairing-plus-particle-hole model\n";
    if (argc > 1) {
        for (int i = 1; i < argc; i++) {
            size_t pos = 0;
            std::string delim = "=";
            std::string line(argv[i]);
            std::string token, param_str, param_val;

            pos = line.find(delim);
            token = line.substr(0, pos);
            line.erase(0, pos + delim.length());
            
            param_str = token;
            param_val = line;

            if (param_str == "numStates")
                numStates = std::stoi(param_val);
            else if (param_str == "nHoles")
                nHoles = std::stoi(param_val);
            else if (param_str == "nParticles")
                nParticles = std::stoi(param_val);
            else if (param_str == "d")
                d = std::stof(param_val);
            else if (param_str == "g")
                g = std::stof(param_val);
            else if (param_str == "pb")
                pb = std::stof(param_val);
            else if (param_str == "t0")
                t0 = std::stof(param_val);
            else if (param_str == "t1")
                t1 = std::stof(param_val);
            else if (param_str == "dt")
                dt = std::stof(param_val);
            else if (param_str == "generator_id")
                generator_id = std::stof(param_val);
            else if (param_str == "reference_type")
                reference_type = std::stof(param_val);
            else if (param_str == "solver")
                solver = std::stof(param_val);
            else if (param_str == "BACKEND_id")
                BACKEND_id = std::stof(param_val);
            else if (param_str == "dmdOn")
                dmdOn = std::stoi(param_val);
            else if (param_str == "dmdObs")
                dmdObs = std::stoi(param_val);
            else if (param_str == "dmdTrunc")
                dmdTrunc = std::stoi(param_val);
            else
                std::cout << param_str << " not recognized as input parameter.\n" << std::endl;
        }
    }

    if (reference_type == 1)
        printf("Parameters:\n\t# s.p. states = %d\n\td=%f\n\tg=%f\n\tpb=%f\n", numStates, d, g, pb);
    else if (reference_type == 0) {
        printf("Parameters:\n\t# holes = %d\n\t# particles = %d\n\td=%0.3f\n\tg=%0.3f\n\tpb=%0.3f\n", nHoles, nParticles, d, g, pb);
        numStates = nHoles + nParticles;
    }
    else {
        std::cout << "REFERENCE TYPE NOT RECOGNIZED. EXITING...";
        exit(1);
    }
        
    std::cout << "IMSRG generator = ";
    if (generator_id == 0) {
        std::cout << "White\n";
    } else if (generator_id == 1) {
        std::cout << "White (atan)\n";
    } else if (generator_id == 2) {
        std::cout << "Brillouin (MR)\n";
    } else {
        std::cout << "GENERATOR INPUT NOT RECOGNIZED. EXITING...";
        exit(1);
    }

    std::cout << "ODE Stepper = ";    
    if (solver == 0) {
        std::cout << "Adams-Bashforth-Moulton\n";
    } else if (solver == 1) {
        std::cout << "Runge-Kutta 4 point\n";
    } else {
        std::cout << "STEPPER DOES NOT EXIST. EXITING...";
        exit(1);            
    }

    printf("Solver range t0=%0.3f, t1=%0.3f, dt=%0.3f\n", t0, t1, dt);

    std::cout << "BACKEND flow computer = ";
    if (BACKEND_id == 0) {
        std::cout << "uBLAS vector loops\n";
    } else if (BACKEND_id == 1) {
        std::cout << "TACO tensor compiler\n";
    } else {
        std::cout << "FLOW COMPUTER DOES NOT EXIST. EXITING...";
        exit(1);
    }

    if (dmdOn == 1) {
        std::cout << "***DMD EMULATION ON***\n"; 
        std::cout << "\tEmulation will BEGIN at " << dmdObs << " solver interations, or s = " << (float)dmdObs*dt << "\n";
        std::cout << "\tEmulating in r = " << dmdTrunc << " truncated measurement subspace\n";
    }
        

}

void write_log(std::ofstream *out_file_imsrg, std::vector<state_type> data_log) {

    if((*out_file_imsrg).is_open()) {

        for (int i = 0; i < data_log.size(); i++) {
            state_type x = data_log[i];

            // for (int j = 0; j < x.size(); j++)
            //     out_file_imsrg->write((char *) &x[j], sizeof(double));

            // for (int j = 0; j < x.size(); j++)
            //     *out_file_imsrg << x[j];

            for (int j = 0; j < x.size()-1; j++)
                *out_file_imsrg << x[j] << ',';
            *out_file_imsrg << x[x.size()-1];
            *out_file_imsrg << std::fixed << std::setprecision(13) << "\n";
        }
    } else {
        std::cout << "Log file not created for some reason" << std::endl;
    }


    // std::ofstream out_file(file_path);
    // if (out_file.is_open()) {
    //     for (int i = 0; i < data.size(); i++) {
    //         state_type x = data[i];
            
    //         for (int j = 0; j < x.size(); j++)
    //             out_file << x[j] << ',';

    //         out_file << '\n';
    //     }

    //     out_file.close();
    // }

}

void write_rho(std::ofstream *out_file_rho, boost::numeric::ublas::vector<double> rho) {
    if((*out_file_rho).is_open()) {

        for (int i = 0; i < rho.size(); i++) {
            double x = rho[i];
            *out_file_rho << x;
            *out_file_rho << std::fixed << std::setprecision(13) << "\n";
        }
    } else {
        std::cout << "rho log file not created for some reason" << std::endl;
    }
}


int main(int argc, char **argv) {
    omp_set_num_threads(4);

    using namespace taco;
    using namespace boost::numeric::odeint;
    using namespace boost::numeric::ublas;
    adams_bashforth_moulton<5, state_type> abm_stepper;
    runge_kutta4<state_type> rk_stepper;
    std::string basis_path;
    int numWeights;
    //using namespace std;

    // int nholes = 4;
    // int nparticles = 4;
    int numStates, nHoles, nParticles, generator_id, reference_type, solver, BACKEND_id, dmdOn, dmdObs, dmdTrunc;
    double d, g, pb, t0, t1, dt, t1_integrate;

    // for IMSRG
    OccupationFactors occ;

    Backend *backend; 
    Backend_UBLAS cb_ublas;
    Backend_TACO cb_taco;

    Generator *generator;
    White white;
    WhiteAtan whiteAtan;
    Brillouin brillouin;

    PairingHamiltonian H;
    Flow_IMSRG2 flow;

    double E;
    vector<double> f, Gamma, W;

    System sys;
    state_type y0;
    std::vector<state_type> data_log;
    state_type test;

    vector<double> ref;
    

    // for Reference Ensemble
    matrix<int> basis;
    std::ifstream in_file;
    vector<double> weights;
    vector<double> rho1b, rho2b, rho3b, lambda2b, lambda3b;
    vector<double> currentVec;

    // for logging
    std::string log_dir;
    std::ofstream out_file_vac;
    std::ofstream out_file_imsrg;

    std::string rho_dir;
    std::ofstream out_file_rho1b;
    std::ofstream out_file_rho2b;

    // for DMD
    DMD dmd;

    process_user_input(argc, argv, 
                       numStates,
                       nHoles, nParticles,
                       d, g, pb,
                       t0, t1, dt,
                       generator_id, reference_type,
                       solver, BACKEND_id,
                       dmdOn, dmdObs, dmdTrunc);
    
    ref = vector<double>(numStates);
    if (reference_type == 0) {

        // build the single reference filling to Fermi surface        
        for (int i = 0; i < numStates; i++) {
            if (i < nHoles)
                ref[i] = 1.0;
            else
                ref[i] = 0.0;
        }
        
    } else {
        basis_path = "sd"+std::to_string(numStates)+".basis";

        std::cout << "Reading SD basis from " << basis_path << std::endl;

        //matrix<int> basis;
        basis = readBasisFromFile(basis_path);

        in_file = std::ifstream(basis_path);

        if (in_file.good()) {
            size_t pos = 0;
            std::string line;
            std::getline(in_file, line);

            std::string delim = " ";
            std::string token;

            pos = line.find(delim);
            token = line.substr(0, pos);
            line.erase(0, pos + delim.length());
            numWeights = std::stoi(line);
        } else {
            std::cout << "Failure to read basis file at " << basis_path << " exiting..." << std::endl;
            exit(1);
        }
        
        weights = vector<double>(numWeights);

        if (in_file.good()) {

            for (int i = 0; i < numWeights; i++) {
                size_t pos = 0;
                std::string line;
                std::getline(in_file, line);

                std::string delim = ",";
                std::string w;

                pos = line.find(delim);
                w = line.substr(0, pos);            
                //line.erase(0, pos + delim.length());
                weights[i] = std::stof(w);
            }

        } else {
            std::cout << "Failure to read basis file at " << basis_path << " exiting..." << std::endl;
            exit(1);
        }

        // weights[0] = sqrt(0.8);
        // weights[1] = sqrt(0.2);
        // for (int i = 0; i < numWeights; i++) {
        //     if (i == 0)
        //         weights[i] = (0.6);
        //     else if (i == 1)
        //         weights[i] = (0.1);
        //     else if (i == 2)
        //         weights[i] = (0.1);
        //     else if (i == 3)
        //         weights[i] = (0.1);
        //     else if (i == 4) 
        //         weights[i] = (0.1);                 
        //     else
        //         weights[i] = 0.0;
        // }
    
        std::cout << "Generating 1b/2b density matrices from SD basis..." << std::endl;

        // need to change arguments here; should be nholes, nparticles; make sure we can generate basis without half-filling restriction
        // rho1b = density_1b((int)numStates/2,(int)numStates/2, weights, 14, basis_path);
        // rho2b = density_2b((int)numStates/2,(int)numStates/2, weights, 14, basis_path);

        rho1b = density_1b(nHoles, nParticles, weights, 14, basis_path);
        rho2b = density_2b(nHoles, nParticles, weights, 14, basis_path);
        rho3b = density_3b(nHoles, nParticles, weights, 14, basis_path);
        lambda2b = density_2b_irr(nHoles+nParticles, rho1b, rho2b, 14);
        lambda3b = density_3b_irr(nHoles+nParticles, rho1b, lambda2b, rho3b, 14);

        //int numStates = nholes + nparticles;
        for (int i = 0; i < numWeights; i++) {
            currentVec = matrix_row<matrix<int>>(basis, i);
            vector_range<vector<double>> multiply_vec (currentVec, range(1,currentVec.size()));

            ref += weights[i]*multiply_vec;
        }
    }

    // ref[0] = 1.0;
    // ref[1] = 1.0;
    // ref[2] = 0.8;
    // ref[3] = 0.8;
    // ref[4] = 0.2;
    // ref[5] = 0.2;
    // ref[6] = 0.0;
    // ref[7] = 0.0;

    // int idx;
    // double val = 0.0;
    // for (int p = 0; p < numStates; p++) {

    //     idx = p*numStates +p;
    //     ref[p] = rho1b[idx];

    //     // val = (p < (int)numStates/2) ? 1.0 : 0.0;
    //     // //printf("%0.2f\n",val);
    //     // //ref.insert({p}, val);
    //     // ref[p] = val;
    // }

    std::cout << "Reference state = " << ref << std::endl;
    //ref.pack();

    occ = OccupationFactors(numStates, reference_type, ref);

    cb_ublas = Backend_UBLAS(numStates, &occ);
    cb_taco = Backend_TACO(numStates, &occ);

    switch (BACKEND_id) {
    case 0:
        backend = &cb_ublas;
        break;        
    case 1:
        backend = &cb_taco;
        break;
    default:
        backend = &cb_ublas;
    }

    white = White(numStates, ref);
    whiteAtan = WhiteAtan(numStates, ref);
    brillouin = Brillouin(numStates, ref, rho1b, lambda2b, lambda3b, &occ);

    switch (generator_id) {
    case 0:
        generator = &white;
        break;
    case 1:
        generator = &whiteAtan;
        break;
    case 2:
        generator = &brillouin; 
        break;
    default:
        generator = &white;
    }

    switch (reference_type) {
    case 0:
        H = PairingHamiltonian(nHoles, nParticles, ref, d,g,pb);
        break;
    case 1:
        H = PairingHamiltonian(numStates, rho1b, rho2b, d,g,pb);
        break;
    default:
        H = PairingHamiltonian(nHoles, nParticles, ref, d,g,pb);
    }

    flow = Flow_IMSRG2(occ, backend);

    E = H.get_E();
    f = H.get_f();
    Gamma = H.get_Gamma();
    //W;

    // Make the rho directory and write the rhos to files
    // Make the rho log directory
    rho_dir = "./rho/";
    const char* path_char = &rho_dir[0];
    int check = mkdir(path_char, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    if (check == -1) {
        std::cout << "Rho dir already exists at ./rho/" << std::endl;
    }

    boost::format rho_suffix = boost::format("%02d.out")%numStates;
    out_file_rho1b = std::ofstream(rho_dir+"rho1b_"+rho_suffix.str());
    out_file_rho2b = std::ofstream(rho_dir+"rho2b_"+rho_suffix.str());
    write_rho(&out_file_rho1b, rho1b);
    write_rho(&out_file_rho2b, rho2b);


    // Make the flow log directory and set up flow logging files
    log_dir = "./flow/";
    path_char = &log_dir[0];
    check = mkdir(path_char, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    if (check == -1) {
        std::cout << "Log dir already exists at ./flow/" << std::endl;
    }

    std::string log_path = log_dir+"H";
    if (reference_type == 0)
        log_path += "S";
    else if (reference_type == 1)
        log_path += "E";
    else 
        log_path += "S";

    boost::format fmt = boost::format("%02d-%f-%f-%f-%0.2f-%0.2f-%0.2f")%numStates%d%g%pb%t0%t1%dt;
    log_path += fmt.str();
    // log_path += std::to_string(numStates);
    // log_path += "-"+std::to_string(d);
    // log_path += "-"+std::to_string(g);
    // log_path += "-"+std::to_string(pb);
    // log_path += "-"+std::to_string(t0);
    // log_path += "-"+std::to_string(t1);
    // log_path += "-"+std::to_string(dt);
    log_path += ".log";
    std::cout << "Writing vacuum coefficient flow data log to " << log_path << std::endl;
    out_file_vac = std::ofstream(log_path); 

    log_path += ".imsrg";
    std::cout << "Writing IMSRG coefficient flow data log to " << log_path << std::endl;
    out_file_imsrg = std::ofstream(log_path); 

    out_file_vac << "#numStates," << numStates << std::endl;
    out_file_vac << "#d," << d << std::endl;
    out_file_vac << "#g," << g << std::endl;
    out_file_vac << "#pb," << pb << std::endl;
    out_file_vac << "#t0," << t0 << std::endl;
    out_file_vac << "#t1," << t1 << std::endl;
    out_file_vac << "#dt," << dt << std::endl;

    out_file_imsrg << "#numStates," << numStates << std::endl;
    out_file_imsrg << "#d," << d << std::endl;
    out_file_imsrg << "#g," << g << std::endl;
    out_file_imsrg << "#pb," << pb << std::endl;
    out_file_imsrg << "#t0," << t0 << std::endl;
    out_file_imsrg << "#t1," << t1 << std::endl;
    out_file_imsrg << "#dt," << dt << std::endl;


    sys = System(numStates, rho1b, lambda2b, lambda3b, E, f, Gamma, W, generator, &flow, &out_file_vac, &out_file_imsrg, reference_type);
    y0 = state_type(1+f.size()+Gamma.size());

    sys.system2vector(E, f, Gamma, y0);

    //size_t steps = integrate_n_steps(stepper, sys, y0, 0.0, 0.1, 500);
    //size_t steps = integrate_n_steps(stepper, sys, y0, 0.0, 0.1, 500, sys);

    if (dmdOn == 1)
        t1_integrate = (float)dmdObs*dt;
    else
        t1_integrate = t1;

    auto start = std::chrono::high_resolution_clock::now();

    size_t steps;
    switch (solver) {
    case 0:
        steps = integrate_adaptive(abm_stepper, sys, y0, t0, t1_integrate, dt, boost::ref(sys));
        break;
    case 1:
        steps = integrate_adaptive(rk_stepper, sys, y0, t0, t1_integrate, dt, boost::ref(sys));
        break;
    default:
        steps = integrate_adaptive(abm_stepper, sys, y0, t0, t1_integrate, dt, boost::ref(sys));
    }

    data_log = sys.getFlowData();

    if (dmdOn == 1) {
        double alpha = 1.0, beta = 0.0;
        MKL_INT incx = 1, incy = 1;

        int nRows = data_log[0].size();
        int nColumns = data_log.size();

        double* dmdData = (double*)malloc(nRows*nColumns*sizeof(double));
        double* ys = (double*)malloc(nRows*sizeof(double));  
        double* work_phi_eig = (double*)malloc(nRows*dmdTrunc*sizeof(double));

        state_type ys_to_log(nRows);

        double numItersTotal = (t1-t0)/dt;
        
        std::cout << "numItersTotal, " << numItersTotal << std::endl;

        for (int j = 0; j < nColumns; j++) {
            for (int i = 0; i < nRows; i++) {
                dmdData[i*nColumns + j] = ceil(data_log[j][i]*1e13)/1e13;
                work_phi_eig[i] = 0.0;
                ys[i] = 0.0;
            }
        }
       
        compute_dmd(&dmd, dmdData, nColumns, nRows, dmdTrunc);

        // set background
        dmd.eig[0] = 1.0;
        for (int i = 1; i < dmdTrunc; i++) {
            if (dmd.eig[i] > 1 || dmd.eig[i] < 0) {
                printf("%6.8f [%d], ", dmd.eig[i], i);
                dmd.eig[i] = 1.;
            }
        }
        printf(" ARE UNSTABLE EIGENVALUES; SETTING TO 1\n");

        printf("\tDMD eigenvalues: ");
        for (int i = 0; i < dmdTrunc; i++)
            printf("%6.4f,", dmd.eig[i]);
        printf("\n");

        printf("\tPhi values: ");
        for (int i = 0; i < dmdTrunc; i++)
            printf("%6.4f,", dmd.phi[i]);
        printf("\n");


        for (int p = steps+1; p < numItersTotal+1; p++) {

            for (int i = 0; i < nRows; i++) {
                for (int j = 0; j < dmdTrunc; j++) {
                    work_phi_eig[i*dmdTrunc + j] = dmd.phi[i*dmdTrunc + j]*pow(dmd.eig[j], p);
                }
            }
            cblas_dgemv(CblasRowMajor, CblasNoTrans, nRows, dmdTrunc, alpha, work_phi_eig, dmdTrunc, dmd.b, incx, beta, ys, incy);

            //state_type ys_to_log(ys);
            std::copy(ys,ys+nRows, ys_to_log.begin());

            data_log.push_back(ys_to_log);

            printf("%0.4f  %0.8f\n", p*dt, ys[0]);
        }


        free(ys);
        free(dmdData);
        free(dmd.phi);
        free(dmd.eig);
        free(dmd.b);
    }
        
    //std::cout << sys.getFlowData() << std::endl;

    auto stop = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    
    std::cout << "\nDone " << duration.count() << " microseconds" << std::endl;

    write_log(&out_file_imsrg, data_log);
    std::cout << data_log.size() << " lines from data log successfully written." << std::endl;

    out_file_vac << "#elapsed," << duration.count() << std::endl;
    out_file_imsrg << "#elapsed," << duration.count() << std::endl;

    out_file_vac.close();
    out_file_imsrg.close();    

    out_file_rho1b.close();
    out_file_rho2b.close();

    std::system(("gzip -f " + log_path).c_str());
    std::remove((log_path).c_str());
    
    //std::cout << sys.getFlowData() << std::endl;
    //write_log(log_path, sys.getFlowData());

    //delete observer;

    // for(size_t i = 0; i<=steps; i++)
    //     cout << times[i] << '\t' << x_vec[i][0] << '\n';

    // delete &H;
    // delete &occ;
    // delete &white;
}
