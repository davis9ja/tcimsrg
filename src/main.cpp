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

#include "white.cpp"
#include "white_atan.cpp"
#include "BACKEND_ublas.cpp"
#include "BACKEND_taco.cpp"
//#include "derivative.hpp"

void process_user_input(int argc, char **argv, 
                        int &numStates, 
                        int &nHoles, int &nParticles,
                        double &d, double &g, double &pb, 
                        double &t0, double &t1, double &dt,
                        int &generator_id,
                        int &reference_type,
                        int &solver, int &BACKEND_id) {

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
            else
                std::cout << param_str << " not recognized as input parameter.\n" << std::endl;
        }
    }

    if (reference_type == 1)
        printf("Parameters:\n\t# s.p. states = %d\n\td=%0.3f\n\tg=%0.3f\n\tpb=%0.3f\n", numStates, d, g, pb);
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
        

}

void write_log(std::string file_path, std::vector<state_type> data) {

    std::ofstream out_file(file_path);
    if (out_file.is_open()) {
        for (int i = 0; i < data.size(); i++) {
            state_type x = data[i];
            
            for (int j = 0; j < x.size(); j++)
                out_file << x[j] << ',';

            out_file << '\n';
        }

        out_file.close();
    }

}

int main(int argc, char **argv) {
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
    int numStates, nHoles, nParticles, generator_id, reference_type, solver, BACKEND_id;
    double d, g, pb, t0, t1, dt;

    // for IMSRG
    OccupationFactors occ;

    Backend *backend; 
    Backend_UBLAS cb_ublas;
    Backend_TACO cb_taco;

    Generator *generator;
    White white;
    WhiteAtan whiteAtan;

    PairingHamiltonian H;
    Flow_IMSRG2 flow;

    double E;
    vector<double> f, Gamma, W;

    System sys;
    state_type y0;

    vector<double> ref;
    

    // for Reference Ensemble
    matrix<int> basis;
    std::ifstream in_file;
    vector<double> weights;
    vector<double> rho1b, rho2b;
    vector<double> currentVec;

    // for logging
    std::string log_dir;
    std::ofstream out_file_vac;
    std::ofstream out_file_imsrg;

    process_user_input(argc, argv, 
                       numStates,
                       nHoles, nParticles,
                       d, g, pb,
                       t0, t1, dt,
                       generator_id, reference_type,
                       solver, BACKEND_id);
    
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

        }
        
        weights = vector<double>(numWeights);
        // weights[0] = sqrt(0.8);
        // weights[1] = sqrt(0.2);
        for (int i = 0; i < numWeights; i++) {
            if (i == 0)
                weights[i] = (0.8);
            else if (i == 1)
                weights[i] = (0.2);
            else
                weights[i] = 0.0;
        }
    
        std::cout << "Generating 1b/2b density matrices from SD basis..." << std::endl;

        rho1b = density_1b((int)numStates/2,(int)numStates/2, weights, 14, basis_path);
        rho2b = density_2b((int)numStates/2,(int)numStates/2, weights, 14, basis_path);

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

    switch (generator_id) {
    case 0:
        generator = &white;
        break;
    case 1:
        generator = &whiteAtan;
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

    log_dir = "./flow/";
    const char* path_char = &log_dir[0];
    int check = mkdir(path_char, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    if (check == -1) {
        std::cout << "Log dir already exists at ./flow/" << std::endl;
    }

    std::string log_path = "./flow/H";
    if (reference_type == 0)
        log_path += "S";
    else if (reference_type == 1)
        log_path += "E";
    else 
        log_path += "S";

    boost::format fmt = boost::format("%02d-%0.2f-%0.2f-%0.2f-%0.2f-%0.2f-%0.2f")%numStates%d%g%pb%t0%t1%dt;
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


    sys = System(numStates, rho1b, rho2b, E, f, Gamma, W, generator, &flow, &out_file_vac, &out_file_imsrg, reference_type);
    y0 = state_type(1+f.size()+Gamma.size());

    sys.system2vector(E, f, Gamma, y0);

    //size_t steps = integrate_n_steps(stepper, sys, y0, 0.0, 0.1, 500);
    //size_t steps = integrate_n_steps(stepper, sys, y0, 0.0, 0.1, 500, sys);

    auto start = std::chrono::high_resolution_clock::now();

    size_t steps;
    switch (solver) {
    case 0:
        steps = integrate_adaptive(abm_stepper, sys, y0, t0, t1, dt, sys);
        break;
    case 1:
        steps = integrate_adaptive(rk_stepper, sys, y0, t0, t1, dt, sys);
        break;
    default:
        steps = integrate_adaptive(abm_stepper, sys, y0, t0, t1, dt, sys);        
    }

    auto stop = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    
    std::cout << "\nDone " << duration.count() << " microseconds" << std::endl;

    out_file_vac << "#elapsed," << duration.count() << std::endl;
    out_file_imsrg << "#elapsed," << duration.count() << std::endl;

    out_file_vac.close();
    out_file_imsrg.close();    

    //std::cout << sys.getFlowData() << std::endl;
    //write_log(log_path, sys.getFlowData());

    //delete observer;

    // for(size_t i = 0; i<=steps; i++)
    //     cout << times[i] << '\t' << x_vec[i][0] << '\n';

    // delete &H;
    // delete &occ;
    // delete &white;
}
