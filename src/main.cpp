#include <stdio.h>
#include <stdlib.h>
#include <ctime>
#include <boost/numeric/odeint.hpp>
#include <chrono>

#include "taco.h"
#include "pairinghamiltonian.hpp"
#include "occupation_factors.hpp"
#include "white.hpp"
#include "flow_imsrg2.hpp"
#include "system.hpp"
#include "BACKEND_ublas.cpp"
#include "BACKEND_taco.cpp"
//#include "derivative.hpp"

void process_user_input(int argc, char **argv, 
                        int &numStates, 
                        double &d, double &g, double &pb, 
                        double &t0, double &t1, double &dt,
                        int &solver, int &BACKEND_id) {

    // SET DEFAULT PARAMS
    numStates = 8;
    d = 1.0;
    g = 0.5;
    pb = 0.0;
    t0 = 0.0;
    t1 = 5.0;
    dt = 0.1;
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
            else if (param_str == "solver")
                solver = std::stof(param_val);
            else if (param_str == "BACKEND_id")
                BACKEND_id = std::stof(param_val);
            else
                std::cout << param_str << " not recognized as input parameter.\n" << std::endl;
        }
    }


    printf("Parameters:\n\t# s.p. states = %d\n\td=%0.3f\n\tg=%0.3f\n\tpb=%0.3f\n", numStates, d, g, pb);

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

int main(int argc, char **argv) {
    using namespace taco;
    using namespace boost::numeric::odeint;
    using namespace boost::numeric::ublas;
    adams_bashforth_moulton<5, state_type> abm_stepper;
    runge_kutta4<state_type> rk_stepper;

    //using namespace std;

    // int nholes = 4;
    // int nparticles = 4;
    int numStates, solver, BACKEND_id;
    double d, g, pb, t0, t1, dt;

    process_user_input(argc, argv, 
                       numStates,
                       d, g, pb,
                       t0, t1, dt,
                       solver, BACKEND_id);

    //int numStates = nholes + nparticles;
    vector<double> ref(numStates);
    
    double val = 0.0;
    for (int p = 0; p < numStates; p++) {
        val = (p < (int)numStates/2) ? 1.0 : 0.0;
        //printf("%0.2f\n",val);
        //ref.insert({p}, val);
        ref[p] = val;
    }
    //ref.pack();

    OccupationFactors occ(numStates, ref);

    Backend *backend; 
    Backend_UBLAS cb_ublas(numStates, &occ);
    Backend_TACO cb_taco(numStates, &occ);

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

    White white(numStates, ref);
    PairingHamiltonian H(numStates, ref, d,g,pb);
    Flow_IMSRG2 flow(occ, backend);

    double E = H.get_E();
    vector<double> f = H.get_f();
    vector<double> Gamma = H.get_Gamma();
    vector<double> W;

    // std::cout << "ref, " << ref << std::endl;
    // std::cout << "E, " << E << std::endl;
    // std::cout << "F, " << f << std::endl;
    // std::cout << "Gamma, " << Gamma << std::endl;

    // for (int a = 0; a < numStates; a++) {
    //     for (int b = 0; b < numStates; b++) {
    //         for (int c = 0; c < numStates; c++) {
    //             for (int d = 0; d < numStates; d++) {
    //                 double val = Gamma(a,b,c,d);
    //                 if (val != 0.0)
    //                     printf("%d%d%d%d, %0.8f\n", a,b,c,d,val);
    //             }
    //         }
    //     }
    // }


    // time_t ti, tf;

    // ti = clock();
    // vector<double> eta1b = white.compute_1b(f, Gamma, W);
    // tf = clock();
    // printf("calculate eta1b, %.4e\n", double(tf-ti)/CLOCKS_PER_SEC);

    // ti = clock();
    // vector<double> eta2b = white.compute_2b(f, Gamma, W);
    // tf = clock();
    // printf("calculate eta2b, %.4e\n", double(tf-ti)/CLOCKS_PER_SEC);

    // ti = clock();
    // double flow_0b = flow.flow_0b(f, Gamma, eta1b, eta2b);
    // tf = clock();
    // printf("calculate flow_0b, %.4e\n", double(tf-ti)/CLOCKS_PER_SEC);

    // ti = clock();
    // vector<double> flow_1b = flow.flow_1b(f, Gamma, eta1b, eta2b);
    // tf = clock();
    // printf("calculate flow_1b, %.4e\n", double(tf-ti)/CLOCKS_PER_SEC);

    // ti = clock();
    // vector<double> flow_2b = flow.flow_2b(f, Gamma, eta1b, eta2b);
    // tf = clock();
    // printf("calculate flow_2b, %.4e\n", double(tf-ti)/CLOCKS_PER_SEC);


    // // for (int i = 0; i < numStates; i++)
    // //     for (int j = 0; j < numStates; j++)
    // //         std::cout << i << j << " " << f[i*numStates+j] << std::endl;
    

    // //eta1b = white.compute_1b(flow_1b, flow_2b, W);


    // printf("ETA1B ********************************\n");
    // std::cout << eta1b << std::endl;

    // printf("ETA2B ********************************\n");
    // std::cout << eta2b << std::endl;

    // printf("FLOW0B ********************************\n");
    // std::cout << flow_0b << std::endl;
    
    // printf("FLOW1B ********************************\n");
    // std::cout << flow_1b << std::endl;

    // printf("FLOW2B ********************************\n");
    // std::cout << flow_2b << std::endl;

    // //vector<double> ref = H.getReference();
    // //double result = static_cast<const double*>(E.getStorage().getValues().getData())[0];

    // std::cout << E  << std::endl;
    //std::cout << ref << std::endl;

    // state_type c;
    // c.E = E;
    // c.f = f;
    // c.Gamma = Gamma;
    // c.W = W;
    // c.generator = white;
    // c.flow = flow;

    // container dcdt;
    // double t = 1.0;
    // derivative(c, dcdt, t);
    
    // vector<double> result;
    // result() = c.E() + dcdt.E();
    // cout << result << endl;

    // std::vector<state_type> x_vec;
    // std::vector<double> times;

    //static SystemObserver *observer; //= observer.getInstance(); //new SystemObserver(x_vec, times);
    //SystemObserver *observer = observer->getInstance();

    System sys(numStates, E, f, Gamma, W, &white, &flow );
    state_type y0(1+f.size()+Gamma.size());

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
    

    //delete observer;

    // for(size_t i = 0; i<=steps; i++)
    //     cout << times[i] << '\t' << x_vec[i][0] << '\n';

    // delete &H;
    // delete &occ;
    // delete &white;
}
