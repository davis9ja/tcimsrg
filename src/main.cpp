#include <stdio.h>
#include <ctime>
#include <boost/numeric/odeint.hpp>

#include "taco.h"
#include "pairinghamiltonian.hpp"
#include "occupation_factors.hpp"
#include "white.hpp"
#include "flow_imsrg2.hpp"
#include "system.hpp"
#include "BACKEND_ublas.cpp"

//#include "derivative.hpp"


struct push_back_system {
    std::vector<state_type> &m_states;
    std::vector<double> &m_times;
    System *sys;

    push_back_system( std::vector< state_type > &states , std::vector< double > &times , System *sys)
        : m_states( states ) , m_times( times ) { this->sys = sys; }

    void operator()( const state_type &x , double t )
    {
        m_states.push_back( x );
        m_times.push_back( t );
        std::cout << t << "\t" << x[0] << "\t"  << (*sys).getEta2bNorm() << std::endl;
    }
};


int main() {
    using namespace taco;
    using namespace boost::numeric::odeint;
    using namespace boost::numeric::ublas;
    //using namespace std;

    int nholes = 4;
    int nparticles = 4;
    double d = 1.0;
    double g = 0.5;
    double pb = 0.0;

    int numStates = nholes + nparticles;
    vector<double> ref(numStates);

    double val = 0.0;
    for (int p = 0; p < numStates; p++) {
        val = (p < nholes) ? 1.0 : 0.0;
        //printf("%0.2f\n",val);
        //ref.insert({p}, val);
        ref[p] = val;
    }
    //ref.pack();

    OccupationFactors occ(nholes, nparticles, ref);

    Backend *backend; 
    Backend_UBLAS cb_ublas(numStates, &occ);
    backend = &cb_ublas;

    White white(nholes, nparticles, ref);
    PairingHamiltonian H(nholes, nparticles, ref, d,g,pb);
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


    time_t ti, tf;

    ti = clock();
    vector<double> eta1b = white.compute_1b(f, Gamma, W);
    tf = clock();
    printf("calculate eta1b, %.4e\n", double(tf-ti)/CLOCKS_PER_SEC);

    ti = clock();
    vector<double> eta2b = white.compute_2b(f, Gamma, W);
    tf = clock();
    printf("calculate eta2b, %.4e\n", double(tf-ti)/CLOCKS_PER_SEC);

    ti = clock();
    double flow_0b = flow.flow_0b(f, Gamma, eta1b, eta2b);
    tf = clock();
    printf("calculate flow_0b, %.4e\n", double(tf-ti)/CLOCKS_PER_SEC);

    ti = clock();
    vector<double> flow_1b = flow.flow_1b(f, Gamma, eta1b, eta2b);
    tf = clock();
    printf("calculate flow_1b, %.4e\n", double(tf-ti)/CLOCKS_PER_SEC);

    ti = clock();
    vector<double> flow_2b = flow.flow_2b(f, Gamma, eta1b, eta2b);
    tf = clock();
    printf("calculate flow_2b, %.4e\n", double(tf-ti)/CLOCKS_PER_SEC);


    // for (int i = 0; i < numStates; i++)
    //     for (int j = 0; j < numStates; j++)
    //         std::cout << i << j << " " << f[i*numStates+j] << std::endl;
    

    //eta1b = white.compute_1b(flow_1b, flow_2b, W);


    printf("ETA1B ********************************\n");
    std::cout << eta1b << std::endl;

    printf("ETA2B ********************************\n");
    std::cout << eta2b << std::endl;

    printf("FLOW0B ********************************\n");
    std::cout << flow_0b << std::endl;
    
    printf("FLOW1B ********************************\n");
    std::cout << flow_1b << std::endl;

    printf("FLOW2B ********************************\n");
    std::cout << flow_2b << std::endl;

    //vector<double> ref = H.getReference();
    //double result = static_cast<const double*>(E.getStorage().getValues().getData())[0];

    std::cout << E  << std::endl;
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

    
    System sys(numStates, E, f, Gamma, W, &white, &flow);                                                     //
                                                                                                              //
    std::vector<state_type> x_vec;                                                                            //
    std::vector<double> times;                                                                                //
    printf("geting ready for integration\n");                                                                 //
                                                                                                              //
    state_type y0(1+f.size()+Gamma.size());
    sys.system2vector(E, f, Gamma, y0);

    printf("starting integration\n");                                                                         //
                                                                                                              //
    runge_kutta4<state_type> stepper;                                                                         //
    size_t steps = integrate_n_steps(stepper, sys, y0, 0.0, 0.1, 500, push_back_system(x_vec, times, &sys)); //
    

    // for(size_t i = 0; i<=steps; i++)
    //     cout << times[i] << '\t' << x_vec[i][0] << '\n';

    // delete &H;
    // delete &occ;
    // delete &white;
}
