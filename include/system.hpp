#ifndef SYSTEM_HPP_
#define SYSTEM_HPP_

//#include <boost/operators.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <math.h>
#include <omp.h>
#include <iostream>
#include <fstream>

#include "taco.h"
#include "white.hpp"
#include "flow_imsrg2.hpp"
#include "state_type.hpp"


class System {

private:
    double eta1b_norm = 1.0;
    double eta2b_norm = 1.0;

    double dE;
    boost::numeric::ublas::vector<double> df, dGamma;

    std::vector<state_type> data_log;
    std::ofstream *out_file;

    boost::numeric::ublas::vector<double> rho1b, rho2b;

    //SystemObserver *observer;

public:
    
    int numStates;
    double E;
    boost::numeric::ublas::vector<double> f, Gamma, W;
    //state_type sys_vec;
    White *white;
    Flow_IMSRG2 *flow;


    double getEta1bNorm() { return eta1b_norm; }
    double getEta2bNorm() { return eta2b_norm; }

    System(int numStates, boost::numeric::ublas::vector<double> &rho1b, boost::numeric::ublas::vector<double> &rho2b,
           double &E, boost::numeric::ublas::vector<double> &f, boost::numeric::ublas::vector<double> &Gamma, boost::numeric::ublas::vector<double> &W,
           White *white, 
           Flow_IMSRG2 *flow,
           std::ofstream *out_file
           );
    ~System();

    void system2vector(double &E, 
                       boost::numeric::ublas::vector<double> &f,
                       boost::numeric::ublas::vector<double> &Gamma,
                       state_type &x);

    void vector2system(const state_type &x, int fSize, 
                       double &E, 
                       boost::numeric::ublas::vector<double> &f,
                       boost::numeric::ublas::vector<double> &Gamma);

    // void reinitSystem(state_type x);
    void operator() (const state_type &x, state_type &dxdt, const double t);
    // void write_step();

    void operator()(const state_type &x , double t);

    std::vector<state_type> getFlowData() {
        return data_log;
    }
};


#endif
