#include <stdio.h>

#include <boost/numeric/ublas/io.hpp>


#include "flow_imsrg2.hpp"

using namespace taco;
using namespace boost::numeric::ublas;

Flow_IMSRG2::Flow_IMSRG2(){}

Flow_IMSRG2::Flow_IMSRG2(OccupationFactors occFact, Backend *backend) {


    //std::string factor_path = occFact.getPath();
    //this->numStates = occFact.getNumStates();

    // occFact.readOccTensors(factor_path, 
    //                        occA_a, occA_b,
    //                        occB_a, occB_b,
    //                        occC_a, occC_b, occC_c,
    //                        occD_a, occD_b, occD_c, occD_d
    //                        );

    this->backend = backend;
    
}

// SINGLE REFERENCE EQUATIONS
double Flow_IMSRG2::flow_0b(vector<double> &f, vector<double> &Gamma, 
                            vector<double> &eta1b, vector<double> &eta2b) {
    // Return
    double dE;

    // RUN FLOW ON BACKEND
    dE = backend->flow_0b(f, Gamma, eta1b, eta2b);

    return dE;
}

vector<double> Flow_IMSRG2::flow_1b(vector<double> &f, vector<double> &Gamma,
                                    vector<double> &eta1b, vector<double> &eta2b) {
    // Return
    vector<double> df;

    // RUN FLOW ON BACKEND
    df = backend->flow_1b(f, Gamma, eta1b, eta2b);    

    return df;
}

vector<double> Flow_IMSRG2::flow_2b(vector<double> &f, vector<double> &Gamma, 
                                    vector<double> &eta1b, vector<double> &eta2b) {

    // Return
    vector<double> dGamma;

    // RUN FLOW ON BACKEND
    dGamma = backend->flow_2b(f, Gamma, eta1b, eta2b);

    return dGamma;
}

// MULTI-REFERENCE EQUATIONS INCLUDING IRREDUCIBLE 1B/2B DENSITIES
double Flow_IMSRG2::flow_0b(vector<double> &f, vector<double> &Gamma, 
                            vector<double> &eta1b, vector<double> &eta2b,
                            vector<double> &rho1b, vector<double> &rho2b,
                            vector<double> &dGamma) {

    // Return
    double dE;

    // RUN FLOW ON BACKEND
    dE = backend->flow_0b(f, Gamma, eta1b, eta2b, rho1b, rho2b, dGamma);

    return dE;
}

vector<double> Flow_IMSRG2::flow_1b(vector<double> &f, vector<double> &Gamma,
                                    vector<double> &eta1b, vector<double> &eta2b,
                                    vector<double> &rho1b, vector<double> &rho2b) {
    // Return
    vector<double> df;

    // RUN FLOW ON BACKEND
    df = backend->flow_1b(f, Gamma, eta1b, eta2b, rho1b, rho2b);  

    return df;
}

vector<double> Flow_IMSRG2::flow_2b(vector<double> &f, vector<double> &Gamma, 
                                    vector<double> &eta1b, vector<double> &eta2b,
                                    vector<double> &rho1b, vector<double> &rho2b) {

    // Return
    vector<double> dGamma;

    // RUN FLOW ON BACKEND
    dGamma = backend->flow_2b(f, Gamma, eta1b, eta2b, rho1b, rho2b);

    return dGamma;
}

