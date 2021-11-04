#include <stdio.h>

#include <boost/numeric/ublas/io.hpp>


#include "flow_imsrg2.hpp"

using namespace taco;
using namespace boost::numeric::ublas;

Flow_IMSRG2::Flow_IMSRG2(OccupationFactors occFact, Backend *backend) {


    std::string factor_path = occFact.getPath();
    this->numStates = occFact.getNumStates();

    occFact.readOccTensors(factor_path, 
                           occA_a, occA_b,
                           occB_a, occB_b,
                           occC_a, occC_b, occC_c,
                           occD_a, occD_b, occD_c, occD_d
                           );

    this->backend = backend;
    
}

double Flow_IMSRG2::flow_0b(vector<double> &f, vector<double> &Gamma, 
                            vector<double> &eta1b, vector<double> &eta2b) {

    // Return
    double dE;

    // TRANSFORM OCCUPATION TENSORS
    std::vector<int> shape_r2 = {numStates,numStates};
    std::vector<int> shape_r4 = {numStates,numStates,numStates,numStates};
    Format format_r2({Dense,Dense});
    Format format_r4({Dense,Dense,Dense,Dense});

    Tensor<double> occA_t(shape_r2, format_r2);
    Tensor<double> occD_t(shape_r4, format_r4);
    vector<double> occA(numStates*numStates);
    vector<double> occD(numStates*numStates*numStates*numStates);

    IndexVar a,b,c,d,p,q;
    occA_t(a,b) = occA_a(a,p)*occA_b(p,b);    
    occD_t(a,b,c,d) = occD_a(a,p)*occD_b(p,b)*occD_c(c,q)*occD_d(q,d);

    occA_t.evaluate();
    occD_t.evaluate();

    tensor2vector(occA_t, occA);
    tensor2vector(occD_t, occD);

    // RUN FLOW ON BACKEND
    dE = backend->flow_0b(occA, occD, f, Gamma, eta1b, eta2b);

    return dE;
}

vector<double> Flow_IMSRG2::flow_1b(vector<double> &f, vector<double> &Gamma,
                                    vector<double> &eta1b, vector<double> &eta2b) {
    // Return
    vector<double> df;

    // TRANSFORM OCCUPATION TENSORS
    std::vector<int> shape_r2 = {numStates,numStates};
    std::vector<int> shape_r3 = {numStates,numStates,numStates};
    Format format_r2({Dense,Dense});
    Format format_r3({Dense,Dense,Dense});

    Tensor<double> occA_t(shape_r2, format_r2);
    Tensor<double> occC_t(shape_r3, format_r3);
    vector<double> occA(numStates*numStates);
    vector<double> occC(numStates*numStates*numStates);

    IndexVar a,b,c,p,q;
    occA_t(a,b) = occA_a(a,p)*occA_b(p,b);    
    occC_t(a,b,c) = occC_a(a,p)*occC_b(p,b,q)*occC_c(q,c);

    occA_t.evaluate();
    occC_t.evaluate();

    tensor2vector(occA_t, occA);
    tensor2vector(occC_t, occC);

    // RUN FLOW ON BACKEND
    df = backend->flow_1b(occA, occC, f, Gamma, eta1b, eta2b);    

    return df;
}

vector<double> Flow_IMSRG2::flow_2b(vector<double> &f, vector<double> &Gamma, 
                                    vector<double> &eta1b, vector<double> &eta2b) {

    // Return
    vector<double> dGamma;

    // TRANSFORM OCCUPATION TENSORS
    std::vector<int> shape_r2 = {numStates,numStates};
    Format format_r2({Dense,Dense});

    Tensor<double> occA_t(shape_r2, format_r2);
    Tensor<double> occB_t(shape_r2, format_r2);
    vector<double> occA(numStates*numStates);
    vector<double> occB(numStates*numStates);

    IndexVar a,b,c,p,q;
    occA_t(a,b) = occA_a(a,p)*occA_b(p,b);    
    occB_t(a,b) = occB_a(a,p)*occB_b(p,b);    

    occA_t.evaluate();
    occB_t.evaluate();

    tensor2vector(occA_t, occA);
    tensor2vector(occB_t, occB);

    dGamma = backend->flow_2b(occA, occB, f, Gamma, eta1b, eta2b);

    return dGamma;
}
