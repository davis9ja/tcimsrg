#ifndef FL_IMSRG2_HPP_
#define FL_IMSRG2_HPP_

#include <math.h>
#include "occupation_factors.hpp"
#include "BACKEND.hpp"
#include "imsrg_utils.hpp"

//using namespace taco;
//using namespace std;

class Flow_IMSRG2 {
    
private:
    
    taco::Tensor<double> occA_a, occA_b;
    taco::Tensor<double> occB_a, occB_b;
    taco::Tensor<double> occC_a, occC_b, occC_c;
    taco::Tensor<double> occD_a, occD_b, occD_c, occD_d;
    
    int numStates;

    Backend *backend;

public:
    Flow_IMSRG2();
    Flow_IMSRG2(OccupationFactors occFact, Backend *backend);

    double flow_0b(boost::numeric::ublas::vector<double> &f, 
                   boost::numeric::ublas::vector<double> &Gamma, 
                   boost::numeric::ublas::vector<double> &eta1b, 
                   boost::numeric::ublas::vector<double> &eta2b);
    boost::numeric::ublas::vector<double> flow_1b(boost::numeric::ublas::vector<double> &f, 
                                                  boost::numeric::ublas::vector<double> &Gamma, 
                                                  boost::numeric::ublas::vector<double> &eta1b, 
                                                  boost::numeric::ublas::vector<double> &eta2b);
    boost::numeric::ublas::vector<double> flow_2b(boost::numeric::ublas::vector<double> &f, 
                                                  boost::numeric::ublas::vector<double> &Gamma, 
                                                  boost::numeric::ublas::vector<double> &eta1b, 
                                                  boost::numeric::ublas::vector<double> &eta2b);
    
};

#endif
