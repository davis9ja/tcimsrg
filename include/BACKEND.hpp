#ifndef BACKEND_HPP_
#define BACKEND_HPP_

#include <boost/numeric/ublas/vector.hpp>
#include "occupation_factors.hpp"
#include "imsrg_utils.hpp"

class Backend {


public:

    OccupationFactors *occFact;
    
    
    virtual double flow_0b(boost::numeric::ublas::vector<double> &occA,
                           boost::numeric::ublas::vector<double> &occD,
                           boost::numeric::ublas::vector<double> &f, 
                           boost::numeric::ublas::vector<double> &Gamma, 
                           boost::numeric::ublas::vector<double> &eta1b, 
                           boost::numeric::ublas::vector<double> &eta2b) = 0;

    virtual boost::numeric::ublas::vector<double> flow_1b(boost::numeric::ublas::vector<double> &occA,
                                                          boost::numeric::ublas::vector<double> &occC,                                                          
                                                          boost::numeric::ublas::vector<double> &f, 
                                                          boost::numeric::ublas::vector<double> &Gamma, 
                                                          boost::numeric::ublas::vector<double> &eta1b, 
                                                          boost::numeric::ublas::vector<double> &eta2b) = 0;

    virtual boost::numeric::ublas::vector<double> flow_2b(boost::numeric::ublas::vector<double> &occA,
                                                          boost::numeric::ublas::vector<double> &occB,
                                                          boost::numeric::ublas::vector<double> &f, 
                                                          boost::numeric::ublas::vector<double> &Gamma, 
                                                          boost::numeric::ublas::vector<double> &eta1b, 
                                                          boost::numeric::ublas::vector<double> &eta2b) = 0;


};


#endif
