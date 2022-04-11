#ifndef GENERATOR_HPP_
#define GENERATOR_HPP_

//#include "taco.h"
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
//using namespace taco;

class Generator {

public:
    boost::numeric::ublas::vector<double> eta1b;
    boost::numeric::ublas::vector<double> eta2b;
    boost::numeric::ublas::vector<double> eta3b;
    
    int numStates;

    virtual boost::numeric::ublas::vector<double> compute_1b(boost::numeric::ublas::vector<double> &f, boost::numeric::ublas::vector<double> &Gamma, boost::numeric::ublas::vector<double> &W) = 0;
    virtual boost::numeric::ublas::vector<double> compute_2b(boost::numeric::ublas::vector<double> &f, boost::numeric::ublas::vector<double> &Gamma, boost::numeric::ublas::vector<double> &W) = 0;
    virtual boost::numeric::ublas::vector<double> compute_3b(boost::numeric::ublas::vector<double> &f, boost::numeric::ublas::vector<double> &Gamma, boost::numeric::ublas::vector<double> &W) = 0;

};

#endif
