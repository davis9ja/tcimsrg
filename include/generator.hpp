#ifndef GENERATOR_HPP_
#define GENERATOR_HPP_

//#include "taco.h"
#include <boost/numeric/ublas/vector.hpp>
//using namespace taco;

class Generator {

public:
    boost::numeric::ublas::vector<double> eta1b;
    boost::numeric::ublas::vector<double> eta2b;
    boost::numeric::ublas::vector<double> eta3b;
    
    int n_holes;
    int n_particles;

    virtual boost::numeric::ublas::vector<double> compute_1b(boost::numeric::ublas::vector<double> E, boost::numeric::ublas::vector<double> f, boost::numeric::ublas::vector<double> Gamma) = 0;
    virtual boost::numeric::ublas::vector<double> compute_2b(boost::numeric::ublas::vector<double> E, boost::numeric::ublas::vector<double> f, boost::numeric::ublas::vector<double> Gamma) = 0;
    virtual boost::numeric::ublas::vector<double> compute_3b(boost::numeric::ublas::vector<double> E, boost::numeric::ublas::vector<double> f, boost::numeric::ublas::vector<double> Gamma) = 0;

};

#endif
