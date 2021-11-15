#ifndef WHITE_HPP_
#define WHITE_HPP_

#include "taco.h"
#include "generator.hpp"
//using namespace taco;

class White: public Generator {

private:
    int numStates;
    boost::numeric::ublas::vector<double> reference;
    //taco::Format f1b, f2b, f3b;

public:
    White(int numStates, boost::numeric::ublas::vector<double> &reference);

    boost::numeric::ublas::vector<double> compute_1b(boost::numeric::ublas::vector<double> &f, boost::numeric::ublas::vector<double> &Gamma, boost::numeric::ublas::vector<double> &W);
    boost::numeric::ublas::vector<double> compute_2b(boost::numeric::ublas::vector<double> &f, boost::numeric::ublas::vector<double> &Gamma, boost::numeric::ublas::vector<double> &W);
    boost::numeric::ublas::vector<double> compute_3b(boost::numeric::ublas::vector<double> &f, boost::numeric::ublas::vector<double> &Gamma, boost::numeric::ublas::vector<double> &W);

};


#endif
