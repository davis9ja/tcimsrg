#ifndef SYSTEM_OBSERVER_HPP_
#define SYSTEM_OBSERVER_HPP_

#include "state_type.hpp"

class SystemObserver {

private:

    std::vector<state_type> states;
    std::vector<double> times;
    double eta2b_norm = 0.0;
    

public:

    SystemObserver(std::vector<state_type> &states, std::vector< double > &times);

    void operator()(const state_type &x , double t);
    void setEta2bNorm(double norm);
};
#endif
