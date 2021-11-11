#ifndef SYSTEM_OBSERVER_HPP_
#define SYSTEM_OBSERVER_HPP_

#include "state_type.hpp"

class SystemObserver {
    
private:
    static SystemObserver *instance;
    std::vector<state_type> states;
    std::vector<double> times;
    double eta2b_norm;

    //SystemObserver(std::vector<state_type> &states, std::vector< double > &times);    
    SystemObserver() {
        eta2b_norm = 0.0;        
    }

public:
    static SystemObserver *getInstance() {
        if(!instance)
            instance = new SystemObserver();
        return instance;
    }

    //SystemObserver(SystemObserver const&) = delete;
    //void operator=(SystemObserver const&) = delete;

    void operator()(const state_type &x , double t);
    void setEta2bNorm(double norm);
    double getEta2bNorm();
};
#endif
