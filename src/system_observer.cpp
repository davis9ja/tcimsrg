#include "system_observer.hpp"


SystemObserver::SystemObserver( std::vector< state_type > &states , std::vector< double > &times){ 
    this->states = states;
    this->times = times;
}

void SystemObserver::operator()( const state_type &x , double t )
{
    states.push_back( x );
    times.push_back( t );
    std::cout << t << "\t" << x[0] << "\t"  << eta2b_norm << std::endl;
}

void SystemObserver::setEta2bNorm(double norm) {
    eta2b_norm = norm;
}
