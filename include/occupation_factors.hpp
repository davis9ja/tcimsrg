#ifndef OCC_FACT_HPP_
#define OCC_FACT_HPP_

#include <boost/numeric/ublas/vector.hpp>
#include "taco.h"

class OccupationFactors {


private:
    
    int n_holes;
    int n_particles;
    boost::numeric::ublas::vector<double> ref;

    void writeA();
    void writeB();
    void writeC();
    void writeD();

    std::string path = "occ_storage/";

public:
    
    OccupationFactors(int n_holes, int n_particles, boost::numeric::ublas::vector<double> ref);
    std::string getPath();
    //~OccupationFactors();

    int getNumStates() { return n_holes+n_particles; }
};


#endif
