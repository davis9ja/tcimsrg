#ifndef OCC_FACT_HPP_
#define OCC_FACT_HPP_

#include "taco.h"

using namespace taco;
using namespace std;

class OccupationFactors {


private:
    
    int n_holes;
    int n_particles;
    Tensor<double> ref;

    void writeA();
    void writeB();
    void writeC();
    void writeD();

    string path = "occ_storage/";

public:

    OccupationFactors(int n_holes, int n_particles, Tensor<double> ref);
    //~OccupationFactors();

};


#endif
