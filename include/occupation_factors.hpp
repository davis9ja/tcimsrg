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

    void readOccTensors(std::string factor_path,
                    taco::Tensor<double> &occA_a, taco::Tensor<double> &occA_b, 
                    taco::Tensor<double> &occB_a, taco::Tensor<double> &occB_b,
                    taco::Tensor<double> &occC_a, taco::Tensor<double> &occC_b, 
                    taco::Tensor<double> &occC_c,
                    taco::Tensor<double> &occD_a, taco::Tensor<double> &occD_b,
                    taco::Tensor<double> &occD_c, taco::Tensor<double> &occD_d);

};


#endif
