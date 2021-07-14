#include "taco.h"

using namespace taco;

class Generator {

public:
    Tensor<double> eta1b;
    Tensor<double> eta2b;
    Tensor<double> eta3b;
    
    int n_holes;
    int n_particles;

    virtual Tensor<double> compute_1b(Tensor<double> E, Tensor<double> f, Tensor<double> Gamma) = 0;
    virtual Tensor<double> compute_2b(Tensor<double> E, Tensor<double> f, Tensor<double> Gamma) = 0;
    virtual Tensor<double> compute_3b(Tensor<double> E, Tensor<double> f, Tensor<double> Gamma) = 0;

};








