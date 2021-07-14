#ifndef GENERATOR_HPP_
#define GENERATOR_HPP_

#include "generator.hpp"
using namespace taco;

class White: public Generator {

private:
    int n_holes, n_particles, numStates;
    Tensor<double> eta1b, eta2b, reference;

public:
    White(int n_holes, int n_particles, Tensor<double> reference);

    Tensor<double> compute_1b(Tensor<double> f, Tensor<double> Gamma, Tensor<double> W);
    Tensor<double> compute_2b(Tensor<double> f, Tensor<double> Gamma, Tensor<double> W);
    Tensor<double> compute_3b(Tensor<double> f, Tensor<double> Gamma, Tensor<double> W);

};


#endif