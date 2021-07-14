#ifndef PHAM_HPP_
#define PHAM_HPP_

#include "taco.h"
#include <string>
using namespace taco;

class PairingHamiltonian {

private:
    int delta2B(int p, int q, int r, int s);
    int deltaPB(int p, int q, int r, int s);

    Tensor<double> H1B; //matrix representation of 1B operator
    Tensor<double> H2B; //matrix representation of 2B operator

    double E;//Tensor<double> E;     //zero body piece of normal ordered H
    Tensor<double> f;     //one body piece of normal ordered H
    Tensor<double> Gamma; //two body piece of normal ordered H

    double d;
    double g;
    double pb;

    int* spBasis;
    int* holes;
    int* particles;
public:

    int n_holes;
    int n_particles;
    //Tensor<double> ref; //= std::make_unique<int[]>(n_holes+n_particles);

    PairingHamiltonian(int n_holes, int n_particles, Tensor<double> ref, double d, double g, double pb);
    ~PairingHamiltonian();

    double getSpacing() { return d; }
    double getStrength() { return g; }
    double getPBStrength() { return pb; }
    //Tensor<double> getReference() { return ref; }
    
    int* get_basisIdx() { return spBasis; }
    int* get_holesIdx() { return holes; }
    int* get_particlesIdx() { return particles; }

    Tensor<double> get_1B() { return H1B; }
    Tensor<double> get_2B() { return H2B; }
    double get_E() {return E; }//Tensor<double> get_E() { return E; }
    Tensor<double> get_f() { return f; }
    Tensor<double> get_Gamma() { return Gamma; }

    //string get_1B_str() {};
};

#endif
