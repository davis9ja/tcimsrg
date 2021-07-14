#include <stdio.h>
#include <ctime>

#include "taco.h"
#include "pairinghamiltonian.hpp"
#include "occupation_factors.hpp"
#include "white.hpp"

using namespace taco;

int main() {
    int nholes = 4;
    int nparticles = 4;
    double d = 1.0;
    double g = 0.5;
    double pb = 0.0;

    int numStates = nholes + nparticles;
    Tensor<double> ref("reference", {numStates}, Format(Dense));

    double val = 0.0;
    for (int p = 0; p < numStates; p++) {
        val = (p < nholes) ? 1.0 : 0.0;
        //printf("%0.2f\n",val);
        ref.insert({p}, val);
    }
    ref.pack();

    
    White* white = new White(nholes, nparticles, ref);
    OccupationFactors* occ = new OccupationFactors(nholes, nparticles, ref);
    PairingHamiltonian* H = new PairingHamiltonian(nholes, nparticles, ref, d,g,pb);

    double E = H->get_E();//taco::Tensor<double> E = H->get_E();
    taco::Tensor<double> f = H->get_f();
    taco::Tensor<double> Gamma = H->get_Gamma();
    taco::Tensor<double> W;

    std::time_t ti, tf;

    ti = std::clock();
    Tensor<double> eta1b = white->compute_1b(f, Gamma, W);
    tf = std::clock();
    printf("%.4e\n", double(tf-ti)/CLOCKS_PER_SEC);

    ti = std::clock();
    Tensor<double> eta2b = white->compute_2b(f, Gamma, W);
    tf = std::clock();
    printf("%.4e\n", double(tf-ti)/CLOCKS_PER_SEC);

    std::cout << eta1b << std::endl;
    std::cout << eta2b << std::endl;
    //taco::Tensor<double> ref = H->getReference();

    //std::cout << E << std::endl;
    //std::cout << ref << std::endl;
    
    delete H;
    delete occ;
    delete white;
}
