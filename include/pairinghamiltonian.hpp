#ifndef PHAM_HPP_
#define PHAM_HPP_

#include "taco.h"
#include <string>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

//using namespace taco;

class PairingHamiltonian {

private:
    int delta2B(int p, int q, int r, int s);
    int deltaPB(int p, int q, int r, int s);

    boost::numeric::ublas::vector<double> H1B; //matrix representation of 1B operator
    boost::numeric::ublas::vector<double> H2B; //matrix representation of 2B operator

    double E;     //zero body piece of normal ordered H
    boost::numeric::ublas::vector<double> f;     //one body piece of normal ordered H
    boost::numeric::ublas::vector<double> Gamma; //two body piece of normal ordered H

    double d;
    double g;
    double pb;

    // int* spBasis;
    // int* holes;
    // int* particles;
public:

    // int n_holes;
    // int n_particles;
    //boost::numeric::ublas::vector<double> ref; //= std::make_unique<int[]>(n_holes+n_particles);

    int numStates;

    PairingHamiltonian(int numStates, boost::numeric::ublas::vector<double> ref, double d, double g, double pb);
    //~PairingHamiltonian();

    double getSpacing() { return d; }
    double getStrength() { return g; }
    double getPBStrength() { return pb; }
    //boost::numeric::ublas::vector<double> getReference() { return ref; }
    
    // int* get_basisIdx() { return spBasis; }
    // int* get_holesIdx() { return holes; }
    // int* get_particlesIdx() { return particles; }

    boost::numeric::ublas::vector<double> get_1B() { return H1B; }
    boost::numeric::ublas::vector<double> get_2B() { return H2B; }
    double get_E() { return E; }
    boost::numeric::ublas::vector<double> get_f() { return f; }
    boost::numeric::ublas::vector<double> get_Gamma() { return Gamma; }

    //string get_1B_str() {};
};

#endif
