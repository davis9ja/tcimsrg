#include "pairinghamiltonian.hpp"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

using namespace taco;
using namespace boost::numeric::ublas;
PairingHamiltonian::PairingHamiltonian(int n_holes, int n_particles, vector<double> ref, double d, double g, double pb) 
{
    this->d = d;
    this->g = g;
    this->pb = pb;

    this->n_holes = n_holes;
    this->n_particles = n_particles;

    int numStates = n_holes+n_particles;

    spBasis = (int*)malloc(numStates*sizeof(int));
    holes = (int*)malloc(n_holes*sizeof(int));
    particles = (int*)malloc(n_particles*sizeof(int));
    
    
    // Declare H1B
    // H1B = Tensor<double>("H1B", 
    //     {numStates, numStates}, 
    //     Format({Dense,Dense})
    //     );

    H1B = vector<double>(numStates*numStates);
    for (int i = 0; i < H1B.size(); i++)
        H1B[i] = 0.;

    // Delcare H2B
    // H2B = Tensor<double>("H2B", 
    //     {numStates, numStates, numStates, numStates}, 
    //     Format({Dense,Dense,Dense,Dense})
    //     );

    H2B = vector<double>(numStates*numStates*numStates*numStates);
    for (int i = 0; i < H2B.size(); i++)
        H2B[i] = 0.;

    // Fill hole and particle index containers
    double ref_i;
    for (int i = 0; i < n_holes; i++) {
        ref_i = ref[i];

        if (ref_i >= 0.5)
            holes[i] = i;
    }
        

    // Construct the vacuum Hamiltonian and fill spBasis
    double occupation = 1.0;
    int idx1b, idx2b;
    for (int p = 0; p < numStates; p++) {
        
        spBasis[p] = p;

        //H1B.insert({p,p}, this->d*(int)(p/2));
        idx1b = p*numStates + p;
        H1B[idx1b] = this->d*(int)(p/2);

        for (int q = 0; q < numStates; q++) {
            for (int r = 0; r < numStates; r++) {
                for (int s = 0; s < numStates; s++) {
                    // H2B.insert({p,q,r,s}, this->g*0.5*this->delta2B(p,q,r,s));
                    // H2B.insert({p,q,r,s}, this->pb*0.5*this->deltaPB(p,q,r,s));
                    
                    idx2b = p*numStates*numStates*numStates + q*numStates*numStates + r*numStates + s;
                    H2B[idx2b] += this->g*0.5*this->delta2B(p,q,r,s);
                    H2B[idx2b] += this->pb*0.5*this->deltaPB(p,q,r,s);

                    // if (H2B[idx2b] != 0)
                    //     printf("%d%d%d%d, %0.2f\n", p,q,r,s, H2B[idx2b]);
                }
            }
        }
    }

    //ref.pack();
    //H1B.pack();
    //H2B.pack();

    // Declare IM-SRG normal-ordered operator coefficients
    
    // One body
    // f = Tensor<double>("f", 
    //     {numStates,numStates},
    //     Format({Dense,Dense})
    //     );
    f = vector<double>(numStates*numStates);
    for(int i = 0; i < f.size(); i++)
        f[i] = 0.;

    // Two body
    Gamma = vector<double>(H2B.size());
    for (int i = 0; i < Gamma.size(); i++)
        Gamma[i] = H2B[i]; // IM-SRG(2) truncation

    // Normal order the Hamiltonian with respect to the reference
    // CANNOT use TACO to trace over indices of same tensor (why not?)

    //Tensor<double> Gamma_piqi("Gamma_piqi", {numStates, numStates}, {Dense, Dense});
    vector<double> Gamma_piqi(numStates*numStates);

    int idxGamma, idxf;
    // double Gamma_val = 0.0;
    // double ref_val = 0.0;
    for (int p = 0; p < numStates; p++) {
        for (int q = 0; q < numStates; q++) {
            double sum = 0.0;
            for (int i = 0; i < numStates; i++) {
                idxGamma = p*numStates*numStates*numStates + i*numStates*numStates + q*numStates + i;

                sum += Gamma[idxGamma]*ref[i];
            }
            //Gamma_piqi.insert({p,q}, sum);
            //Gamma_piqi[p,q] += sum;
            idxf = p*numStates+q;
            f[idxf] += H1B[idxf] + sum;
        }
    }

    //Gamma_piqi.pack();

    // Gamma_piqi.compile();
    // Gamma_piqi.assemble();
    // Gamma_piqi.compute();

    // IndexVar p,q;
    // f(p,q) = H1B(p,q) + Gamma_piqi(p,q);
    // f.evaluate();


    // Computing E zero body
    // Need to trace over 1B and 2B operators
    
    //One body trace
    double h1b_ii = 0.0;
    for (int i = 0; i < numStates; i++) {
        h1b_ii += H1B[i*numStates+i]*ref[i];     
    }

    double Gamma_ijij = 0.0;
    for (int i = 0; i < numStates; i++) {
        for (int j = 0; j < numStates; j++) {
            idxGamma = i*numStates*numStates*numStates + j*numStates*numStates + i*numStates + j;
            Gamma_ijij += Gamma[idxGamma]*ref[i]*ref[j];
        }
    }
    //std::cout << "Gamma_ijij, " << Gamma_ijij << std::endl;
    
    E = h1b_ii + 0.5*Gamma_ijij;//double E_val = h1b_ii + 0.5*Gamma_ijij;

    // E = Tensor<double>(E_val);
    // E.pack();

    // Tensor<double> E(E_val);
    // E.pack();

    // for (int p = 0; p < numStates; p++) {
    //     for (int q = 0; q < numStates; q++) {
    //         f->insert({p,q}, *H1B(p,q) + Gamma_piqi);
    //     }
    // }
}

int PairingHamiltonian::delta2B(int p, int q, int r, int s) {

    int pp = floor(p/2.);
    int qp = floor(q/2.);
    int rp = floor(r/2.);
    int sp = floor(s/2.);

    int ps = (p%2 == 0) ? 1 : -1;   
    int qs = (q%2 == 0) ? 1 : -1;
    int rs = (r%2 == 0) ? 1 : -1;
    int ss = (s%2 == 0) ? 1 : -1;

    if ( pp != qp || rp != sp ) 
        return 0;
    if ( ps == qs || rs == ss ) 
        return 0;
    if ( ps == rs && qs == ss ) 
        return -1;
    if ( ps == ss && qs == rs ) 
        return 1;

    return 0;

}

int PairingHamiltonian::deltaPB(int p, int q, int r, int s) {

    int pp = floor(p/2.);
    int qp = floor(q/2.);
    int rp = floor(r/2.);
    int sp = floor(s/2.);

    int ps = (p%2 == 0) ? 1 : -1;   
    int qs = (q%2 == 0) ? 1 : -1;
    int rs = (r%2 == 0) ? 1 : -1;
    int ss = (s%2 == 0) ? 1 : -1;

    if ( (pp != qp && rp == sp) || (pp == qp && rp != sp) ) {
        if ( ps == qs || rs == ss )
            return 0;
        if ( ps == rs && qs == ss )
            return -1;
        if ( ps == ss && qs == rs )
            return 1;
    }

    return 0;

}


// string PairingHamiltonian::get_1B_str() {
//     int numStates = n_holes + n_particles;
//     for (int p = 0; p < numStates; p++) {
        
//     }
// }

PairingHamiltonian::~PairingHamiltonian() {
    //free(ref);
    free(holes);
    free(particles);
    free(spBasis);
    //free(H1B);
    //free(H2B);
}
