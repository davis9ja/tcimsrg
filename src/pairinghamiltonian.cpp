#include "pairinghamiltonian.hpp"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

using namespace taco;
using namespace boost::numeric::ublas;

PairingHamiltonian::PairingHamiltonian(){}

PairingHamiltonian::PairingHamiltonian(int numStates, vector<double> rho1b, vector<double> rho2b, double d, double g, double pb) 
{
    this->d = d;
    this->g = g;
    this->pb = pb;
    this->numStates = numStates;

    // this->n_holes = n_holes;
    // this->n_particles = n_particles;

    //int numStates = n_holes+n_particles;

    // spBasis = (int*)malloc(numStates*sizeof(int));
    // holes = (int*)malloc(n_holes*sizeof(int));
    // particles = (int*)malloc(n_particles*sizeof(int));

    // Allocate Gamma and f
    this->Gamma = vector<double>(numStates*numStates*numStates*numStates);
    this->f = vector<double>(numStates*numStates);
    this-> E = 0.0;
    
    // Initialize H1B
    H1B = vector<double>(numStates*numStates);
    for (int i = 0; i < H1B.size(); i++)
        H1B[i] = 0.;

    // Initialize H2B
    H2B = vector<double>(numStates*numStates*numStates*numStates);
    for (int i = 0; i < H2B.size(); i++)
        H2B[i] = 0.;

    // Fill hole and particle index containers
    // double ref_i;
    // for (int i = 0; i < n_holes; i++) {
    //     ref_i = ref[i];

    //     if (ref_i >= 0.5)
    //         holes[i] = i;
    // }
        

    // Construct the vacuum Hamiltonian and fill spBasis
    double occupation = 1.0;
    int idx1b, idx2b;
    for (int p = 0; p < numStates; p++) {
        
        //spBasis[p] = p;

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
    
    normal_order_ensRef(rho1b, rho2b);

}

PairingHamiltonian::PairingHamiltonian(int nHoles, int nParticles, vector<double> singleRef, double d, double g, double pb) 
{
    this->d = d;
    this->g = g;
    this->pb = pb;
    this->numStates = nHoles+nParticles;

    // this->n_holes = n_holes;
    // this->n_particles = n_particles;

    //int numStates = n_holes+n_particles;

    // spBasis = (int*)malloc(numStates*sizeof(int));
    // holes = (int*)malloc(n_holes*sizeof(int));
    // particles = (int*)malloc(n_particles*sizeof(int));
    
    vector<int> holes(nHoles), particles(nParticles);

    // Allocate Gamma and f
    Gamma = vector<double>(numStates*numStates*numStates*numStates);
    f = vector<double>(numStates*numStates);
    
    
    // Initialize H1B
    H1B = vector<double>(numStates*numStates);
    for (int i = 0; i < H1B.size(); i++)
        H1B[i] = 0.;

    // Initialize H2B
    H2B = vector<double>(numStates*numStates*numStates*numStates);
    for (int i = 0; i < H2B.size(); i++)
        H2B[i] = 0.;

    // Fill hole and particle index containers
    int ref_i;
    for (int i = 0; i < numStates; i++) {
        ref_i = (int)singleRef[i];

        if (ref_i == 1)
            holes[i] = i;
        else 
            particles[i-nHoles] = i;
                    
    }

    // Construct the vacuum Hamiltonian and fill spBasis
    double occupation = 1.0;
    int idx1b, idx2b;
    for (int p = 0; p < numStates; p++) {
        
        //spBasis[p] = p;

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
    
    normal_order_singleRef(holes, particles);

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

void PairingHamiltonian::normal_order_ensRef(vector<double> rho1b, vector<double> rho2b) {

    int idx1b, idxGamma, idxf, idxRho1b;;

    // Declare IM-SRG normal-ordered operator coefficients

    // One body
    for(int i = 0; i < f.size(); i++)
        f[i] = 0.;

    // Two body
    for (int i = 0; i < Gamma.size(); i++)
        Gamma[i] = H2B[i]; // IM-SRG(2) truncation

    // Normal order the Hamiltonian with respect to the reference
    // CANNOT use TACO to trace over indices of same tensor (why not?)

    //Tensor<double> Gamma_piqi("Gamma_piqi", {numStates, numStates}, {Dense, Dense});
    //vector<double> Gamma_piqi(numStates*numStates);

    // double Gamma_val = 0.0;
    // double ref_val = 0.0;
    for (int p = 0; p < numStates; p++) {
        for (int q = 0; q < numStates; q++) {
            double sum = 0.0;
            for (int i = 0; i < numStates; i++) {
                for (int j = 0; j < numStates; j++) {
                    idxGamma = p*numStates*numStates*numStates + i*numStates*numStates + q*numStates + j;
                    idxRho1b = i*numStates + j;
                    sum += Gamma[idxGamma]*rho1b[idxRho1b];
                }
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
    double h1b_ij = 0.0;
    for (int i = 0; i < numStates; i++) {
        for (int j = 0; j < numStates; j++) {
            idx1b = i*numStates+j;
            h1b_ij += H1B[idx1b]*rho1b[idx1b];     
        }
    }

    double Gamma_ijkl = 0.0;
    for (int i = 0; i < numStates; i++) {
        for (int j = 0; j < numStates; j++) {
            for (int k = 0; k < numStates; k++) {
                for (int l = 0; l < numStates; l++) {
                    idxGamma = i*numStates*numStates*numStates + j*numStates*numStates + k*numStates + l;
                    Gamma_ijkl += Gamma[idxGamma]*rho2b[idxGamma];
                }
            }
        }
    }
    //std::cout << "Gamma_ijij, " << Gamma_ijij << std::endl;
    
    E = h1b_ij + 0.25*Gamma_ijkl;//double E_val = h1b_ii + 0.5*Gamma_ijij;
    
}

void PairingHamiltonian::normal_order_singleRef(vector<int> holes, vector<int> particles) {
    for(int i = 0; i < f.size(); i++)
        f[i] = H1B[i];
    for (int i = 0; i < H2B.size(); i++)
        Gamma[i] = H2B[i];

    E = 0.0;
    for (int i : holes)
        E += H1B[i*numStates + i];

    for (int i : holes)
        for (int j : holes)
            E += 0.5*H2B[i*pow(numStates,3) + j*pow(numStates,2) + i*numStates + j];

    for (int i = 0; i < numStates; i++)
        for (int j = 0; j < numStates; j++) 
            for (int h : holes)  
                f[i*numStates+j] += H2B[i*pow(numStates,3) + h*pow(numStates,2) + j*numStates + h];
            

}
