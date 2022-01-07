#include "generator.hpp"
#include <stdlib.h>
#include <math.h>
#include <omp.h>

using namespace boost::numeric::ublas;

class White : public Generator {

private:
    boost::numeric::ublas::vector<double> reference;

public:

    White(int numStates, vector<double> &ref) {
        // this->n_holes = n_holes;
        // this->n_particles = n_particles;

        this->reference = ref;
        this->numStates = numStates;

        // f1b = Format({Dense,Dense});
        // f2b = Format({Dense,Dense,Dense,Dense});
        // f3b = Format({Dense,Dense,Dense,Dense,Dense,Dense});

        // this->eta1b = vector<double>("eta1b", {numStates,numStates}, f1b);
        // this->eta2b = vector<double>("eta2b", {numStates,numStates,numStates,numStates}, f2b);
        //eta3b = vector<double>("eta2b", {numStates,numStates,numStates,numStates,numStates,numStates}, f3b);
    }

    int sgn(double x) {
        return (x > 0) ? 1 : ((x < 0) ? -1 : 0);
    }

    vector<double> compute_1b(vector<double> &f, vector<double> &Gamma, vector<double> &W) {

        // vector<double> eta1b_tensor("eta1b", {numStates,numStates}, f1b);
        // eta1b_tensor.pack();

        // double* eta1b_arr = (double*)eta1b_tensor.getStorage().getValues().getData();

        vector<double> eta1b_arr(numStates*numStates);

        double nbara, ni, nbara_sq, ni_sq;
        double denom, result;
        int a, i, idxai, idxii, idxaa, idxia, idxaiai;

        //#pragma omp parallel for private(a,i,nbara, ni, nbara_sq, ni_sq, denom, result)
        for (a = 0; a < numStates; a++) {
            for (i = 0; i <= a; i++) {
                idxai = a*numStates + i;
                idxii = i*numStates + i;
                idxaa = a*numStates + a;
                idxia = i*numStates + a;
                idxaiai = a*numStates*numStates*numStates + i*numStates*numStates + a*numStates + i;

                nbara = 1 - reference[a];
                ni = reference[i];
                nbara_sq = nbara*nbara;
                ni_sq = ni*ni;
                denom = f[idxaa]*nbara_sq - f[idxii]*ni_sq + Gamma[idxaiai]*nbara_sq*ni_sq;

                if ( abs(denom) < 1e-10 ) 
                    result = 0.25 * M_PI * sgn(f[idxai]*nbara*ni) * sgn(denom);
                else
                    result = f[idxai] * nbara*ni / denom;
            
                if (result > 100) {
                    std::cout << a << " " << i << " " << idxai << " " << result << " DENOM " << denom << " f[idxai]*nbara*ni " << f[idxai] << std::endl;
                }

                //std::cout << a << "," << i << " " << result << std::endl;
                //std::cout << result << std::endl;

                // std::cout << a << i << " RESULT " << result << std::endl;


                eta1b_arr[idxai] = result;
                eta1b_arr[idxia] = -result;

                // std::cout << a << i << " ETA 1B " << eta1b_arr[idxia] << std::endl;

            
                //std::cout << "gamma " << Gamma[idxaiai] << std::endl;

                // eta1b.insert({a,i}, result);
                // eta1b.insert({i,a}, -result);
            
            }
        }

        // for (a = 0; a < numStates; a++)
        //     for (i = 0; i < numStates; i++)
        //         std::cout << a << i << " ETA 1B " << eta1b_arr[a*numStates+i] << std::endl;
    
        return eta1b_arr;
    }

    vector<double> compute_2b(vector<double> &f, vector<double> &Gamma, vector<double> &W) {

        // vector<double> eta2b_tensor("eta2b", {numStates,numStates,numStates,numStates}, f2b);
        // eta2b_tensor.pack();

        // double* eta2b_arr = (double*)eta2b_tensor.getStorage().getValues().getData();
        
        //double* eta2b = (double*) malloc(sizeof(double)*numStates*numStates*numStates*numStates);
        //std::cout << f << "\n" << Gamma << "\n" << "\n" << reference << std::endl;

        vector<double> eta2b_arr(numStates*numStates*numStates*numStates);

        double nbara, nbarb, ni, nj, nbara_sq, nbarb_sq, ni_sq, nj_sq;
        double denom, result;
    
        int a,b,i,j, idxaa, idxbb, idxii, idxjj, idxabab, idxijij, idxaiai, idxbjbj, idxajaj, idxbibi, idxabij, idxijab;
        //#pragma omp parallel for private(a,b,i,j,nbara, nbarb, ni, nj, nbara_sq, nbarb_sq, ni_sq, nj_sq, denom, result) //shared(eta2b_arr,numStates,reference,f,Gamma) private(a,b,i,j,nbara, nbarb, ni, nj, nbara_sq, nbarb_sq, ni_sq, nj_sq, denom, result)
        for (a = 0; a < numStates; a++) {
            for (b = 0; b < numStates; b++) {
                for (i = 0; i < numStates; i++) {
                    for (j = 0; j < numStates; j++) { 

                        idxaa = a*numStates + a;
                        idxbb = b*numStates + b;
                        idxii = i*numStates + i;
                        idxjj = j*numStates + j;
                        idxabab = a*numStates*numStates*numStates + b*numStates*numStates + a*numStates + b;
                        idxijij = i*numStates*numStates*numStates + j*numStates*numStates + i*numStates + j;
                        idxaiai = a*numStates*numStates*numStates + i*numStates*numStates + a*numStates + i;
                        idxbjbj = b*numStates*numStates*numStates + j*numStates*numStates + b*numStates + j;
                        idxajaj = a*numStates*numStates*numStates + j*numStates*numStates + a*numStates + j;
                        idxbibi = b*numStates*numStates*numStates + i*numStates*numStates + b*numStates + i;
                        idxabij = a*numStates*numStates*numStates + b*numStates*numStates + i*numStates + j;
                        idxijab = i*numStates*numStates*numStates + j*numStates*numStates + a*numStates + b;
                    
                        nbara = 1 - reference[a];
                        nbarb = 1 - reference[b];                   
                        ni = reference[i];
                        nj = reference[j];
                        nbara_sq = nbara*nbara;
                        nbarb_sq = nbarb*nbarb;
                        ni_sq = ni*ni;
                        nj_sq = nj*nj;


                        denom = f[idxaa]*nbara_sq + f[idxbb]*nbarb_sq - f[idxii]*ni_sq - f[idxjj]*nj_sq \
                            + Gamma[idxabab]*nbara_sq*nbarb_sq + Gamma[idxijij]*ni_sq*nj_sq \
                            - Gamma[idxaiai]*nbara_sq*ni_sq - Gamma[idxbjbj]*nbarb_sq*nj_sq \
                            - Gamma[idxajaj]*nbara_sq*nj_sq - Gamma[idxbibi]*nbarb_sq*ni_sq;
                    
                        if ( abs(denom) < 1e-2 ) 
                            result = 0.25 * M_PI * sgn(Gamma[idxabij]*nbarb*nbara*ni*nj) * sgn(denom);
                        else
                            result = Gamma[idxabij]*nbarb*nbara*ni*nj / denom;

                        eta2b_arr[idxabij] = result;
                        //if (result != 0.0)
                        //printf("%d%d%d%d\n", a,b,i,j);

                        eta2b_arr[idxijab] = -result;

                        //printf("%d\n",a*numStates*numStates*numStates + b*numStates*numStates + i*numStates + j);
                        // eta2b_tensor.insert({a,b,i,j}, eta2b[a*numStates + b*numStates + i*numStates + j]);
                        // eta2b_tensor.insert({i,j,a,b}, eta2b[a*numStates + b*numStates + i*numStates + j]);
                    }
                }            
            }
        }
    
    
        // for(int i = 0; i < numStates*numStates*numStates*numStates; i++) {
        //     eta2b_arr[i] = eta2b[i];
        //     //printf("eta2b[%d] = %0.8f\n",i, eta2b[i]);
        // }

        //std::cout <<eta2b_arr[4*numStates*numStates*numStates + 5*numStates*numStates + 0*numStates + 1] << " " << eta2b_arr[0*numStates*numStates*numStates + 1*numStates*numStates + 4*numStates + 5] << std::endl;

        //std::cout << eta2b_tensor << std::endl;
        return eta2b_arr;

    }

    vector<double> compute_3b(vector<double> &f, vector<double> &Gamma, vector<double> &W) {
    
        return W;
    }

};
