#include "BACKEND.hpp"
#include <stdio.h>
//#include <math.h>
#include <boost/numeric/ublas/io.hpp>

using namespace boost::numeric::ublas;

class Backend_UBLAS : public Backend {

private:
    int numStates;
    
    int idxr2(int i, int j) {
        return i*numStates + j;
    }
    
    int idxr3(int i, int j, int k) {
        return i*numStates*numStates + j*numStates + k;
    }

    int idxr4(int i, int j, int k, int l) {
        return i*numStates*numStates*numStates + j*numStates*numStates + k*numStates + l;
    }

    vector<double> occA, occB, occC, occD;

public:
    Backend_UBLAS(int numStates, OccupationFactors *occFact) {        

        std::string path = (*occFact).getPath();

        this->numStates = numStates;
        this->occFact = occFact;
        this->occA = vector<double>(numStates*numStates);
        this->occB = vector<double>(numStates*numStates);
        this->occC = vector<double>(numStates*numStates*numStates);        
        this->occD = vector<double>(numStates*numStates*numStates*numStates);        
        
        (*occFact).contractOccTensors(path, this->occA, this->occB, this->occC, this->occD);
        
    }

    double flow_0b(vector<double> &occA, vector<double> &occD, 
                   vector<double> &f, vector<double> &Gamma, 
                   vector<double> &eta1b, vector<double> &eta2b) {
             
        // Return
        double dE = 0;

        // TERM 1 VARS
        //double sum1_0b = 0;
        int ab, ba;

        // TERM 2 VARS
        //double sum2_0b = 0;
        int abcd, cdab;


        // TERM 1 ZERO BODY
        for(int a = 0; a < numStates; a++) {
            for (int b = 0; b < numStates; b++) {
                // ab = a*numStates+b;
                // ba = b*numStates+a;
                ab = idxr2(a,b);
                ba = idxr2(b,a);

                dE += this->occA[ab]*eta1b[ab]*f[ba];                
            }
        }

        // TERM 2 ZERO BODY
        for(int a = 0; a < numStates; a++) {
            for (int b = 0; b < numStates; b++) {
                for (int c = 0; c < numStates; c++) {
                    for (int d = 0; d < numStates; d++) {
                        // abcd = a*numStates*numStates*numStates + b*numStates*numStates + c*numStates + d;
                        // cdab = c*numStates*numStates*numStates + d*numStates*numStates + a*numStates + b;
                        abcd = idxr4(a,b,c,d);
                        cdab = idxr4(c,d,a,b);

                        dE += 0.5*eta2b[abcd]*Gamma[cdab]*this->occD[abcd];
                    }
                }
            }
        }

        //dE = sum1_0b + 0.5*sum2_0b;
        
        return dE;
    }

    vector<double> flow_1b(vector<double> &occA, vector<double> &occC, 
                                          vector<double> &f, vector<double> &Gamma, 
                                          vector<double> &eta1b, vector<double> &eta2b) {

        // Return 
        vector<double> df(numStates*numStates);
        for(int i = 0; i < df.size(); i++)
            df[i] = 0.;

        int ij;


        // TERM 1 VARS
        int ia, aj, ja, ai;

        // TERM 2 VARS
        int ab, biaj;

        // TERM 3 VARS
        int abc, ciab, abcj, cjab, abci;
        
        for (int i = 0; i < numStates; i++) {
            for (int j = 0; j < numStates; j++) {
                ij = i*numStates + j;

                // TERM 1 ONE BODY
                for (int a = 0; a < numStates; a++) {
                    // ia = i*numStates + a;
                    // aj = a*numStates + j;
                    // ja = j*numStates + a;
                    // ai = a*numStates + i;

                    ia = idxr2(i,a);
                    aj = idxr2(a,j);
                    ja = idxr2(j,a);
                    ai = idxr2(a,i);                    

                    df[ij] += eta1b[ia]*f[aj] + eta1b[ja]*f[ai];
                    
                }

                // TERM 2 ONE BODY
                for (int a = 0; a < numStates; a++) {
                    for (int b = 0; b < numStates; b++) {
                        // ab = a*numStates + b;
                        // biaj = b*numStates*numStates*numStates + i*numStates*numStates + a*numStates + j;
                        
                        ab = idxr2(a,b);
                        biaj = idxr4(b,i,a,j);

                        df[ij] += this->occA[ab]*(eta1b[ab]*Gamma[biaj] - f[ab]*eta2b[biaj]);
                    }
                }

                // TERM 3 ONE BODY
                for (int a = 0; a < numStates; a++) {
                    for (int b = 0; b < numStates; b++) {
                        for (int c = 0; c < numStates; c++) {
                            // abc = a*numStates*numStates + b*numStates + c;
                            // ciab = c*numStates*numStates*numStates + i*numStates*numStates + a*numStates + b;
                            // abcj = a*numStates*numStates*numStates + b*numStates*numStates + c*numStates + j;
                            // cjab = c*numStates*numStates*numStates + j*numStates*numStates + a*numStates + b;
                            // abci = a*numStates*numStates*numStates + b*numStates*numStates + c*numStates + i;

                            abc = idxr3(a,b,c);
                            ciab = idxr4(c,i,a,b);
                            abcj = idxr4(a,b,c,j);
                            cjab = idxr4(c,j,a,b);
                            abci = idxr4(a,b,c,i);

                            df[ij] += 0.5*this->occC[abc]*(eta2b[ciab]*Gamma[abcj] + eta2b[cjab]*Gamma[abci]);
                        }
                    }
                }
                // if ( abs(df[ij]) > 100 ) {
                //     std::cout << "FAILURE ON DF "<<i<<","<<j<<" = "<<df[ij] << std::endl;
                //     exit(1);
                // }

            }
        }


        return df;
        
    }

    vector<double> flow_2b(vector<double> &occA, vector<double> &occB, 
                                          vector<double> &f, vector<double> &Gamma, 
                                          vector<double> &eta1b, vector<double> &eta2b) {
        // Return
        vector<double> dGamma(numStates*numStates*numStates*numStates);
        for (int i =0; i < dGamma.size(); i++)
            dGamma[i] = 0.;

        int ijkl;

        // TERM 1 VARS
        int ia, ajkl, ja, aikl, ak, ijal, al, ijak;

        // TERM 2 VARS
        int ab, ijab, abkl;

        // TERM 3 VARS
        int bjal, aibk, bial, ajbk, bjak, aibl, biak, ajbl;

        for (int i = 0; i < numStates; i++) {
            for (int j = 0; j < numStates; j++) {
                for (int k = 0; k < numStates; k++) {
                    for (int l = 0; l < numStates; l++) {
                        ijkl = idxr4(i,j,k,l);
                        
                        // TERM 1 TWO BODY
                        for (int a = 0; a < numStates; a++) {
                            ia = idxr2(i,a);
                            ajkl = idxr4(a,j,k,l);
                            ja = idxr2(j,a);
                            aikl = idxr4(a,i,k,l);
                            ak = idxr2(a,k);
                            ijal = idxr4(i,j,a,l);
                            al = idxr2(a,l);
                            ijak = idxr4(i,j,a,k);

                            dGamma[ijkl] += 
                                (eta1b[ia]*Gamma[ajkl] - f[ia]*eta2b[ajkl]) - 
                                (eta1b[ja]*Gamma[aikl] - f[ja]*eta2b[aikl]) - 
                                (eta1b[ak]*Gamma[ijal] - f[ak]*eta2b[ijal]) + 
                                (eta1b[al]*Gamma[ijak] - f[al]*eta2b[ijak]);
                        }

                        // TERM 2/3 TWO BODY
                        for (int a = 0; a < numStates; a++) {
                            for (int b = 0; b < numStates; b++) {
                                ab = idxr2(a,b);
                                ijab = idxr4(i,j,a,b);
                                abkl = idxr4(a,b,k,l);

                                bjal = idxr4(b,j,a,l);
                                aibk = idxr4(a,i,b,k);
                                bial = idxr4(b,i,a,l);
                                ajbk = idxr4(a,j,b,k);
                                bjak = idxr4(b,j,a,k);
                                aibl = idxr4(a,i,b,l);
                                biak = idxr4(b,i,a,k);
                                ajbl = idxr4(a,j,b,l);

                                //TERM 2
                                dGamma[ijkl] += 0.5*
                                    this->occB[ab]*(eta2b[ijab]*Gamma[abkl] - Gamma[ijab]*eta2b[abkl]);

                                //TERM 3
                                dGamma[ijkl] -=
                                    this->occA[ab]*(eta2b[bjal]*Gamma[aibk] -
                                                    eta2b[bial]*Gamma[ajbk] -
                                                    eta2b[bjak]*Gamma[aibl] +
                                                    eta2b[biak]*Gamma[ajbl]);
                                
                            }
                        }
                    }
                }
            }
        }

        return dGamma;
    }
    
};
