#include "generator.hpp"
#include "taco.h"
#include "imsrg_utils.hpp"
#include "occupation_factors.hpp"

#include <stdlib.h>
#include <math.h>
#include <omp.h>

using namespace boost::numeric::ublas;
using namespace taco;

class Brillouin : public Generator {

private:
    vector<double> reference;
    vector<double> rho1b, rho2b, rho3b;

    Tensor<double> 
        occA_a, occA_b, 
        occB_a, occB_b,
        occC_a, occC_b, occC_c,
        occC2_a, occC2_b, occC2_c,
        occD_a, occD_b, occD_c, occD_d,
        occD2_a, occD2_b, occD2_c, occD2_d;


public:
    Brillouin() {}

    Brillouin(int numStates, vector<double> &ref, vector<double> &rho1b, vector<double> &rho2b, vector<double> &rho3b, OccupationFactors *occFact) {
        this->numStates = numStates;
        this->reference = reference;

        this->rho1b = rho1b;
        this->rho2b = rho2b;
        this->rho3b = rho3b;

        std::string factor_path = occFact->getPath();
        occFact->readOccTensors(factor_path,
                                occA_a, occA_b,
                                occB_a, occB_b,
                                occC_a, occC_b, occC_c,
                                occC2_a, occC2_b, occC2_c,
                                occD_a, occD_b, occD_c, occD_d,
                                occD2_a, occD2_b, occD2_c, occD2_d);


    }

    vector<double> compute_1b(vector<double> &f, vector<double> &Gamma, vector<double> &W) {

        std::vector<int> shape_r2 = {numStates,numStates};
        std::vector<int> shape_r4 = {numStates,numStates,numStates,numStates};
        Format format_r2({Dense,Dense});
        Format format_r4({Dense,Dense,Dense,Dense});


        vector<double> eta1b_arr(numStates*numStates);
        for (int i = 0; i  < eta1b_arr.size(); i++)
            eta1b_arr[i] = 0.0;


        Tensor<double> term1_1b(shape_r2, format_r2), sum2_1b(shape_r2, format_r2);

        Tensor<double> 
            f_t("f_t", shape_r2, format_r2),
            Gamma_t("Gamma_t", shape_r4, format_r4),
            eta1b_t("eta1b_t", shape_r2, format_r2),
            rho2b_t("rho2b_t", shape_r4, format_r4);
        

        f_t.pack();
        Gamma_t.pack();
        eta1b_t.pack();
        rho2b_t.pack();

        vector2tensor(f, f_t);
        vector2tensor(Gamma, Gamma_t);
        vector2tensor(rho2b, rho2b_t);
        //vector2tensor(eta1b, eta1b_t);


        IndexVar i,j,a,b,c,p;    

        term1_1b(i,j) = occA_a(j,p)*occA_b(p,i)*f_t(j,i);
        sum2_1b(i,j) = Gamma_t(j,a,b,c)*rho2b_t(i,a,b,c) - Gamma_t(a,b,i,c)*rho2b_t(a,b,j,c);

        
        eta1b_t(i,j) = term1_1b(i,j); //- 0.5*sum2_1b(i,j);
        
        eta1b_t.evaluate();
        
        tensor2vector(eta1b_t, eta1b_arr);

        return eta1b_arr;


    }


    vector<double> compute_2b(vector<double> &f, vector<double> &Gamma, vector<double> &W) {

        std::vector<int> shape_r2 = {numStates,numStates};
        std::vector<int> shape_r4 = {numStates,numStates,numStates,numStates};
        std::vector<int> shape_r6 = {numStates,numStates,numStates,numStates,numStates,numStates};
        
        Format format_r2({Dense,Dense});
        Format format_r4({Dense,Dense,Dense,Dense});
        Format format_r6({Dense,Dense,Dense,Dense,Dense,Dense});

        vector<double> eta2b_arr(numStates*numStates*numStates*numStates);
        for (int i = 0; i  < eta2b_arr.size(); i++)
            eta2b_arr[i] = 0.0;


        Tensor<double> 
            term1_2b(shape_r4, format_r4), 
            sum2_2b(shape_r4, format_r4),
            sum3_2b(shape_r4, format_r4),
            sum4_2b(shape_r4, format_r4),
            sum5_2b(shape_r4, format_r4),
            LambdaGamma(shape_r4, format_r4),
            GammaLambda(shape_r4, format_r4);
            

        Tensor<double> 
            f_t("f_t", shape_r2, format_r2),
            Gamma_t("Gamma_t", shape_r4, format_r4),
            eta2b_t("eta2b_t", shape_r4, format_r4),
            rho2b_t("rho2b_t", shape_r4, format_r4),
            rho3b_t("rho3b_t", shape_r6, format_r6);        

        f_t.pack();
        Gamma_t.pack();
        eta2b_t.pack();
        rho2b_t.pack();
        rho3b_t.pack();

        vector2tensor(f, f_t);
        vector2tensor(Gamma, Gamma_t);
        vector2tensor(rho2b, rho2b_t);
        vector2tensor(rho3b, rho3b_t);

        vector2tensor(eta2b, eta2b_t);


        IndexVar i,j,k,l,a,b,c,p,q,r,s;    
        
        LambdaGamma(i,j,k,l) = rho2b_t(i,j,a,b)*Gamma_t(a,b,k,l);
        GammaLambda(i,j,k,l) = Gamma_t(i,j,a,b)*rho2b_t(a,b,k,l);

        term1_2b(i,j,k,l) = Gamma_t(k,l,i,j)*occD2_a(i,p)*occD2_b(p,j,q)*occD2_c(q,k,r)*occD2_d(r,l);
        sum2_2b(i,j,k,l) = f_t(a,i)*rho2b_t(a,j,k,l) - f_t(a,j)*rho2b_t(a,i,k,l) - f_t(k,a)*rho2b_t(i,j,a,l) + f_t(l,a)*rho2b_t(i,j,a,k);
        sum3_2b(i,j,k,l) = LambdaGamma(k,l,i,j)*occB_a(i,p)*occB_b(p,j) - GammaLambda(k,l,i,j)*occB_a(k,q)*occB_b(q,l);
        sum4_2b(i,j,k,l) = \
            occA_a(j,p)*occA_b(p,k)*Gamma_t(a,k,c,j)*rho2b_t(a,i,c,l) - \
            occA_a(i,q)*occA_b(q,k)*Gamma_t(a,k,c,i)*rho2b_t(a,j,c,l) - \
            occA_a(j,r)*occA_b(r,l)*Gamma_t(a,l,c,j)*rho2b_t(a,i,c,k) + \
            occA_a(i,s)*occA_b(s,l)*Gamma_t(a,l,c,i)*rho2b_t(a,j,c,k);
        sum5_2b(i,j,k,l) = \
            Gamma_t(k,a,b,c)*rho3b_t(a,i,j,b,c,l) - Gamma_t(l,a,b,c)*rho3b_t(a,i,j,b,c,k) - \
            Gamma_t(a,b,i,c)*rho3b_t(a,b,j,c,k,l) + Gamma_t(a,b,j,c)*rho3b_t(a,b,i,c,k,l);

        eta2b_t(i,j,k,l) = sum2_2b(i,j,k,l); + 0.5*sum3_2b(i,j,k,l) + sum4_2b(i,j,k,l) + 0.5*sum5_2b(i,j,k,l);
            // +                                       \
              //+ ; //+ ; //+ ;

        eta2b_t.evaluate();

        tensor2vector(eta2b_t, eta2b_arr);

        return eta2b_arr;


    }

    vector<double> compute_3b(vector<double> &f, vector<double> &Gamma, vector<double> &W) {
    
        return W;
    }



};
