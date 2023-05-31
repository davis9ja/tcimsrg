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
        occD_a, occD_b, occD_c, occD_d;

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
                                occD_a, occD_b, occD_c, occD_d);


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

        
        eta1b_t(i,j) = term1_1b(i,j) - 0.5*sum2_1b(i,j);

        
        tensor2vector(eta1b_t, eta1b);

        return eta1b;


    }



}
