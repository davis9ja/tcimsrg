#include "BACKEND.hpp"
#include "taco.h"
#include "imsrg_utils.hpp"
#include <omp.h>

using namespace boost::numeric::ublas;
using namespace taco;

class Backend_TACO : public Backend {

private:

    int numStates;

    Tensor<double> 
        occA_a, occA_b, 
        occB_a, occB_b,
        occC_a, occC_b, occC_c,
        occC2_a, occC2_b, occC2_c,
        occD_a, occD_b, occD_c, occD_d,
        occD2_a, occD2_b, occD2_c, occD2_d;

    
    // void vector2tensor(vector<double> &input, Tensor<double> &output) {
    //     double* output_arr = (double*)output.getStorage().getValues().getData();
    //     for (int i = 0; i < input.size(); i++)
    //         output_arr[i] = input[i];
    // }

    // void tensor2vector(Tensor<double> &input, vector<double> &output) {
    //     double* input_arr = (double*)input.getStorage().getValues().getData();
    //     for (int i = 0; i < output.size(); i++)
    //         output[i] = input_arr[i];
    // }

public:

    using Backend::flow_0b;
    using Backend::flow_1b;
    using Backend::flow_2b;

    Backend_TACO(){}
    Backend_TACO(int numStates, OccupationFactors *occFact) {

        std::string factor_path = occFact->getPath();
        this->numStates = numStates;

        occFact->readOccTensors(factor_path,
                                occA_a, occA_b,
                                occB_a, occB_b,
                                occC_a, occC_b, occC_c,
                                occC2_a, occC2_b, occC2_c,
                                occD_a, occD_b, occD_c, occD_d,
                                occD2_a, occD2_b, occD2_c, occD2_d);
    }
    
    //*********** SINGLE REFERENCE FLOW EQUATIONS

    double flow_0b(vector<double> &f, vector<double> &Gamma, 
                   vector<double> &eta1b, vector<double> &eta2b) {
        int numStates = (int)sqrt(f.size());
        std::vector<int> shape_r2 = {numStates,numStates};
        std::vector<int> shape_r3 = {numStates,numStates,numStates};
        std::vector<int> shape_r4 = {numStates,numStates,numStates,numStates};
        Format format_r2({Dense,Dense});
        Format format_r3({Dense,Dense,Dense});
        Format format_r4({Dense,Dense,Dense,Dense});

        double dE;
        Tensor<double> sum1_0b, sum2_0b, dE_t;

        Tensor<double> 
            f_t("f_t", shape_r2, format_r2),
            Gamma_t("Gamma_t", shape_r4, format_r4),
            eta1b_t("eta1b_t", shape_r2, format_r2),
            //eta1b_t_OT("eta1b_t_OT_0b", shape_r2, format_r2),
            eta2b_t("eta2b_t", shape_r4, format_r4);
            //eta2b_t_OT("eta2b_t_OT_0b", shape_r4, format_r4);

        f_t.pack();
        Gamma_t.pack();
        eta1b_t.pack();
        eta2b_t.pack();
        // eta1b_t_OT.pack();
        // eta2b_t_OT.pack();

        vector2tensor(f, f_t);
        vector2tensor(Gamma, Gamma_t);
        vector2tensor(eta1b, eta1b_t);
        vector2tensor(eta2b, eta2b_t);

        IndexVar a,b,c,d,p,q;    
    
        sum1_0b() = occA_a(a,p)*occA_b(p,b)*eta1b_t(a,b)*f_t(a,b);
        sum2_0b() = occD_a(a,p)*occD_b(p,b)*occD_c(c,q)*occD_d(q,d)*eta2b_t(a,b,c,d)*Gamma_t(c,d,a,b);

        // eta1b_t_OT(a,b) = occA_a(a,p)*occA_b(p,b)*eta1b_t(a,b);
        // eta2b_t_OT(a,b,c,d) = occD_a(a,p)*occD_b(p,b)*occD_c(c,q)*occD_d(q,d)*eta2b_t(a,b,c,d);

        // sum1_0b() = eta1b_t_OT(a,b)*f_t(a,b);
        // sum2_0b() = eta2b_t_OT(a,b,c,d)*Gamma_t(c,d,a,b);

        dE_t() = sum1_0b() + 0.5*sum2_0b();

        // sum1_0b.evaluate();
        // sum2_0b.evaluate();
        dE_t.evaluate();
    
        dE = static_cast<const double*>(dE_t.getStorage().getValues().getData())[0];

        return dE;

    }

    vector<double> flow_1b(vector<double> &f, vector<double> &Gamma, 
                           vector<double> &eta1b, vector<double> &eta2b) {

        int numStates = (int)sqrt(f.size());
        std::vector<int> shape_r2 = {numStates,numStates};
        std::vector<int> shape_r4 = {numStates,numStates,numStates,numStates};
        Format format_r2({Dense,Dense});
        Format format_r4({Dense,Dense,Dense,Dense});

        vector<double> df(f.size());

        Tensor<double> 
            sum1_1b(shape_r2, format_r2), 
            sum2_1b(shape_r2, format_r2), 
            sum3_1b(shape_r2, format_r2), 
            df_t(shape_r2, format_r2);

        Tensor<double> 
            f_t("f_t", shape_r2, format_r2),
            Gamma_t("Gamma_t", shape_r4, format_r4),
            eta1b_t("eta1b_t", shape_r2, format_r2),
            eta2b_t("eta2b_t", shape_r4, format_r4);
            // eta1b_t_OT("eta1b_t_OT_1b", shape_r2, format_r2),
            // f_t_OT("f_t_OT_1b", shape_r2, format_r2),
            // eta2b_t_OT("eta2b_t_OT_1b", shape_r4, format_r4);
    
        f_t.pack();
        Gamma_t.pack();
        eta1b_t.pack();
        eta2b_t.pack();
        // eta1b_t_OT.pack();
        // f_t_OT.pack();
        // eta2b_t_OT.pack();

        vector2tensor(f, f_t);
        vector2tensor(Gamma, Gamma_t);
        vector2tensor(eta1b, eta1b_t);
        vector2tensor(eta2b, eta2b_t);
        
        IndexVar i,j,k,l,a,b,c,p,q;
        // eta1b_t_OT(i,j) = occA_a(i,p)*occA_b(p,j)*eta1b_t(i,j);
        // f_t_OT(i,j) =  occA_a(i,p)*occA_b(p,j)*f_t(i,j);
        // eta2b_t_OT(i,j,k,l) = occC_a(i,p)*occC_b(p,j,q)*occC_c(q,k)*eta2b_t(k,l,i,j);

        // sum1_1b(i,j) = eta1b_t(i,a)*f_t(a,j) + eta1b_t(j,a)*f_t(a,i);
        // sum2_1b(i,j) = eta1b_t_OT(a,b)*Gamma_t(b,i,a,j) - f_t_OT(a,b)*eta2b_t(b,i,a,j);
        // sum3_1b(i,j) = eta2b_t_OT(c,i,a,b)*Gamma_t(a,b,c,j) + eta2b_t_OT(c,j,a,b)*Gamma_t(a,b,c,i);

        sum1_1b(i,j) = eta1b_t(i,a)*f_t(a,j) + eta1b_t(j,a)*f_t(a,i);
        sum2_1b(i,j) = occA_a(a,p)*occA_b(p,b)*eta1b_t(a,b)*Gamma_t(b,i,a,j) - occA_a(a,p)*occA_b(p,b)*f_t(a,b)*eta2b_t(b,i,a,j);
        sum3_1b(i,j) = occC_a(a,p)*occC_b(p,b,q)*occC_c(q,c)*eta2b_t(c,i,a,b)*Gamma_t(a,b,c,j) + occC_a(a,p)*occC_b(p,b,q)*occC_c(q,c)*eta2b_t(c,j,a,b)*Gamma_t(a,b,c,i);
    
        df_t(i,j) = sum1_1b(i,j) + sum2_1b(i,j) + 0.5*sum3_1b(i,j);
    
        df_t.evaluate();

        // for (int a = 0; a <numStates; a++) {
        //     for (int i = 0; i < numStates; i++) {
        //         if (df_t(a,i) > 1e10)
        //             std::cout << " BIG F VALUE " << df_t(a,i) << std::endl;
        //     }
        // }

        tensor2vector(df_t, df);

        return df;

    }

    vector<double> flow_2b(vector<double> &f, vector<double> &Gamma, 
                           vector<double> &eta1b, vector<double> &eta2b) {
        int numStates = (int)sqrt(f.size());
        std::vector<int> shape_r2 = {numStates,numStates};
        std::vector<int> shape_r4 = {numStates,numStates,numStates,numStates};
        Format format_r2({Dense,Dense});
        Format format_r4({Dense,Dense,Dense,Dense});

        vector<double> dGamma(Gamma.size());

        Tensor<double> 
            sum1_2b(shape_r4, format_r4), 
            sum2_2b(shape_r4, format_r4), 
            sum3_2b(shape_r4, format_r4), 
            dGamma_t(shape_r4, format_r4);

        Tensor<double> 
            f_t("f_t", shape_r2, format_r2),
            Gamma_t("Gamma_t", shape_r4, format_r4),
            eta1b_t("eta1b_t", shape_r2, format_r2),
            eta2b_t("eta2b_t", shape_r4, format_r4);
            // eta2b_t_OT("eta2b_t_OT_2b", shape_r4, format_r4),
            // Gamma_t_OT("Gamma_t_OT_2b", shape_r4, format_r4);

        f_t.pack();
        Gamma_t.pack();
        eta1b_t.pack();
        eta2b_t.pack();

        vector2tensor(f, f_t);
        vector2tensor(Gamma, Gamma_t);
        vector2tensor(eta1b, eta1b_t);
        vector2tensor(eta2b, eta2b_t);

        IndexVar i,j,k,l,a,b,p;
        // eta2b_t_OT(i,j,a,b) = occB_a(a,p)*occB_b(p,b)*eta2b_t(i,j,a,b);
        // Gamma_t_OT(i,j,a,b) = occB_a(a,p)*occB_b(p,b)*Gamma_t(i,j,a,b);

        sum1_2b(i,j,k,l) = eta1b_t(i,a)*Gamma_t(a,j,k,l) - f_t(i,a)*eta2b_t(a,j,k,l) - eta1b_t(j,a)*Gamma_t(a,i,k,l) + f_t(j,a)*eta2b_t(a,i,k,l) \
            - eta1b_t(a,k)*Gamma_t(i,j,a,l) + f_t(a,k)*eta2b_t(i,j,a,l) + eta1b_t(a,l)*Gamma_t(i,j,a,k) - f_t(a,l)*eta2b_t(i,j,a,k);

        sum2_2b(i,j,k,l) = occB_a(a,p)*occB_b(p,b)*eta2b_t(i,j,a,b)*Gamma_t(a,b,k,l) - occB_a(a,p)*occB_b(p,b)*Gamma_t(i,j,a,b)*eta2b_t(a,b,k,l);

        sum3_2b(i,j,k,l) = occA_a(a,p)*occA_b(p,b)*eta2b_t(b,j,a,l)*Gamma_t(a,i,b,k) - occA_a(a,p)*occA_b(p,b)*eta2b_t(b,i,a,l)*Gamma_t(a,j,b,k) \
            - occA_a(a,p)*occA_b(p,b)*eta2b_t(b,j,a,k)*Gamma_t(a,i,b,l) + occA_a(a,p)*occA_b(p,b)*eta2b_t(b,i,a,k)*Gamma_t(a,j,b,l);

        dGamma_t(i,j,k,l) = sum1_2b(i,j,k,l) + 0.5*sum2_2b(i,j,k,l) - sum3_2b(i,j,k,l);
    
        dGamma_t.evaluate();

        tensor2vector(dGamma_t, dGamma);

        return dGamma;

    }

    //******* MULTI-REFERENCE FLOW EQUATIONS

    //** up to 2B density
    double flow_0b(vector<double> &f, vector<double> &Gamma, 
                   vector<double> &eta1b, vector<double> &eta2b,
                   vector<double> &rho1b, vector<double> &rho2b, vector<double> &dGamma) {

        int numStates = (int)sqrt(f.size());
        std::vector<int> shape_r2 = {numStates,numStates};
        std::vector<int> shape_r3 = {numStates,numStates,numStates};
        std::vector<int> shape_r4 = {numStates,numStates,numStates,numStates};
        Format format_r2({Dense,Dense});
        Format format_r3({Dense,Dense,Dense});
        Format format_r4({Dense,Dense,Dense,Dense});

        double dE;
        Tensor<double> sum1_0b("sum1_0b"), sum2_0b1("sum2_0b1"), sum2_0b2("sum2_0b2"), sum3_0b("sum3_0b"), dE_t("dE");

        Tensor<double> 
            f_t("f_t", shape_r2, format_r2),
            Gamma_t("Gamma_t", shape_r4, format_r4),
            dGamma_t("dGamma_t", shape_r4, format_r4),
            eta1b_t("eta1b_t", shape_r2, format_r2),
            eta2b_t("eta2b_t", shape_r4, format_r4),
            rho1b_t("rho1b_t", shape_r2, format_r2),
            rho2b_t("rho2b_t", shape_r4, format_r4);
        
        
        f_t.pack();
        Gamma_t.pack();
        dGamma_t.pack();
        eta1b_t.pack();
        eta2b_t.pack();
        rho1b_t.pack();
        rho2b_t.pack();

        vector2tensor(f, f_t);
        vector2tensor(Gamma, Gamma_t);
        vector2tensor(dGamma, dGamma_t);
        vector2tensor(eta1b, eta1b_t);
        vector2tensor(eta2b, eta2b_t);
        vector2tensor(rho1b, rho1b_t);
        vector2tensor(rho2b, rho2b_t);


        IndexVar a,b,c,d,p,q;    
        sum1_0b() = occA_a(a,p)*occA_b(p,b)*eta1b_t(a,b)*f_t(a,b);
        sum2_0b1() = occD_a(a,p)*occD_b(p,b)*occD_c(c,q)*occD_d(q,d)*eta2b_t(a,b,c,d)*Gamma_t(c,d,a,b);
        sum2_0b2() = occD_a(a,p)*occD_b(p,b)*occD_c(c,q)*occD_d(q,d)*Gamma_t(a,b,c,d)*eta2b_t(c,d,a,b);
        sum3_0b() = dGamma_t(a,b,c,d)*rho2b_t(a,b,c,d);

        dE_t() = sum1_0b() + 0.25*(sum2_0b1() - sum2_0b2()) + 0.25*sum3_0b();

        dE_t.evaluate();
    
        dE = static_cast<const double*>(dE_t.getStorage().getValues().getData())[0];

        return dE;

    }

    //** up to 3B density
    double flow_0b(vector<double> &f, vector<double> &Gamma, 
                   vector<double> &eta1b, vector<double> &eta2b,
                   vector<double> &rho1b, vector<double> &rho2b, vector<double> &rho3b,
                   vector<double> &dGamma) {

        int numStates = (int)sqrt(f.size());
        std::vector<int> shape_r2 = {numStates,numStates};
        std::vector<int> shape_r3 = {numStates,numStates,numStates};
        std::vector<int> shape_r4 = {numStates,numStates,numStates,numStates};
        std::vector<int> shape_r6 = {numStates,numStates,numStates,numStates,numStates,numStates};
        Format format_r2({Dense,Dense});
        Format format_r3({Dense,Dense,Dense});
        Format format_r4({Dense,Dense,Dense,Dense});
        Format format_r6({Dense,Dense,Dense,Dense,Dense,Dense});

        double dE;
        Tensor<double> sum1_0b("sum1_0b"), sum2_0b1("sum2_0b1"), sum2_0b2("sum2_0b2"), sum3_0b("sum3_0b"), sum4_0b1("sum4_0b1"), sum4_0b2("sum4_0b2"), dE_t("dE");

        Tensor<double> 
            f_t("f_t", shape_r2, format_r2),
            Gamma_t("Gamma_t", shape_r4, format_r4),
            dGamma_t("dGamma_t", shape_r4, format_r4),
            eta1b_t("eta1b_t", shape_r2, format_r2),
            eta2b_t("eta2b_t", shape_r4, format_r4),
            rho1b_t("rho1b_t", shape_r2, format_r2),
            rho2b_t("rho2b_t", shape_r4, format_r4),
            rho3b_t("rho3b_t", shape_r6, format_r6);        
        
        f_t.pack();
        Gamma_t.pack();
        dGamma_t.pack();
        eta1b_t.pack();
        eta2b_t.pack();
        rho1b_t.pack();
        rho2b_t.pack();
        rho3b_t.pack();

        vector2tensor(f, f_t);
        vector2tensor(Gamma, Gamma_t);
        vector2tensor(dGamma, dGamma_t);
        vector2tensor(eta1b, eta1b_t);
        vector2tensor(eta2b, eta2b_t);
        vector2tensor(rho1b, rho1b_t);
        vector2tensor(rho2b, rho2b_t);
        vector2tensor(rho3b, rho3b_t);

        IndexVar a,b,c,d,k,l,m,p,q;    
        sum1_0b() = occA_a(a,p)*occA_b(p,b)*eta1b_t(a,b)*f_t(a,b);
        sum2_0b1() = occD_a(a,p)*occD_b(p,b)*occD_c(c,q)*occD_d(q,d)*eta2b_t(a,b,c,d)*Gamma_t(c,d,a,b);
        sum2_0b2() = occD_a(a,p)*occD_b(p,b)*occD_c(c,q)*occD_d(q,d)*Gamma_t(a,b,c,d)*eta2b_t(c,d,a,b);
        sum3_0b() = dGamma_t(a,b,c,d)*rho2b_t(a,b,c,d);
        sum4_0b1() = eta2b_t(a,b,c,d)*Gamma_t(k,l,a,m)*rho3b_t(b,k,l,c,d,m);
        sum4_0b2() = Gamma_t(a,b,c,d)*eta2b_t(k,l,a,m)*rho3b_t(b,k,l,c,d,m);

        dE_t() = sum1_0b() + 0.25*(sum2_0b1() - sum2_0b2()) + 0.25*sum3_0b() + 0.25*(sum4_0b1() - sum4_0b2());

        dE_t.evaluate();
    
        dE = static_cast<const double*>(dE_t.getStorage().getValues().getData())[0];

        return dE;

    }

    //**  up to 2B density
    vector<double> flow_1b(vector<double> &f, vector<double> &Gamma, 
                           vector<double> &eta1b, vector<double> &eta2b,
                           vector<double> &rho1b, vector<double> &rho2b) {

        int numStates = (int)sqrt(f.size());
        std::vector<int> shape_r2 = {numStates,numStates};
        std::vector<int> shape_r4 = {numStates,numStates,numStates,numStates};
        Format format_r2({Dense,Dense});
        Format format_r4({Dense,Dense,Dense,Dense});

        vector<double> df(f.size());

        Tensor<double> 
            sum1_1b(shape_r2, format_r2), 
            sum2_1b(shape_r2, format_r2), 
            sum3_1b(shape_r2, format_r2), 
            sum4_1b(shape_r2, format_r2), 
            sum5_1b(shape_r2, format_r2), 
            sum6_1b(shape_r2, format_r2), 
            sum7_1b(shape_r2, format_r2), 
            df_t(shape_r2, format_r2);

        Tensor<double> 
            f_t("f_t", shape_r2, format_r2),
            Gamma_t("Gamma_t", shape_r4, format_r4),
            eta1b_t("eta1b_t", shape_r2, format_r2),
            eta2b_t("eta2b_t", shape_r4, format_r4),
            rho1b_t("rho1b_t", shape_r2, format_r2),
            rho2b_t("rho2b_t", shape_r4, format_r4);
    
        f_t.pack();
        Gamma_t.pack();
        eta1b_t.pack();
        eta2b_t.pack();
        rho1b_t.pack();
        rho2b_t.pack();

        vector2tensor(f, f_t);
        vector2tensor(Gamma, Gamma_t);
        vector2tensor(eta1b, eta1b_t);
        vector2tensor(eta2b, eta2b_t);
        vector2tensor(rho1b, rho1b_t);
        vector2tensor(rho2b, rho2b_t);
        
        IndexVar i,j,k,l,a,b,c,d,e,p,q;

        sum1_1b(i,j) = eta1b_t(i,a)*f_t(a,j) - f_t(i,a)*eta1b_t(a,j); //enforcing anti-hermiticity?
        sum2_1b(i,j) = \
            occA_a(a,p)*occA_b(p,b)*eta1b_t(a,b)*Gamma_t(b,i,a,j) - \
            occA_a(a,p)*occA_b(p,b)*f_t(a,b)*eta2b_t(b,i,a,j);
        sum3_1b(i,j) = \
            occC2_b(a,p)*occC2_c(p,b,q)*occC2_a(q,c)*eta2b_t(i,a,b,c)*Gamma_t(b,c,j,a) - \
            occC2_b(a,p)*occC2_c(p,b,q)*occC2_a(q,c)*Gamma_t(i,a,b,c)*eta2b_t(b,c,j,a);
        sum4_1b(i,j) = \
            eta2b_t(i,a,b,c)*Gamma_t(d,e,j,a)*rho2b_t(d,e,b,c) - \
            Gamma_t(i,a,b,c)*eta2b_t(d,e,j,a)*rho2b_t(d,e,b,c);
        sum5_1b(i,j) = \
            eta2b_t(i,a,b,c)*Gamma_t(b,e,j,d)*rho2b_t(a,e,c,d) - \
            Gamma_t(i,a,b,c)*eta2b_t(b,e,j,d)*rho2b_t(a,e,c,d);        
        sum6_1b(i,j) = \
            eta2b_t(i,a,j,b)*Gamma_t(c,d,a,e)*rho2b_t(c,d,b,e) - \
            Gamma_t(i,a,j,b)*eta2b_t(c,d,a,e)*rho2b_t(c,d,b,e);        
        sum7_1b(i,j) = \
            eta2b_t(i,a,j,b)*Gamma_t(b,c,d,e)*rho2b_t(a,c,d,e) -    \
            Gamma_t(i,a,j,b)*eta2b_t(b,c,d,e)*rho2b_t(a,c,d,e);        

    
        df_t(i,j) = sum1_1b(i,j) + sum2_1b(i,j) + 0.5*sum3_1b(i,j) + 0.25*sum4_1b(i,j) + sum5_1b(i,j) - 0.5*sum6_1b(i,j) + 0.5*sum6_1b(i,j);
    
        df_t.evaluate();

        // for (int a = 0; a <numStates; a++) {
        //     for (int i = 0; i < numStates; i++) {
        //         if (df_t(a,i) > 1e10)
        //             std::cout << " BIG F VALUE " << df_t(a,i) << std::endl;
        //     }
        // }

        tensor2vector(df_t, df);

        return df;

    }

    //** up to 2B density
    vector<double> flow_2b(vector<double> &f, vector<double> &Gamma, 
                           vector<double> &eta1b, vector<double> &eta2b,
                           vector<double> &rho1b, vector<double> &rho2b) {
        int numStates = (int)sqrt(f.size());
        std::vector<int> shape_r2 = {numStates,numStates};
        std::vector<int> shape_r4 = {numStates,numStates,numStates,numStates};
        Format format_r2({Dense,Dense});
        Format format_r4({Dense,Dense,Dense,Dense});

        vector<double> dGamma(Gamma.size());

        Tensor<double> 
            sum1_2b(shape_r4, format_r4), 
            sum2_2b(shape_r4, format_r4), 
            sum3_2b(shape_r4, format_r4), 
            dGamma_t(shape_r4, format_r4);

        Tensor<double> 
            f_t("f_t", shape_r2, format_r2),
            Gamma_t("Gamma_t", shape_r4, format_r4),
            eta1b_t("eta1b_t", shape_r2, format_r2),
            eta2b_t("eta2b_t", shape_r4, format_r4),
            rho1b_t("rho1b_t", shape_r2, format_r2),
            rho2b_t("rho2b_t", shape_r4, format_r4);

        f_t.pack();
        Gamma_t.pack();
        eta1b_t.pack();
        eta2b_t.pack();

        vector2tensor(f, f_t);
        vector2tensor(Gamma, Gamma_t);
        vector2tensor(eta1b, eta1b_t);
        vector2tensor(eta2b, eta2b_t);

        IndexVar i,j,k,l,a,b,p;

        sum1_2b(i,j,k,l) = eta1b_t(i,a)*Gamma_t(a,j,k,l) - f_t(i,a)*eta2b_t(a,j,k,l) + eta1b_t(j,a)*Gamma_t(i,a,k,l) - f_t(j,a)*eta2b_t(i,a,k,l) \
            - eta1b_t(a,k)*Gamma_t(i,j,a,l) + f_t(a,k)*eta2b_t(i,j,a,l) - eta1b_t(a,l)*Gamma_t(i,j,k,a) - f_t(a,l)*eta2b_t(i,j,k,a);

        sum2_2b(i,j,k,l) = occB_a(a,p)*occB_b(p,b)*eta2b_t(i,j,a,b)*Gamma_t(a,b,k,l) - occB_a(a,p)*occB_b(p,b)*Gamma_t(i,j,a,b)*eta2b_t(a,b,k,l);

        sum3_2b(i,j,k,l) = occA_a(a,p)*occA_b(p,b)*eta2b_t(i,a,k,b)*Gamma_t(j,b,l,a) - occA_a(a,p)*occA_b(p,b)*eta2b_t(j,a,k,b)*Gamma_t(i,b,l,a) \
            - occA_a(a,p)*occA_b(p,b)*Gamma_t(i,a,k,b)*eta2b_t(j,b,l,a) + occA_a(a,p)*occA_b(p,b)*Gamma_t(j,a,k,b)*eta2b_t(i,b,l,a);

        dGamma_t(i,j,k,l) = sum1_2b(i,j,k,l) + 0.5*sum2_2b(i,j,k,l) - sum3_2b(i,j,k,l);
    
        dGamma_t.evaluate();

        tensor2vector(dGamma_t, dGamma);

        return dGamma;

    }


};
