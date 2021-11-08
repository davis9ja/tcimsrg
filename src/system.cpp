#include "system.hpp"

using namespace taco;
using namespace boost::numeric::ublas;
//using namespace std;

System::System(int numStates,
               double &E, vector<double> &f, vector<double> &Gamma, vector<double> &W,
               White *white, Flow_IMSRG2 *flow, SystemObserver *observer) {

    int fSize = numStates*numStates;
    int GammaSize = fSize*fSize;

    this->numStates = numStates;
    this->E = E;
    this->f = f;
    this->Gamma = Gamma;
    this->W = W;
    this->white = white;
    this->flow = flow;
    this->observer = observer;
    
    //this->sys_vec = vector<double>(1+fSize+GammaSize);

    // this->sys_vec.reserve(1+fSize+GammaSize);
    
    // for(int i = 0; i < fSize+GammaSize+1; i++) {
    //     sys_vec.push_back(0.0);
    // }

    //system2vector();

    //printf("%8s\t%10s\t%10s\t%10s\n", "s", "E", "||eta1b||", "||eta2b||");
}

void System::system2vector(double &E, vector<double> &f, vector<double> &Gamma, state_type &x) {
    
    // state_type y;

    // double val = static_cast<const double*>(E.getStorage().getValues().getData())[0];
    // y.push_back(val);

    // for(int i = 0; i < f.getStorage().getValues().getSize(); i++)
    //     y.push_back(static_cast<const double*>(f.getStorage().getValues().getData())[i]);

    // for(int i = 0; i < Gamma.getStorage().getValues().getSize(); i++)
    //     y.push_back(static_cast<const double*>(Gamma.getStorage().getValues().getData())[i]);
    

    // return y;

    //double E_val, f_val, Gamma_val;
    // int fSize = f.getStorage().getValues().getSize();
    // int GammaSize = Gamma.getStorage().getValues().getSize();
    int fSize = f.size();
    int GammaSize = Gamma.size();
    //state_type x(1+fSize+GammaSize);

    //E_val = static_cast<const double*>(E.getStorage().getValues().getData())[0];
    x[0] = E;

    //#pragma omp parallel for
    for(int i = 0; i < fSize; i++) {
        //f_val = static_cast<const double*>(f.getStorage().getValues().getData())[i];
        x[i+1] = f[i];
    }
    //std::cout << "done 1b" << std::endl;
    //#pragma omp parallel for
    for(int i = 0; i < GammaSize; i++) {
        //Gamma_val = static_cast<const double*>(Gamma.getStorage().getValues().getData())[i];
        x[i+1+fSize] = Gamma[i];
    }
    //std::cout << "done 2b" << std::endl;

}

void System::vector2system(const state_type &x, int fSize, double &E, vector<double> &f, vector<double> &Gamma) {
    
    //double E_val, f_val, Gamma_val;
    // int fSize = f.getStorage().getValues().getSize();
    // int GammaSize = Gamma.getStorage().getValues().getSize();
    int GammaSize = fSize*fSize;
    //state_type x(1+fSize+GammaSize);

    //E_val = static_cast<const double*>(E.getStorage().getValues().getData())[0];
    E = x[0];

    //#pragma omp parallel for
    for(int i = 0; i < fSize; i++) {
        //f_val = static_cast<const double*>(f.getStorage().getValues().getData())[i];
        f[i] = x[i+1];
    }
    //std::cout << "done 1b" << std::endl;
    //#pragma omp parallel for
    for(int i = 0; i < GammaSize; i++) {
        //Gamma_val = static_cast<const double*>(Gamma.getStorage().getValues().getData())[i];
        Gamma[i] = x[i+1+fSize];
    }
    //std::cout << "done 2b" << std::endl;

}

void System::reinitSystem(state_type x) {
    
    this->E = x[0];

    // double* fArr = (double*)this->f.getStorage().getValues().getData();
    // double* GammaArr = (double*)this->Gamma.getStorage().getValues().getData();
    int fSize = numStates*numStates;
    int GammaSize = fSize*fSize;

    //#pragma omp parallel for
    for(int i = 1; i < fSize+1; i++)
        this->f[i-1] = x[i];
    //std::cout << "done 1b" << std::endl;
    //#pragma omp parallel for
    for(int i = fSize+1; i < GammaSize+1; i++)
        this->Gamma[i-fSize-1] = x[i];
    //std::cout << "done 2b" << std::endl;
}

void System::operator() (const state_type &x, state_type &dxdt, const double t) {
    //std::cout << "BEFORE, " << x[0] << std::endl;
    //reinitSystem(x);
    //std::cout << "AFTER, " << x[0] << std::endl;

    //std::cout << "SIZE OF X " << x.size() << std::endl;
    //std::cout << "SIZE OF E+F+GAMMA" << 1+f.size()+Gamma.size() << std::endl;

    //std::cout << "rinit" << std::endl;
    //std::cout << "x[0] " << x[0] << std::endl;
    
    // double E = this->E;
    // vector<double> f = this->f;
    // vector<double> Gamma = this->Gamma;

    //std::cout << "E, " << E << std::endl;
    //std::cout << "set up system variable" << std::endl;
    // calculate generator from f, Gamma

    // Take state vector and extract into E, f, Gamma pointers
    vector2system(x, f.size(), E, f, Gamma);

    // Compute generator from state
    vector<double> eta1b = (*white).compute_1b(f, Gamma, W);
    vector<double> eta2b = (*white).compute_2b(f, Gamma, W);

    // Flow E, f, Gamma
    dE = (*flow).flow_0b(f, Gamma, eta1b, eta2b);
    df = (*flow).flow_1b(f, Gamma, eta1b, eta2b);
    dGamma = (*flow).flow_2b(f, Gamma, eta1b, eta2b);
    
    // Set derivative
    system2vector(dE, df, dGamma, dxdt);

    double norm;

    norm = 0.0;
    for (int i = 0; i < eta1b.size(); i++)
        norm += eta1b[i]*eta1b[i];
    this->eta1b_norm = sqrt(norm);

    norm = 0.0;
    for (int i = 0; i < eta2b.size(); i++)
        norm += eta2b[i]*eta2b[i];
    this->eta2b_norm = sqrt(norm);


    (*observer).setEta2bNorm(this->eta2b_norm);
    // std::cout << eta2b_norm << std::endl;
    
    //std::cout << "eta2b, " << this->getEta2bNorm() << std::endl;

    //std::cout << "F AT THE END " << f[7*numStates+3] << std::endl;
    //std::cout << "de, " << dxdt[0] << std::endl;
    //std::cout << "syste2mvector" << std::endl;
    //reinitSystem();


    
    // for (int i = 0; i  < sys_vec.size(); i++)
    //     cout << sys_vec[i] << endl;

    //printf("%8s\t%10s\t%10s\t%10s\n", t, x[0], this->eta1b_norm, this->eta2b_norm);
}

