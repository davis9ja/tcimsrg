#include "system.hpp"

using namespace taco;
using namespace boost::numeric::ublas;
//using namespace std;

System::System(){}

System::System(int numStates, vector<double> &rho1b, vector<double> &rho2b,
               double &E, vector<double> &f, vector<double> &Gamma, vector<double> &W,
               Generator *generator, Flow_IMSRG2 *flow, 
               std::ofstream *out_file_vac_in, std::ofstream *out_file_imsrg_in,
               int reference_type) {

    int fSize = numStates*numStates;
    int GammaSize = fSize*fSize;

    this->numStates = numStates;
    this->rho1b = rho1b;
    this->rho2b = rho2b;
    this->E = E;
    this->f = f;
    this->Gamma = Gamma;
    this->W = W;
    this->generator = generator;
    this->flow = flow;
    out_file_vac = out_file_vac_in;
    out_file_imsrg = out_file_imsrg_in;
    this->reference_type = reference_type;
}

// System::~System() {

//     int size1b = numStates*numStates;
//     int size2b = numStates*numStates*numStates*numStates;
//     double h0b;
//     vector<double> h1b(size1b);
//     vector<double> h2b(size2b);

//     double E;
//     vector<double> f(size1b);
//     vector<double> Gamma(size2b);

    

//     if ((*out_file_vac).is_open() && reference_type==1) {
//         for (int i = 0; i < data_log.size(); i++) {
//             state_type x = data_log[i];
            
//             E = x[0];
            
//             for (int i = 1; i < size1b+1; i++)
//                 f[i-1] = x[i];

//             for (int i = 1+size1b; i < size2b+1; i++)
//                 Gamma[i-size1b-1] = x[i];

//             for (int i = 0; i < Gamma.size(); i++)
//                 h2b[i] = Gamma[i];
            
//             for (int p = 0; p < numStates; p++) {
//                 for (int q = 0; q < numStates; q++) {

//                     double sum = 0.0;
//                     for (int i = 0; i < numStates; i++)
//                         for (int j = 0; j < numStates; j++)
//                             sum += Gamma[p*numStates*numStates*numStates + i*numStates*numStates + q*numStates + j]*rho1b[i*numStates+j];
                    
//                     h1b[p*numStates+q] = f[p*numStates+q] - sum;
//                 }
//             }

//             double h1b_ij = 0.0;
//             for (int i = 0; i < numStates; i++)
//                 for (int j = 0; j < numStates; j++)
//                     h1b_ij += h1b[i*numStates+j]*rho1b[i*numStates+j];

//             double h2b_ijkl = 0.0;
//             for (int i = 0; i < numStates; i++)
//                 for (int j = 0; j < numStates; j++)
//                     for (int k = 0; k < numStates; k++)
//                         for (int l = 0; l < numStates; l++)
//                             h2b_ijkl += h2b[i*numStates*numStates*numStates + j*numStates*numStates + k*numStates + l]*rho2b[i*numStates*numStates*numStates + j*numStates*numStates + k*numStates + l];

//             h0b = E - h1b_ij - 0.25*h2b_ijkl;
            
//             *out_file_vac << h0b << ',';
//             for (int i = 0; i < h1b.size(); i++)
//                 *out_file_vac << h1b[i] << ',';
//             for (int i = 0; i < h2b.size()-1; i++)
//                 *out_file_vac << h2b[i] << ',';
//             *out_file_vac << h2b[h2b.size()-1];
//             *out_file_vac << std::fixed << std::setprecision(13) << "\n";

//         }
//     }

//     // if((*out_file_imsrg).is_open()) {

//     //     for (int i = 0; i < data_log.size(); i++) {
//     //         state_type x = data_log[i];

//     //         // for (int j = 0; j < x.size(); j++)
//     //         //     out_file_imsrg->write((char *) &x[j], sizeof(double));

//     //         // for (int j = 0; j < x.size(); j++)
//     //         //     *out_file_imsrg << x[j];

//     //         for (int j = 0; j < x.size()-1; j++)
//     //             *out_file_imsrg << x[j] << ',';
//     //         *out_file_imsrg << x[x.size()-1];
//     //         *out_file_imsrg << std::fixed << std::setprecision(13) << "\n";
//     //     }
//     // }
    
// }

void System::system2vector(double &E, vector<double> &f, vector<double> &Gamma, state_type &x) {

    int fSize = f.size();
    int GammaSize = Gamma.size();

    x[0] = E;

    //#pragma omp parallel for
    for(int i = 0; i < fSize; i++) {
        x[i+1] = f[i];
    }
    //std::cout << "done 1b" << std::endl;
    //#pragma omp parallel for
    for(int i = 0; i < GammaSize; i++) {
        x[i+1+fSize] = Gamma[i];
    }

}

void System::vector2system(const state_type &x, int fSize, double &E, vector<double> &f, vector<double> &Gamma) {
    
    int GammaSize = fSize*fSize;

    E = x[0];

    //#pragma omp parallel for
    for(int i = 0; i < fSize; i++) {
        f[i] = x[i+1];
    }

    //#pragma omp parallel for
    for(int i = 0; i < GammaSize; i++) {
        Gamma[i] = x[i+1+fSize];
    }
}


void System::operator() (const state_type &x, state_type &dxdt, const double t) {
    // Take state vector and extract into E, f, Gamma pointers
    vector2system(x, f.size(), E, f, Gamma);

    // Compute generator from state
    vector<double> eta1b = (*generator).compute_1b(f, Gamma, W);
    //std::cout << eta1b << std::endl;
    vector<double> eta2b = (*generator).compute_2b(f, Gamma, W);

    // Flow E, f, Gamma
    switch (reference_type) {
    case 0:
        dE = (*flow).flow_0b(f, Gamma, eta1b, eta2b);
        df = (*flow).flow_1b(f, Gamma, eta1b, eta2b);
        dGamma = (*flow).flow_2b(f, Gamma, eta1b, eta2b);
        break;
    case 1:        
        dGamma = (*flow).flow_2b(f, Gamma, eta1b, eta2b, rho1b, rho2b);
        df = (*flow).flow_1b(f, Gamma, eta1b, eta2b, rho1b, rho2b);
        dE = (*flow).flow_0b(f, Gamma, eta1b, eta2b, rho1b, rho2b, dGamma);
        break;
    default:
        dE = (*flow).flow_0b(f, Gamma, eta1b, eta2b);
        df = (*flow).flow_1b(f, Gamma, eta1b, eta2b);
        dGamma = (*flow).flow_2b(f, Gamma, eta1b, eta2b);
        
    }
    // Set derivative
    system2vector(dE, df, dGamma, dxdt);

}

void System::operator()( const state_type &x , double t )
{
    double norm;
    vector<double> eta1b, eta2b;

    if (t == 0.0)
        printf("\n%-4s\t%-8s\t%-8s\t%-8s\n", "t", "E", "||eta1b||", "||eta2b||");

    vector2system(x, f.size(), E, f, Gamma);

    // Compute generator from state
    eta1b = (*generator).compute_1b(f, Gamma, W);
    eta2b = (*generator).compute_2b(f, Gamma, W);

    norm = 0.0;
    for (int i = 0; i < eta1b.size(); i++)
        norm += eta1b[i]*eta1b[i];
    eta1b_norm = sqrt(norm);

    
    norm = 0.0;
    for (int i = 0; i < eta2b.size(); i++)
        norm += eta2b[i]*eta2b[i];
    eta2b_norm = sqrt(norm);


    printf("%0.4f\t%0.8f\t%0.8f\t%0.8f\n", t, x[0], eta1b_norm, eta2b_norm);

    data_log.push_back(x);
    
}
