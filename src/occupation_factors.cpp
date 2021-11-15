#include "occupation_factors.hpp"
#include <boost/numeric/ublas/io.hpp>

using namespace taco;
using namespace boost::numeric::ublas;

OccupationFactors::OccupationFactors(int numStates, vector<double> ref)
{
    // this->n_holes = n_holes;
    // this->n_particles = n_particles;
    this->numStates = numStates;
    this->ref = ref;

    const char* path_char = &path[0];
    int check = mkdir(path_char, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    if (check == -1) {
        printf("Occ factor files probably exist -- not writing to file\n");
        return;
    }

    writeA();
    writeB();
    writeC();
    writeD();
}

std::string OccupationFactors::getPath() {
    return path;
}

void OccupationFactors::readOccTensors(std::string factor_path,
                                       Tensor<double> &occA_a, Tensor<double> &occA_b, 
                                       Tensor<double> &occB_a, Tensor<double> &occB_b,
                                       Tensor<double> &occC_a, Tensor<double> &occC_b, 
                                       Tensor<double> &occC_c,
                                       Tensor<double> &occD_a, Tensor<double> &occD_b,
                                       Tensor<double> &occD_c, Tensor<double> &occD_d) {
    
    Format format_r2 = Format({Dense, Dense});
    Format format_r3 = Format({Dense, Dense, Dense});
    
    occA_a = read(factor_path+"occA_a.tns", format_r2);
    occA_b = read(factor_path+"occA_b.tns", format_r2);

    occB_a = read(factor_path+"occB_a.tns", format_r2);
    occB_b = read(factor_path+"occB_b.tns", format_r2);
    
    occC_a = read(factor_path+"occC_a.tns", format_r2);
    occC_b = read(factor_path+"occC_b.tns", format_r3);
    occC_c = read(factor_path+"occC_c.tns", format_r2);

    occD_a = read(factor_path+"occD_a.tns", format_r2);
    occD_b = read(factor_path+"occD_b.tns", format_r2);
    occD_c = read(factor_path+"occD_c.tns", format_r2);
    occD_d = read(factor_path+"occD_d.tns", format_r2);
}

void OccupationFactors::contractOccTensors(std::string factor_path,
                                           vector<double> &occA, 
                                           vector<double> &occB, 
                                           vector<double> &occC,
                                           vector<double> &occD) {
    //int numStates = n_holes+n_particles;

    Tensor<double> 
        occA_a, occA_b, 
        occB_a, occB_b,
        occC_a, occC_b, occC_c,
        occD_a, occD_b, occD_c, occD_d;

    // TRANSFORM OCCUPATION TENSORS
    std::vector<int> shape_r2 = {numStates,numStates};
    std::vector<int> shape_r3 = {numStates,numStates,numStates};
    std::vector<int> shape_r4 = {numStates,numStates,numStates,numStates};
    Format format_r2({Dense,Dense});
    Format format_r3({Dense,Dense,Dense});
    Format format_r4({Dense,Dense,Dense,Dense});

    Tensor<double> occA_t(shape_r2, format_r2);
    Tensor<double> occB_t(shape_r2, format_r2);
    Tensor<double> occC_t(shape_r3, format_r3);
    Tensor<double> occD_t(shape_r4, format_r4);
    // vector<double> occA(numStates*numStates);
    // vector<double> occD(numStates*numStates*numStates*numStates);

    readOccTensors(factor_path, 
                   occA_a, occA_b, 
                   occB_a, occB_b,
                   occC_a, occC_b, occC_c,
                   occD_a, occD_b, occD_c, occD_d);

    IndexVar a,b,c,d,p,q;
    occA_t(a,b) = occA_a(a,p)*occA_b(p,b);
    occB_t(a,b) = occB_a(a,p)*occB_b(p,b);
    occC_t(a,b,c) = occC_a(a,p)*occC_b(p,b,q)*occC_c(q,c);
    occD_t(a,b,c,d) = occD_a(a,p)*occD_b(p,b)*occD_c(c,q)*occD_d(q,d);

    occA_t.evaluate();
    occB_t.evaluate();
    occC_t.evaluate();
    occD_t.evaluate();

    tensor2vector(occA_t, occA);
    tensor2vector(occB_t, occB);
    tensor2vector(occC_t, occC);
    tensor2vector(occD_t, occD);
    
}

void OccupationFactors::writeA() {
    //int numStates = n_holes+n_particles;

    Tensor<double> occA_a({numStates,2}, Format({Dense,Dense}));
    Tensor<double> occA_b({2,numStates}, Format({Dense,Dense}));

    for (int a = 0; a < numStates; a++) {
        double val = ref[a];
        occA_a.insert({a,0}, val);
        occA_a.insert({a,1}, 1.0);

        occA_b.insert({0,a}, 1.0);
        occA_b.insert({1,a}, -val);
    }

    occA_a.pack();
    occA_b.pack();

    write(path+"occA_a.tns", occA_a);
    write(path+"occA_b.tns", occA_b);

}

void OccupationFactors::writeB() {
    //int numStates = n_holes+n_particles;

    Tensor<double> occB_a({numStates,2}, Format({Dense,Dense}));
    Tensor<double> occB_b({2,numStates}, Format({Dense,Dense}));

    for (int a = 0; a < numStates; a++) {
        double val = ref[a];
        occB_a.insert({a,0}, 1-val);
        occB_a.insert({a,1}, 1.0);
        
        occB_b.insert({0,a}, 1.0);
        occB_b.insert({1,a}, -val);
    }

    occB_a.pack();
    occB_b.pack();
       
    write(path+"occB_a.tns", occB_a);
    write(path+"occB_b.tns", occB_b);    
}

void OccupationFactors::writeC() {
    //int numStates = n_holes+n_particles;

    Tensor<double> occC_a({numStates,4}, Format({Dense,Dense}));
    Tensor<double> occC_b({4,numStates,4}, Format({Dense,Dense,Dense}));
    Tensor<double> occC_c({4, numStates}, Format({Dense,Dense}));

    for (int a = 0; a < numStates; a++) {
        double val = ref[a];
        occC_a.insert({a,0}, val);
        occC_a.insert({a,1}, 1.0);
        occC_a.insert({a,2}, 1.0);
        occC_a.insert({a,3}, val);

        occC_b.insert({0,a,0}, val);
        occC_b.insert({1,a,1}, 1.0);
        occC_b.insert({2,a,2}, -val);
        occC_b.insert({3,a,3}, 1.0);

        occC_c.insert({0,a}, 1.0);
        occC_c.insert({1,a}, val);
        occC_c.insert({2,a}, val);
        occC_c.insert({3,a}, val);        
    }

    occC_a.pack();
    occC_b.pack();
    occC_c.pack();

    write(path+"occC_a.tns", occC_a);
    write(path+"occC_b.tns", occC_b);
    write(path+"occC_c.tns", occC_c);    
}

void OccupationFactors::writeD() {
    //int numStates = n_holes+n_particles;

    Tensor<double> occD_a({numStates,2}, Format({Dense,Dense}));
    Tensor<double> occD_b({2,numStates}, Format({Dense,Dense}));
    Tensor<double> occD_c({numStates,2}, Format({Dense,Dense}));
    Tensor<double> occD_d({2,numStates}, Format({Dense,Dense}));

    for (int a = 0; a < numStates; a++) {
        double val = ref[a];

        occD_a.insert({a,0}, val);
        occD_b.insert({0,a}, val);
        occD_c.insert({a,0}, 1-val);
        occD_d.insert({0,a}, 1-val);
    }

    occD_a.pack();
    occD_b.pack();
    occD_c.pack();
    occD_d.pack();

    write(path+"occD_a.tns", occD_a);
    write(path+"occD_b.tns", occD_b);
    write(path+"occD_c.tns", occD_c);
    write(path+"occD_d.tns", occD_d);        
}

