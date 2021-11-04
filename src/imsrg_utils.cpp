#include "imsrg_utils.hpp"

using namespace taco;
using namespace boost::numeric::ublas;

void readOccTensors(std::string factor_path,
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

void vector2tensor(vector<double> &input, Tensor<double> &output) {
        double* output_arr = (double*)output.getStorage().getValues().getData();
        for (int i = 0; i < input.size(); i++)
            output_arr[i] = input[i];
}

void tensor2vector(Tensor<double> &input, vector<double> &output) {
        double* input_arr = (double*)input.getStorage().getValues().getData();
        for (int i = 0; i < output.size(); i++)
            output[i] = input_arr[i];
}
