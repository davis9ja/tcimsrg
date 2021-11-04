#ifndef IMSRG_UTILS_HPP_
#define IMSRG_UTILS_HPP_

#include <boost/numeric/ublas/vector.hpp>
#include "taco.h"
#include "occupation_factors.hpp"

void readOccTensors(std::string factor_path,
                    taco::Tensor<double> &occA_a, taco::Tensor<double> &occA_b, 
                    taco::Tensor<double> &occB_a, taco::Tensor<double> &occB_b,
                    taco::Tensor<double> &occC_a, taco::Tensor<double> &occC_b, 
                    taco::Tensor<double> &occC_c,
                    taco::Tensor<double> &occD_a, taco::Tensor<double> &occD_b,
                    taco::Tensor<double> &occD_c, taco::Tensor<double> &occD_d);

void vector2tensor(boost::numeric::ublas::vector<double> &input, taco::Tensor<double> &output);
void tensor2vector(taco::Tensor<double> &input, boost::numeric::ublas::vector<double> &output);


#endif
