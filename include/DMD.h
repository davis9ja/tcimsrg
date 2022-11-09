#ifndef DMD_H_
#define DMD_H_

#include <mkl.h>
#include <mkl_types.h>
#include <mkl_lapacke.h>
#include <stdio.h>
#include <complex.h>
#include <string.h>
#include <math.h>

//#define MKL_Complex16 double __complex__

typedef struct {
    double* phi;
    double* eig;
    double* b;
} DMD;

void compute_dmd(DMD* dmd, double* data, int nObs, int nVariables, int svdRank); 

/* void complex_divide(double* u, double* v, double x, double y); */
/* void readDataFromLog(FILE *rfp, MKL_Complex16* data, int nVariables, int nObs); */

#endif
