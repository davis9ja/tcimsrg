/************************************************************************/
/* Calculate the Dynamic Mode Decomposition using LAPACK/BLAS routines, */
/* implemented in the Intel MKL.                                        */
/*                                                                      */
/*     Author: Jacob Davison                                            */
/*     Date:   02/17/2022                                               */
/************************************************************************/

#include "DMD.h"

/***********************************************************************/
/* Divide complex numbers (u + vi)/(x + yi) and write output into u,v. */
/***********************************************************************/
void print_matrix(double* mat, int rows, int columns) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < columns; j++) {
            //printf("%8.3f+%.2fi  ", mat[i*columns+j].real, mat[i*columns+j].imag);//creal(mat[i*columns+j]), cimag(mat[i*columns+j]));
            printf("%8.3f  ", mat[i*columns+j]);
        }
        printf("\n\n");
    }
}

void print_vector(double* vec, int dim) {
    for (int i = 0; i < dim; i++)
        printf("%8.3f  ", vec[i]);
    printf("\n");
}

void readDataFromLog(char* logname, double* data, int nVariables, int nObs) {
    
    FILE *rfp = fopen(logname, "r");



    fclose(rfp);
}


/**************************************************************************/
/* Compute the Dynamic Mode Decomposition for data matrix collected from  */
/* IMSRG flow.                                                            */
/**************************************************************************/
void compute_dmd(DMD* dmd, double* data, int nObs, int nVariables, int svdRank) {
    
    int dataBytes = sizeof(double)*nVariables*nObs;
    printf("\tDATA MATRIX SIZE = %6.4f GB\n", (double)dataBytes/1024/1024/1024);
    //printf("Initialize variables for DMD\n");
    
    MKL_INT m = nVariables, n = nObs-1, lda = nVariables, ldu = nVariables, ldvt = nObs-1, r, info;
    double alpha = 1.0, beta = 0.0;
    double cx, cy, cu, cv; // for complex division

    double superb[nObs-2];
    double s[n],  vt[n*n], u[1]; //u[m*n], vt[n*n];

    //printf("Malloc the data matrices\n");    
    double X[m*n]; //= malloc(m*n*2*sizeof(double));
    double Xp[m*n]; //= malloc(m*n*2*sizeof(double));
    double X0[m]; //= malloc(m*2*sizeof(double));



    //printf("Assigning data to shifted matrices\n");

    // Assign data to shifted matrices
    for (int i = 0; i < nVariables; i++) {
        for (int j = 0; j < nObs-1; j++)
            X[i*(nObs-1) + j] = data[i*(nObs) + j];

        for (int k = 0; k < nObs-1; k++)
            Xp[i*(nObs-1) + k] = data[i*(nObs) + (k+1)];

        X0[i] = data[i*(nObs) + 0];
    }

    /* double checkX, checkD; */
    /* for (int i = 0; i < nVariables; i++) { */
    /*     checkX = Xp[i*(nObs-1) + 3].real; */
    /*     checkD = data[i*nObs + 4].real; */
        
    /*     if (abs(checkX - checkD) > 1e-6) */
    /*         printf("caught data inconsistency in col 3., index %i\n", i); */
    /* } */
        

    /* printf("Data first five rows\n"); */
    /* print_matrix(data, 5, n+1); */

    /* printf("X first five rows\n"); */
    /* print_matrix(X, 5, n); */

    /* printf("Xp first five rows\n"); */
    /* print_matrix(Xp, 5, n); */

    /* printf("SVD the data matrix X\n"); */

    // SVD data matrix X (OVERWRITE LEFT SINGULAR VECTORS INTO X)
    info = LAPACKE_dgesdd(LAPACK_ROW_MAJOR, 'O', m, n, X, n, s, u, n, vt, n);
    //info = LAPACKE_dgesvd(LAPACK_ROW_MAJOR, 'S', 'S', m, n, X, n, s, u, n, vt, n, superb);

    printf("\tDGESDD complete with code %i\n", info);

    /* double tmp = 0.; */
    /* for (int i = 0; i < n; i++) { */
    /*     tmp = fabs(s[i]); */
    /*     if (tmp > 1e-3) { */
    /*         r = i+1; */
    /*     } */
    /* } */

    //*svdRank = r;

    r = svdRank;
    printf("\tTRUNCATING AT %i SINGULAR VALUE\n", r);
    
    // Last three terms in the sum are from malloc, which will be freed OUTSIDE of the function
    long totalBytes = sizeof(double)*( nObs-2 + n + 1 + n*n + m*n + m*n + m + r + m*r + r*n + r*n + r*r + r + r + 1 + r*r + r + m*r + r + r );

    printf("\tTOTAL SIZE TO ALLOCATE = %6.4f GB\n", (double)totalBytes/1024/1024/1024);
    
    (*dmd).phi = (double*)malloc(nVariables*r*sizeof(double));
    (*dmd).b = (double*)malloc(r*sizeof(double));
    (*dmd).eig = (double*)malloc(r*sizeof(double));


    double s_trunc[r], u_trunc[m*r], vt_trunc[r*n]; 
    double A_1[r*n], A_2[r*r];
    //MKL_Complex16 A_1[n*m], A_2[n*n];
    //double w[r];
    double wr[r], wi[r];
    double vl[1];
    double vr[r*r];

    //double D[r], E[r-1], tau[r-1], work_dsytrd[r*r];
    //double work_dsteqr[r*r];
    
    double b[r];//phi[m*r], b[r];

    for (int i = 0; i < m*r; i++)
        u_trunc[i] = 0.0;

    //printf("Truncate the SVD result\n");
    // Truncate single vectors at desired rank
    //r = min(r, nObs-1);

    if (r > nObs-1)
        r = nObs-1;
    
    for (int i = 0; i < r; i++)
        s_trunc[i] = s[i];

    for (int i = 0; i < nVariables; i++)
        for (int j = 0; j < r; j++)
            u_trunc[i*r + j] = X[i*n + j];

    for (int i = 0; i < r; i++)
        for (int j = 0; j < (nObs-1); j++)
            vt_trunc[i*(nObs-1) + j] = vt[i*(nObs-1) + j];


    //s_trunc[r-1] = 0.000244;

    for (int i = 0; i < r; i++)
        printf("%0.6f, ", s_trunc[i]);
    printf("\n");


    /* printf("vt full matrix\n"); */
    /* print_matrix(vt, n, n); */

    /* printf("Singular values up to truncation\n"); */
    /* print_vector(s, n); */

    /* printf("Right singular vectors up to truncation\n"); */
    /* print_matrix(vt_trunc, r, n); */

    /* printf("Calculate reduced DMD operator A\n"); */

    // Calculate reduced DMD operator via Xp projection into truncated space
    // A = U*.T @ Xp @ Vt*.T @ 1/s
    // Use CBLAS zgemm
    // *** IMPORTANT: when CblasConjTrans, use OTHER dimension for LD

    // A_1 = U*.T @ Xp = (r x m)(m x n) = (r x n)
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, r, n, m, alpha, u_trunc, r, Xp, n, beta, A_1, n);


    // A_2 = A_1 @ Vt*.T = (r x n)(n x r) = (r x r)
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans, r, r, n, alpha, A_1, n, vt_trunc, n, beta, A_2, r);

    // multiply each row j by 1/s(j,:)
    for (int i = 0; i < r; i++)
        for (int j = 0; j < r; j++) {
            /* cu = A_2[i*r + j].real; */
            /* cv = A_2[i*r + j].imag; */
            /* cx = s_trunc[j]; */
            /* cy = 0.0; */

            /* cu = (cu*cx + cv*cy)/(cx*cx + cy*cy); */
            /* cv = (cv*cx - cu*cy)/(cx*cx + cy*cy); */

            /* A_2[i*r + j].real = cu; */
            /* A_2[i*r + j].imag = cv; */

            A_2[i*r + j] = A_2[i*r + j]/s_trunc[j];
        }

    // Get eigenvalues from A (OVERWRITE EIGENVECTORS INTO A_2)
    info = LAPACKE_dgeev(LAPACK_ROW_MAJOR, 'N', 'V', r, A_2, r, wr, wi, vl, r, A_2, r);
    printf("\tDGEEV complete with code %i\n", info);

    // Compute Phi (OVERWRITE INTO u_trunc)
    // phi = Xp @ Vt*.T @ 1/s @ v @ 1/w

    // phi = Xp @ Vt*.T = (m x n)(n x r) = (m x r)
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans, m, r, n, alpha, Xp, n, vt_trunc, n, beta, u_trunc, r);
    
    // multply each row j by 1/s(j,:)
    for (int i = 0; i < nVariables; i++)
        for (int j = 0; j < r; j++) {
            /* cu = phi[i*r + j].real; */
            /* cv = phi[i*r + j].imag; */
            /* cx = s_trunc[j]; */
            /* cy = 0.0; */

            /* cu = (cu*cx + cv*cy)/(cx*cx + cy*cy); */
            /* cv = (cv*cx - cu*cy)/(cx*cx + cy*cy); */

            /* phi[i*r + j].real = cu; */
            /* phi[i*r + j].imag = cv; */
            
            u_trunc[i*r + j] = u_trunc[i*r + j]/s_trunc[j];

        }

    // phi = phi @ v = (m x r)(r x r) = (m x r)
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, r, r, alpha, u_trunc, r, A_2, r, beta, u_trunc, r);

    // multiply each row j by 1/w(j,:)
    for (int i = 0; i < nVariables; i++) {
        for (int j = 0; j < r; j++) {
            /* cu = phi[i*r + j].real; */
            /* cv = phi[i*r + j].imag; */
            /* cx = wr[j];//.real; */
            /* cy = wi[j];//.imag; */

            /* cu = (cu*cx + cv*cy)/(cx*cx + cy*cy); */
            /* cv = (cv*cx - cu*cy)/(cx*cx + cy*cy); */

            /* phi[i*r + j].real = cu; */
            /* phi[i*r + j].imag = cv; */

            //u_trunc[i*r + j] = u_trunc[i*r + j]/wr[j];
            
            (*dmd).phi[i*r + j] = u_trunc[i*r + j]/wr[j];
        }
    }
    
    printf("u_trunc, ");
    for (int i = 0; i < r; i++)
        printf("%0.6f, ", u_trunc[i]);
    printf("\n");

    // Set the DMD struct vals
    // (*dmd).phi = phi, (*dmd).sigma = sigma; (*dmd).b = b;
    /* for (int i = 0; i < nVariables; i++) */
    /*     for (int j = 0; j < r; j++) */
    /*         (*dmd).phi[i*r + j] = u_trunc[i*r + j]; */
    //memcpy((*dmd).phi, u_trunc, sizeof((*dmd).phi));

    printf("SIZE OF WR %i\n", (int)(sizeof(wr)/sizeof(double)));

    for (int i = 0; i < r; i++)
        (*dmd).eig[i] = wr[i];

    for (int i = 0; i < r; i++)
        printf("%0.6f, ", wi[i]);
    printf("\n");

    //memcpy((*dmd).eig, wr, sizeof((*dmd).eig));

    // Compute DMD mode amplitudes via least squares
    // b = LSTSQ(Phi, H0)
    info = LAPACKE_dgels(LAPACK_ROW_MAJOR, 'N', m, r, 1, u_trunc, r, X0, 1);
    printf("\tDGELS complete with code %i\n", info);

    for (int i = 0; i < r; i++)
        (*dmd).b[i] = X0[i];
    //memcpy((*dmd).b, X0, sizeof((*dmd).b));
    

}
