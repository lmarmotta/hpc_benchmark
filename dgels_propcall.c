#include <stdlib.h>
#include <stdio.h>

#define M 6
#define N 4
#define NRHS 2
#define LDA M
#define LDB M

extern void dgels( char* trans, int* m, int* n, int* nrhs, double* a, int* lda,
                double* b, int* ldb, double* work, int* lwork, int* info );

int main(){

    double a[LDA*N] = {
        1.44, -9.96, -7.55,  8.34,  7.08, -5.45,
       -7.84, -0.28,  3.24,  8.09,  2.52, -5.70,
       -4.39, -3.24,  6.27,  5.28,  0.74, -1.19,
        4.53,  3.83, -6.64,  2.06, -2.47,  4.70
    };

    double b[LDB*NRHS] = {
        8.58,  8.26,  8.48, -5.28,  5.72,  8.93,
        9.35, -4.43, -0.70, -0.26, -7.36, -2.52
    };

    int m = M, n = N, nrhs = NRHS, lda = LDA, ldb = LDB, info, lwork;

    double wkopt;
    double* work;

    lwork = -1;

    dgels( "No transpose", &m, &n, &nrhs, a, &lda, b, &ldb, &wkopt, &lwork, &info );

    lwork = (int)wkopt;
    work = (double*)malloc( lwork*sizeof(double) );

    dgels( "No transpose", &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, &info );

    free(work);

    return 0;

}
