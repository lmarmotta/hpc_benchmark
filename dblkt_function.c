#include <stdlib.h>
#include <stdio.h>

extern void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);
extern void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);

int ind2d(int i, int j, int tam){
    return (i)*(tam)+j; 
}

void inv(double * A, int N)
{
    int * IPIV = calloc((N+1),sizeof(int));
    int LWORK = N*N;
    double * WORK = calloc(LWORK,sizeof(double));
    int INFO;

    dgetrf_(&N,&N,A,&N,IPIV,&INFO);
    dgetri_(&N,A,&N,IPIV,WORK,&LWORK,&INFO);

    free(IPIV);
    free(WORK);
}

double ** alloc_double_matrix(int imax, int jmax){

    double ** matrix = NULL;

    matrix = (double**)calloc(imax, sizeof(double*));

    if (matrix == NULL) {
        printf("ERROR: Memory Allocation\n"); exit(1); 
    }

    for (int i = 0; i<imax; i++){

        matrix[i] = (double*)calloc(jmax, sizeof(double));

        if (matrix[i] == NULL) {
            printf("ERROR: Memory Allocation\n"); exit(1); 
        }
    }

    return matrix;
}

void blk_tri(double *** lower, double *** main, double *** upper, int size_m, int num_m, double ** XB, double *X){

}

int main(){

    double ld[2][2][3];
    double lm[2][2][3];
    double lu[2][2][3];
    double xb[2][3];

    /* Fill the test matrices. */

    for (int i = 0; i<3; i++){

        /* Fill the lower diagonal. */

        ld[0][0][i] = 4.0;
        ld[0][1][i] = 3.0;
        ld[1][0][i] = 3.0;
        ld[1][1][i] = 2.0;

        /* Fill the main matrix. */

        lm[0][0][i] = 8.0;
        lm[0][1][i] = 6.0;
        lm[1][0][i] = 6.0;
        lm[1][1][i] = 4.0;

        /* Fill the upper matrix. */

        lm[0][0][i] = 16.0;
        lm[0][1][i] = 9.0;
        lm[1][0][i] = 9.0;
        lm[1][1][i] = 6.0;

    }

    /* Fill the solution vector. */

    for (int i = 0; i<2; i++)
        for (int j = 0; j<2; j++)
            xb[i][j] = 3.0;

    return 0;
}
