#include <stdlib.h>
#include <stdio.h>

/* gcc -Ddgetrf=dgetrf_ -Ddgetri=dgetri_ m_inversion.c -lblas -llapack */

extern void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);
extern void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);

/* Allocate a double matrix. */

double ** alloc_dmatrix(int imax, int jmax){

    double ** out = NULL;

    out = (double**)calloc(imax,sizeof(double));

    if (out == NULL) {printf("ERROR: double matrix allocation.\n"); exit(1);}

    for (int i = 0; i<imax; i++){
        out[i] = (double*)calloc(jmax,sizeof(double));
        if (out[i] == NULL) {printf("ERROR: double matrix allocation.\n"); exit(1);}
    }

    /* Initialize the matrix. */

    for (int i = 0; i<imax; i++)
        for (int j = 0; j<jmax; j++) out[i][j] = 0.0;

    return out;
}

/* Maps an matrix to an array. Where tam is the dimension in i direction
 * following the C notation standard. Do not work with this guy thinking in
 * Fortran cache efficiency. If so, your code will work fast, but it will be
 * wrong. */

int ind2d(int i, int j, int tam){
    return (i)*(tam)+j; 
    
}

void free_dmatrix(double ** matrix, int imax, int jmax){

    for (int i = 0; i<imax; i++)
        free(matrix[i]);

    free(matrix);
}

void inv(double ** matrix, int N)
{
    int * IPIV = calloc((N+1),sizeof(int));
    int LWORK = N*N;
    double * WORK = calloc(LWORK,sizeof(double));
    int INFO;

    double A[N*N];

    for (int i = 0; i<N; i++)
        for (int j = 0; j<N; j++)
            A[(i)*(N)+j] = matrix[i][j];

    dgetrf_(&N,&N,A,&N,IPIV,&INFO);
    dgetri_(&N,A,&N,IPIV,WORK,&LWORK,&INFO);

    for (int i = 0; i<N; i++)
        for (int j = 0; j<N; j++)
            matrix[i][j] = A[(i)*(N)+j];

    free(IPIV);
    free(WORK);
}

int main(){

    double ** A = alloc_dmatrix(2,2);

    A[0][0] = 1.0; A[0][1] = 2.0;
    A[1][0] = 3.0; A[1][1] = 4.0;

    inv(A, 2);

    printf("\nOriginal Inversion results\n");

    printf("%lf %lf\n", A[0][0], A[0][1]);
    printf("%lf %lf\n", A[1][0], A[1][1]);

    free_dmatrix(A, 2, 2);

    return 0;
}
