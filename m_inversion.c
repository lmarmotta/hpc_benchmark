#include <stdlib.h>
#include <stdio.h>

/* gcc -Ddgetrf=dgetrf_ -Ddgetri=dgetri_ m_inversion.c -lblas -llapack */

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

int main(){

    double A [2*2] = {
        1,2,
        3,4
    };

    inv(A, 2);

    printf("\nOriginal Inversion results\n");

    printf("%f %f\n", A[0], A[1]);
    printf("%f %f\n", A[2], A[3]);

    A[0] = 0.0, A[1] = 0.0;
    A[2] = 0.0, A[3] = 0.0;

    printf("\nCheck the index\n\n");

    for (int i = 0; i<2; i++)
        for (int j = 0; j<2; j++)
            printf("ind2d(%d,%d) = %d\n",i,j,ind2d(i,j,2));

    A[ind2d(0,0,2)] = 1.0;
    A[ind2d(0,1,2)] = 2.0;
    A[ind2d(1,0,2)] = 3.0;
    A[ind2d(1,1,2)] = 4.0;

    printf("\nInd2d inversion results.\n");

    inv(A, 2);

    printf("%f %f\n", A[0], A[1]);
    printf("%f %f\n", A[2], A[3]);

    return 0;
}
