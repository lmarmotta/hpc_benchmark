#include <stdlib.h>
#include <stdio.h>

/* Multiply two double squared matrices. */

void dsmm(double ** mA, double ** mB, double ** mC, int size){

    for (int i=0; i<size; i++)
        for (int j=0; j<size; j++) mC[i][j] = 0.0;

    for (int i=0; i<size; i++)
        for (int j=0; j<size; j++) 
            for (int k=0; k<size; k++)
                mC[i][j] = mC[i][j] + mA[i][k]*mB[k][j];
}

/* Matrix multiplication any size. Remember what I always forget:
 * [m x n] X [n x p] = [m x p] Where [a,b] are lines and columns. 
 *    A    x    B    =    C */ 

void dmm(double ** mA, double ** mB, double ** mC, int m, int n, int p){

    for (int i=0; i<m; i++)
        for (int j=0; j<p; j++) mC[i][j] = 0.0;

    for (int i = 0; i < m; i++)
        for (int j = 0; j < p; j++)
            for (int k = 0; k < n; k++)
                mC[i][j] = mC[i][j] + mA[i][k]*mB[k][j];
}

int main(){

    int size = 2;

    printf("Testing the square size matrix multiplication.\n");

    double ** A = (double**)calloc(size,sizeof(double*));
    for (int i = 0; i<size; i++) A[i] = (double*)calloc(size,sizeof(double));

    double ** B = (double**)calloc(size,sizeof(double*));
    for (int i = 0; i<size; i++) B[i] = (double*)calloc(size,sizeof(double));

    double ** C = (double**)calloc(size,sizeof(double*));
    for (int i = 0; i<size; i++) C[i] = (double*)calloc(size,sizeof(double));

    A[0][0] = 2.0;
    A[0][1] = 2.0;
    A[1][0] = 2.0;
    A[1][1] = 2.0;

    B[0][0] = 2.0;
    B[0][1] = 2.0;
    B[1][0] = 2.0;
    B[1][1] = 2.0;

    C[0][0] = 0.0;
    C[0][1] = 0.0;
    C[1][0] = 0.0;
    C[1][1] = 0.0;

    dsmm(A, B, C, 2);

    printf("%lf %lf\n",C[0][0], C[0][1]);
    printf("%lf %lf\n",C[1][0], C[1][1]);

    printf("Testing the any size matrix multiplication.\n");

    int m = 2;
    int n = 3;
    int p = 2;

    double ** mA = (double**)calloc(n,sizeof(double*));
    for (int i = 0; i<m; i++) mA[i] = (double*)calloc(size,sizeof(double));

    double ** mB = (double**)calloc(p,sizeof(double*));
    for (int i = 0; i<n; i++) mB[i] = (double*)calloc(size,sizeof(double));

    double ** mC = (double**)calloc(p,sizeof(double*));
    for (int i = 0; i<m; i++) mC[i] = (double*)calloc(size,sizeof(double));

    /* Matrix A. */

    mA[0][0] = 1.0;
    mA[0][1] = 2.0;
    mA[0][2] = 3.0;

    mA[1][0] = 4.0;
    mA[1][1] = 5.0;
    mA[1][2] = 6.0;

    /* Matrix B. */

    mB[0][0] = 7.0;
    mB[0][1] = 8.0;

    mB[1][0] = 9.0;
    mB[1][1] = 10.0;

    mB[2][0] = 11.0;
    mB[2][1] = 12.0;

    dmm(mA, mB, mC, m, n, p);

    printf("%lf %lf\n",mC[0][0], mC[0][1]);
    printf("%lf %lf\n",mC[1][0], mC[1][1]);

    return 0;
}
