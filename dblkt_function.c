#include <stdlib.h>
#include <stdio.h>

/* gcc -O0 -ggdb -Wall -std=c99 -Ddgetrf=dgetrf_ -Ddgetri=dgetri_ dblkt_function.c -lblas -llapack */

/* To do...
 * [ ] - Make an cube allocator.
 * [ ] - Test the argument passing before implementation.
 * [ ] - Copy solver function. */

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

/* Alloc a double cube. */

double *** alloc_dcube(int imax, int jmax, int kmax){

    double *** out = NULL;

    /* Allocate the first line. */

    out = (double***)calloc(imax,sizeof(double**));

    /* Check the allocation pointer position. */

    if (out == NULL) {printf("ERROR: double matrix allocation.\n"); exit(1);}

    /* Allocate the matrix. */

    for (int i = 0; i<imax; i++){
        out[i] = (double**)calloc(jmax,sizeof(double*));
        if (out[i] == NULL) {printf("ERROR: double matrix allocation.\n"); exit(1);}
    }

    /* Allocate the block. */

    for (int i = 0; i<imax; i++){
        for (int j = 0; j<jmax; j++){
            out[i][j] = (double*)calloc(kmax,sizeof(double));
            if (out[i][j] == NULL) {printf("ERROR: double matrix allocation.\n"); exit(1);}
        }
    }

    /* Initialize the cube. */

    for (int i = 0; i<imax; i++)
        for (int j = 0; j<jmax; j++) 
            for (int k = 0; k<kmax; k++) out[i][j][k] = 0.0;

    return out;
}

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

void free_dmatrix(double ** matrix, int imax, int jmax){

    for (int i = 0; i<imax; i++)
        free(matrix[i]);

    free(matrix);
}

void free_dcube(double *** cube, int imax, int jmax, int kmax){

    for (int i = 0; i<imax; i++)
        for (int j = 0; j<jmax; j++) free(cube[i][j]);
            
    for (int i = 0; i<imax; i++) free(cube[i]);

    free(cube);

}

void blk_tri(double *** lower, double *** main, double *** upper, int size_m, int num_m, double ** XB, double **X){

    /* Fill the test matrices. */

    for (int i = 0; i<3; i++){

        /* Fill the lower diagonal. */

        printf("-- ld: %lf\n", lower[0][0][i]);
        printf("-- ld: %lf\n", lower[0][1][i]);
        printf("-- ld: %lf\n", lower[1][0][i]);
        printf("-- ld: %lf\n", lower[1][1][i]);

        /* Fill the main matrix. */

        printf("-- lm: %lf\n", main[0][0][i]);
        printf("-- lm: %lf\n", main[0][1][i]);
        printf("-- lm: %lf\n", main[1][0][i]);
        printf("-- lm: %lf\n", main[1][1][i]);

        /* Fill the upper matrix. */

        printf("-- lu: %lf\n", upper[0][0][i]);
        printf("-- lu: %lf\n", upper[0][1][i]);
        printf("-- lu: %lf\n", upper[1][0][i]);
        printf("-- lu: %lf\n", upper[1][1][i]);

    }

    /* Fill the solution vector. */

    for (int i = 0; i<2; i++)
        for (int j = 0; j<2; j++)
            printf("-- xb: %lf\n",XB[i][j]);

}

int main(){

    double *** ld = alloc_dcube(2,2,3);
    double *** lm = alloc_dcube(2,2,3);
    double *** lu = alloc_dcube(2,2,3);
    double **  xb = alloc_dmatrix(2,3);
    double **  x  = alloc_dmatrix(2,3);

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

        lu[0][0][i] = 16.0;
        lu[0][1][i] = 9.0;
        lu[1][0][i] = 9.0;
        lu[1][1][i] = 6.0;

    }

    /* Fill the solution vector. */

    for (int i = 0; i<2; i++)
        for (int j = 0; j<2; j++)
            xb[i][j] = 3.0;

    /* Solve the system. */

    blk_tri(ld, lm, lu, 4, 3, xb, x);

    free_dcube(ld,2,2,3);
    free_dcube(lm,2,2,3);
    free_dcube(lu,2,2,3);
    free_dmatrix(xb,2,3);
    free_dmatrix(x,2,3);

    return 0;
}
