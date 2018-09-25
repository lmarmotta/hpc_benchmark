#include <stdlib.h>
#include <stdio.h>

/* gcc -O0 -ggdb -Wall -std=c99 -Ddgetrf=dgetrf_ -Ddgetri=dgetri_ dblkt_function.c -lblas -llapack */

/* To do...
 * [ ] - Copy solver function. */

extern void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);
extern void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);

void dmcpy(double ** mA, double ** mB, int size){
    for (int i = 0; i<size; i++)
        for (int j = 0; j<size; j++)
            mA[i][j] = mB[i][j];
}

void dmadd(double ** mA, double ** mB, double ** mC, int size){
    for (int i = 0; i<size; i++)
        for (int j = 0; j<size; j++)
            mC[i][j] = mA[i][j] + mB[i][j];
}

int ind2d(int i, int j, int tam){
    return (i)*(tam)+j; 
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

/* Matrix multiplication any size. Remember what I always forget:
 * [m x n] X [n x p] = [m x p] Where [a,b] are lines and columns. 
 *    A    x    B    =    C */ 

void dmgss(double ** mA, double ** mB, double ** mC, int m, int n, int p){

    for (int i=0; i<m; i++)
        for (int j=0; j<p; j++) mC[i][j] = 0.0;

    for (int i = 0; i < m; i++)
        for (int j = 0; j < p; j++)
            for (int k = 0; k < n; k++)
                mC[i][j] = mC[i][j] + mA[i][k]*mB[k][j];
}

/* Square matrix multiplication. */

void dmuls(double ** mA, double ** mB, double ** mC, int size){

    for (int i=0; i<size; i++)
        for (int j=0; j<size; j++) mC[i][j] = 0.0;

    for (int i=0; i<size; i++)
        for (int j=0; j<size; j++) 
            for (int k=0; k<size; k++)
                mC[i][j] = mC[i][j] + mA[i][k]*mB[k][j];
}

/* Alloc double vector. */

double * alloc_dvector(int imax){

    double * out = NULL;

    out = (double*)calloc(imax,sizeof(double));

    if (out == NULL) {printf("ERROR: double matrix allocation.\n"); exit(1);}

    return out;
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

        /* Allocate the block. */

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

void free_vector(double * vector, int imax){
    free(vector);
}

void blk_tri(double *** lower, double *** main, double *** upper, int size_m, int num_m, double ** XB, double ** X){

    /* Define auxiliar variables. */

    double *** gamm = alloc_dcube(size_m,size_m,num_m);

    double **  beta = alloc_dmatrix(size_m,num_m);

    double ** aux_copy = alloc_dmatrix(size_m,size_m);

    /* Get the first gamma. */

    for (int i = 0; i<size_m; i++)
        for (int j = 0; j<size_m; j++) aux_copy[i][j] = main[i][j][0];

    /* Invert matrix. */
            
    inv(aux_copy, size_m);

    double ** auxm1 = alloc_dmatrix(size_m,size_m);
    double ** auxm2 = alloc_dmatrix(size_m,size_m);

    for (int i = 0; i<size_m; i++)
        for (int j = 0; j<size_m; j++)
            auxm1[i][j] = upper[i][j][0];

    /* Multiply to obtain the first gamma. */

     dmuls(aux_copy, auxm1, auxm2, size_m);

    /* Get the gamma, aleluia !! */

    for (int i = 0; i<size_m; i++)
        for (int j = 0; j<size_m; j++)
            gamm[i][j][0] = auxm2[i][j];

    /* Get the rest of the gammas. */

    double ** aux_mult = alloc_dmatrix(size_m,size_m);
    double ** aux_summ = alloc_dmatrix(size_m,size_m);

    for (int m = 1; m<num_m-1; m++){

        for (int i = 0; i<size_m; i++){
            for (int j = 0; j<size_m; j++){ 
                auxm1[i][j] = lower[i][j][m];
                auxm2[i][j] = gamm[i][j][m-1];
            }
        }

        dmuls(auxm1, auxm2, aux_mult, size_m);

        for (int i = 0; i<size_m; i++){
            for (int j = 0; j<size_m; j++){
                aux_summ[i][j] = main[i][j][m] - aux_mult[i][j];
            }
        }

        inv(aux_summ, size_m);

        for (int i = 0; i<size_m; i++){
            for (int j = 0; j<size_m; j++){ 
                auxm1[i][j] = upper[i][j][m];
            }
        }

        dmuls(aux_summ, auxm1, auxm2, size_m);

        for (int i = 0; i<size_m; i++){
            for (int j = 0; j<size_m; j++){ 
                gamm[i][j][m] = auxm2[i][j];
            }
        }
    }

    /* Zero out our arrays. */

    for (int i = 0; i<size_m; i++)
        for (int j = 0; j<size_m; j++) aux_copy[i][j] = main[i][j][0];

    inv(aux_copy,size_m);

    double ** auxv1 = alloc_dmatrix(size_m,1);
    double ** auxv2 = alloc_dmatrix(size_m,1);

    for (int i = 0; i<size_m; i++) auxv1[i][0] = XB[i][0];

    dmgss(aux_copy, auxv1, auxv2, size_m, size_m, 1);

    double * aux_dumm = alloc_dvector(size_m);

    /* Form the first beta. */

    for (int i = 0; i<size_m; i++) beta[i][0] = auxv2[i][0]; 

    /* Go and grab the rest. */

    for (int m = 1; m<num_m; m++){

        for (int i = 0; i<size_m; i++)
            for (int j = 0; j<size_m; j++){
                aux_mult[i][j] = 0.0;
                aux_summ[i][j] = 0.0;
        }

        for (int i = 0; i<size_m; i++){
            for (int j = 0; j<size_m; j++){ 
                auxm1[i][j] = lower[i][j][m];
                auxm2[i][j] = gamm[i][j][m-1];
            }
        }

        dmuls(auxm1, auxm2, aux_mult, size_m);

        for (int i = 0; i<size_m; i++){
            for (int j = 0; j<size_m; j++){
                aux_summ[i][j] = main[i][j][m] - aux_mult[i][j];
            }
        }

        inv(aux_summ, size_m);

        for (int i = 0; i<size_m; i++){
            for (int j = 0; j<size_m; j++){
                auxm1[i][j] = lower[i][j][m];
            }
        }

        for (int i = 0; i<size_m; i++) auxv1[i][0] = beta[i][m-1];

        dmgss(auxm1, auxv1, auxv2, size_m, size_m, 1);

        for (int i = 0; i<size_m; i++)
            aux_dumm[i] = XB[i][m] - auxv2[i][0];

        for (int i = 0; i<size_m; i++) auxv1[i][0] = aux_dumm[i];

        dmgss(aux_summ, auxv1, auxv2, size_m, size_m, 1);

        for (int i = 0; i<size_m; i++) beta[i][m] = auxv2[i][0];

    }

    /* Start backward sweep. */

    for (int i = 0; i<size_m; i++)
        for (int j = 0; j<num_m; j++) X[i][j] = 0.0;

    for (int i = 0; i<size_m; i++) X[i][num_m-1] = beta[i][num_m-1];

    for (int m = num_m-2; m >= 0; m--){

        for (int i = 0; i<size_m; i++) 
            for (int j = 0; j<size_m; j++) auxm1[i][j] = gamm[i][j][m];

        for (int i = 0; i<size_m; i++) auxv1[i][0] = X[i][m+1];

        dmgss(auxm1, auxv1, auxv2, size_m, size_m, 1);

        for (int i = 0; i<size_m; i++) 
            X[i][m] = beta[i][m] - auxv2[i][0];
    }

    /* Avoid leaks baby. */

    free_dcube(gamm,size_m,size_m,num_m);

    free_dmatrix(beta,size_m,num_m);

    free_dmatrix(aux_copy,size_m,size_m);

    free_dmatrix(auxm1,size_m,size_m);

    free_dmatrix(auxm2,size_m,size_m);

    free_dmatrix(aux_mult,size_m,size_m);

    free_dmatrix(aux_summ,size_m,size_m);

    free_dmatrix(auxv1,size_m,1);

    free_dmatrix(auxv2,size_m,1);

    free_vector(aux_dumm,size_m);
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
        for (int j = 0; j<3; j++)
            xb[i][j] = 3.0;

    /* Solve the system. Here I am repeating the solution in order to catch any memory leak. */

    int repeat = 1;

    for (int r = 0; r<repeat; r++){
        printf("RUN: Entering sys solve\n");
        blk_tri(lm, ld, lu, 2, 3, xb, x);
        printf("SUCCESS: sys solve\n");
    }

    /* Print out hard coded and correct solution. */

    printf("x(1) = 5.8571429\n");
    printf("x(2) = -8.9220779\n");
    printf("x(3) = 0.5714286\n");
    printf("x(4) = -0.3116883\n");
    printf("x(5) = 1.8571429\n");
    printf("x(6) = -2.3766234\n");

    printf("\n");
    printf("Computed solution\n");
    printf("\n");

    for (int j = 0; j<3; j++){
        for (int i = 0; i<2; i++){
            printf("x(%d,%d) = %lf\n",i,j,x[i][j]);
        }
    }


    free_dcube(ld,2,2,3);
    free_dcube(lm,2,2,3);
    free_dcube(lu,2,2,3);
    free_dmatrix(xb,2,3);
    free_dmatrix(x,2,3);

    return 0;
}
