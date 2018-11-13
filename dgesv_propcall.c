#include <stdlib.h>
#include <stdio.h>

extern void dgesv(int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b, int* ldb, int* info);

/* Allocate a double matrix. */

double ** alloc_dmatrix(int imax, int jmax){

    double ** out = NULL;

    out = (double**)calloc(imax,sizeof(double*));

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

/* Free double matrix. */

void free_dmatrix(double ** matrix, int imax){

    for (int i = 0; i<imax; i++)
        free(matrix[i]);

    free(matrix);
}

void axb_sol(double ** a_matrix, double * b_vector, double * x_vector, int sys_order){

    /* Auxiliary variables. */

    int n = sys_order;
    int nrhs = 1;
    int lda = n;
    int ldb = n;
    int info = 0;

    /* Alloc auxiliary array. */

    int * ipiv = (int*)calloc(n,sizeof(int));

    /* Linearize the matrix.*/

    double * lin_A = (double*)calloc(sys_order*sys_order,sizeof(double));

    for (int i = 0; i < sys_order; i++)
        for (int j = 0; j < sys_order; j++) 
            lin_A[(i)*(sys_order)+j] = a_matrix[j][i];

    /* Calling the system solver. */

    dgesv(&sys_order, &nrhs, lin_A, &lda, ipiv, b_vector, &ldb, &info);

    /* Check status. */

    if( info > 0 ) {
        printf( "The diagonal element of the triangular factor of A,\n" );
        printf( "U(%i,%i) is zero, so that A is singular;\n", info, info );
        printf( "the solution could not be computed.\n" );
        exit(1);
    }

    /* Pass the vector. */

    for (int i = 0; i<sys_order; i++) x_vector[i] = b_vector[i];

    /* Free the matrix. */

    free(ipiv);
    free(lin_A);

}

int main() {

    double ** A = alloc_dmatrix(3,3);
    double  * x = (double*)calloc(3,sizeof(double));
    double  * B = (double*)calloc(3,sizeof(double));

    /* Testing for memory leaks.*/

    for (int i = 0; i<10000; i++){

        A[0][0] = 1.0;
        A[0][1] = 1.0;
        A[0][2] = 1.0;

        A[1][0] = 5.0;
        A[1][1] = 3.0;
        A[1][2] = 2.0;

        A[2][0] =   0.0;
        A[2][1] =   1.0;
        A[2][2] = - 1.0;

        B[0] = 25.0;
        B[1] =  0.0;
        B[2] =  6.0;

        x[0] = 0.0;
        x[1] = 0.0;
        x[2] = 0.0;

        axb_sol(A, B, x, 3);

    }

    for (int i = 0; i<3; i++) printf("x[%d] = %lf\n",i,x[i]);

} 


