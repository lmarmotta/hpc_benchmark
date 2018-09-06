#include <stdlib.h>
#include <stdio.h>

/* Allocate a double vector. */

double * alloc_dvector(int imax){

    double ** out = NULL;

    out = (double**)calloc(imax,sizeof(double));

    if (out == NULL) {printf("ERROR: double matrix allocation.\n"); exit(1);}

    return 0;
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

/* Main */

int main(){

    double ** matrix = alloc_dmatrix(3,3);

    for(int i = 0; i<3; i++){
        for(int j = 0; j<3; j++){
            matrix[i][j] = -1.0;
            printf("%lf ",matrix[i][j]);
        }

        printf("\n");

    }

    return 0;
}
