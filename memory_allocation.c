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

void free_dmatrix(double ** matrix, int imax, int jmax){

    for (int i = 0; i<imax; i++)
        free(matrix[i]);

    free(matrix);
}

/* Alloc a double cube. */

double *** alloc_dcube(int imax, int jmax, int kmax){

    double *** out = NULL;

    out = (double***)calloc(imax,sizeof(double**));

    if (out == NULL) {printf("ERROR: double matrix allocation.\n"); exit(1);}

    for (int i = 0; i<imax; i++){
        out[i] = (double**)calloc(jmax,sizeof(double*));
        if (out[i] == NULL) {printf("ERROR: double matrix allocation.\n"); exit(1);}
    }

    for (int i = 0; i<imax; i++)
        for (int j = 0; j<jmax; j++)
            out[i][j] = (double*)calloc(kmax,sizeof(double));

    /* Initialize the matrix. */

    for (int i = 0; i<imax; i++)
        for (int j = 0; j<jmax; j++) 
            for (int k = 0; k<kmax; k++) out[i][j][k] = 0.0;

    return out;
}

void free_dcube(double *** cube, int imax, int jmax, int kmax){

    for (int i = 0; i<imax; i++)
        for (int j = 0; j<jmax; j++) free(cube[i][j]);
            
    for (int i = 0; i<imax; i++) free(cube[i]);

    free(cube);

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

    free_dmatrix(matrix,3,3);

    double *** cube = alloc_dcube(3,3,3);

    for(int i = 0; i<3; i++)
        for(int j = 0; j<3; j++)
            for(int k = 0; k<3; k++)
                cube[i][j][k] = 24.0;

    for(int i = 0; i<3; i++)
        for(int j = 0; j<3; j++)
            for(int k = 0; k<3; k++)
                printf("%lf\n",cube[i][j][k]);

    free_dcube(cube,3,3,3);
    

    return 0;
}
