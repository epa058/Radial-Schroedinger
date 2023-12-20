// File vector_mtx.c
#include <stdio.h>
#include <stdlib.h>
#include "vector_mtx.h"

// Allocate memory space for a 1D array

double *vector_malloc(int nmax)
{
    double *pt;
    int n;

    pt = (double *)malloc(sizeof(double)*nmax);

    // Initialize all entries to zero
    for(n=0; n<nmax; n++)
    {
        pt[n] = 0.0;
    }

    return pt;
}
// End of vector_malloc

// Allocate memory space for a 2D array
double **mtx_malloc(int mmax, int nmax)
{
    double **pt;
    int m, n;

    pt = (double **)malloc(sizeof(double *)*mmax);

    for(m=0; m<mmax; m++)
    {
        pt[m] = (double *)malloc(sizeof(double)*nmax);
    }
    // End of m-loop

    // Initialize all entries to zero
    for(m=0; m<mmax; m++)
    {
        for(n=0; n<nmax; n++)
        {
            pt[m][n] = 0.0;
        }
    }
    // End of n-loop

    return pt;
}
// End of mtx_malloc

// Free memory space allocated byb mtx_malloc
void mtx_free(double **mtx, int mmax)
{
    int m;

    for(m=0; m<mmax; m++)
    {
        free(mtx[m]);
    }
    free(mtx);

    return;
}
// End of mtx_free