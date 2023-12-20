// File vector_mtx.h

#ifndef VECTOR_MTX_H // if VECTOR_MTX_H is not defined,
#define VECTOR_MTX_H // we define VECTOR_MTX_H

double *vector_malloc(int nmax); // allocates memory space for 1D arrays
double **mtx_malloc(int mmax, int nmax); // allocates memory space for 2D arrays
void mtx_free(double **mtx, int mmax); // frees the memory space allocated by mtx_malloc

#endif // end if