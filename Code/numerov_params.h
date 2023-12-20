// File numerov_params.h

#ifndef NUMEROV_PARAMS_H // if NUMEROV_PARAMS_H is not defined,
#define NUMEROV_PARAMS_H // we define NUMEROV_PARAMS_H
#include "params.h" // contains DynamicVars

typedef double (*Func_1D)(double, DynamicVars *);
// This defines Func_1D to be a pointer to functions that
// takes a double variable and a DynamicVars type variable
// and returns a double precision number.

typedef struct numerov_params
{
    double x_f; // maximum x
    double x_i; // minimum x
    double y_0; // y(x_i)
    double y_1; // y(x_i+h)
    int nmax; // number of sampling points
    double h; // step size: (x_f-x_i)/nmax

    // The function F in y''=Fy
    Func_1D NumerovFunc_F;
}
NumerovParams; // Define a data type named NumerovParams

#endif // end if