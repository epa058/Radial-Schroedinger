// File extremum.c
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "derivatives.h"
#include "extremum.h"
#include "solve.h"

typedef double (*FuncPT)(double); 
// FuncPT is a pointer to a function that takes in and outputs a double.
// We use this to pass an arbitrary function to another function.

FuncPT ORIG_FUNC; // A common variable, only valid within this file

double Extremum_DF(double x); // Used only within this file
// Extremum _GetExtremum finds the min or max near x_init
// This function returns the value x of the extremum
// The variable curvature has the value of the second derivative at the extremum

double Extremum_GetExtremum(FuncPT func, double x_init, double *curvature)
{
    double x, tol, ddf;
    int count;

    ORIG_FUNC = func; // To communicate with Extremum_DF

    tol = 1.0e-10;
    x = Solve_Newton(0.0, Extremum_DF, x_init, tol, &count);

    ddf = Derivative_SecondD(x, ORIG_FUNC); // ???

    *curvature = ddf; // Note that curvature is a pointer

    return x;
}
// End of Extremum_GetExtremum

// We use FuncPT func -> ORIG_FUNC because Solve_Newton only takes in double func(double) type functions

double Extremum_DF(double x)
{
    double f, df;

    df = Derivative_FirstD(x, ORIG_FUNC);

    return df;
}
// End of Extremum_DF

