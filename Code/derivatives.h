// File derivatives.h

#ifndef DERIVATIVES_H // if DERIVATIVES_H is not defined,
#define DERIVATIVES_H // we define DERIVATIVES_H

double Derivative_FirstD
(double x, double (*func)(double));

double Derivative_SecondD
(double x, double (*func)(double));

#endif // end if