// File extremum.h

#ifndef EXTREMUM_H // if EXTREMUM_H is not defined,
#define EXTREMUM_H // we define EXTREMUM_H

double Extremum_GetExtremum
(double (*func)(double), double x_init, double *curvature);

#endif // end if