// File numerov.h
#ifndef NUMEROV_H // if NUMEROV_H is not defined,
#define NUMEROV_H // we define NUMEROV_H
#include "numerov_params.h"

void Numerov_Advance
(double *y, NumerovParams *Num_Params, DynamicVars *Dyn_Vars);

#endif // end if