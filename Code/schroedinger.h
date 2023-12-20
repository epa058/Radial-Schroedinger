// File schroedinger.h

#ifndef SCHROEDINGER_H // if SCHROEDINGER_H is not defined,
#define SCHROEDINGER_H // we define SCHROEDINGER_H
#include "params.h"
#include "numerov_params.h"

void Schroedinger_GetBoundState
(DynamicVars *Dyn_Vars, NumerovParams *Num_Params_f, 
NumerovParams *Num_Params_b, double *yf, double *yb);

#endif // end if