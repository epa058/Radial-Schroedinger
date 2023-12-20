// File radial_eq_functions.h

#ifndef RADIAL_EQ_FUNCTIONS_H // if RADIAL_EQ_FUNCTIONS_H is not defined,
#define RADIAL_EQ_FUNCTIONS_H // we define RADIAL_EQ_FUNCTIONS_H
#include "params.h" // contains DynamicVars


// Potential energy in r
double RadialEqFunctions_V(double r);

// V_eff=V+ell(ell+1)/(2mr^2)
double RadialEqFunctions_Veff(double r);

// The F function in the right hand side of the ODE y’’ = Fy:
// Forward evolution from x=x_i
double RadialEqFunctions_F_Forward(double x, DynamicVars *Dyn_Vars);
// Backward evolution from x=x_f
double RadialEqFunctions_F_Backward(double x, DynamicVars *Dyn_Vars);

#endif // end if