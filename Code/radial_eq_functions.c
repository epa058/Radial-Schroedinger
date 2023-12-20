// File radial_eq_functions.c
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "radial_eq_functions.h"
#include "numerov_params.h"
#include "params.h"

// Define the fine structure constant
#define alpha_EM (1.0 / 137.0)

// User defined function
// RadialEqFunctions_V(double r) is the only function one should change.

double RadialEqFunctions_V(double r) // 1/fm or 1/nm
{
    double f;
    double A, R_A, a, R0, V0;

    if(PARAM_DATA.nucA == 0.0)
    {
        // Use Coulomb potential
        f = -alpha_EM / r;
    } 
    else 
    {
        // Use Nuclear Woods-Saxon potential
        V0 = 50.0 / hbarc;  // MeV to 1/fm
        a = 0.7;
        R_A = PARAM_DATA.r0;
        f = -V0 / (1.0 + exp((r - R_A) / a));
    }
    
    return f;
}
// End of RadialEqFunctions_V

// No user serviceable parts from here on

// Veff(r) = V(r) + ell(ell+1)/(2 m r^2)
double RadialEqFunctions_Veff(double r)
{
    double f, ell, mass, V;

    ell = (double) PARAM_DATA.ell; // Convert ell from integer (in PARAM_DATA) to double
    mass = PARAM_DATA.mass;
    V = RadialEqFunctions_V(r); // Cleaner notation

    f = V + ell * (ell + 1) / (2 * mass * r * r);

    return f;
}
// End of RadialEqFunctions_Veff

// Forward evolution of F in y''=Fy
double RadialEqFunctions_F_Forward(double x, DynamicVars *Dyn_Vars)
{
    double x0, ka, r, f, g, Ea, Et, ell, eps, V;

    ell = (double) PARAM_DATA.ell;
    x0 = PARAM_DATA.x0;
    ka = PARAM_DATA.ka;
    Ea = PARAM_DATA.Ea;
    Et = Dyn_Vars->Et;

    // Small number to prevent x=0
    eps = 1.0e-15;
    x += eps;
    r = x / ka;
    V = RadialEqFunctions_V(r); // Cleaner notation

    f = ell * (ell + 1) / (x * x) + V / Ea + Et;

    return f;
}
// End of RadialEqFunctions_F_Forward

// Backward evolution of F in y''=Fy where y=x_f-x
double RadialEqFunctions_F_Backward(double y, DynamicVars *Dyn_Vars)
{
    double x0, ka, r, f, g, ell, Ea, Et, x, V;

    ell = (double) PARAM_DATA.ell;
    x0 = PARAM_DATA.x0;
    ka = PARAM_DATA.ka;
    Ea = PARAM_DATA.Ea;
    Et = Dyn_Vars->Et;

    x = Dyn_Vars->xf - y;
    r = x / ka;
    V = RadialEqFunctions_V(r); // Cleaner notation

    f = ell * (ell + 1) / (x * x) + V / Ea + Et;

    return f;
}
// End of RadialEqFunctions_F_Backward