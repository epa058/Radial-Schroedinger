// File init.c
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "params.h"
#include "init.h"
#include "radial_eq_functions.h"
#include "solve.h"
#include "extremum.h"
#include "derivatives.h"

// The function to minimize to get Ea and ka
// Equals V(r) + 1/(2 m r^2)
double Init_Rmin_Function(double r);

void Init_CalcScales(void)
{
    double r_min, r_init, E_min, Eb, mass, r_c, curvature;
    double ka, x0;
    double R_A;
    FILE *output;
    
    mass = PARAM_DATA.mass;

    if(PARAM_DATA.nucA == 0.0) // Not a nuclear potential
    {
        r_init = 0.01; // Start the iteration
        r_min = Extremum_GetExtremum(Init_Rmin_Function, r_init, &curvature);

        PARAM_DATA.r0 = r_min;
        E_min = Init_Rmin_Function(r_min);
        PARAM_DATA.Ea = fabs(E_min);
    }
    else // Nuclear potential. We know scales in this case.
    {
        PARAM_DATA.r0 = 1.3 * pow(PARAM_DATA.nucA, 1.0 / 3.0); // Nuclear radius
        PARAM_DATA.Ea = 50.0 / hbarc; // About 50 MeV
    }

    ka = sqrt(2.0*mass*PARAM_DATA.Ea);
    PARAM_DATA.ka = ka;

    x0 = ka*PARAM_DATA.r0;
    PARAM_DATA.x0 = x0;

    return;
}
// End of Init_CalcScales

double Init_Rmin_Function(double r) // Function to minimize = V(r) + 1/(2 m r^2)
{
    double f, mass;
    double V = RadialEqFunctions_V(r); // Cleaner notation

    mass = PARAM_DATA.mass;
    
    f = V + 1.0 / (2.0 * mass * r * r);

    return f;
}
// End of Init_Rmin_Function