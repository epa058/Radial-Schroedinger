// File params.h

#ifndef PARAMS_H // if PARAMS_H is not defined,
#define PARAMS_H // we define PARAMS_H

#define hbarc (197.3) // MeV.fm = eV.nm

// Collection of parameters that do not change during the calculation.
// This will be set once in the beginning of the calculation and accessed via PARAM_DATA below.

typedef struct params
{
    double mass; // Particle mass
    double Ea; // Energy scale
    double ka; // Momentum scale
    double r0; // Length scale
    double x0; // k*r0
    int ell;

    char *mass_unit; // Unit of mass: either MeV or eV
    char *length_unit; // Unit of length: either fm or nm respectively

    double nucA; // Atomic mass
    double nucZ; // Atomic number

    // To bracket the energy eigenvalue search
    double Et_min;
    double Et_max;
}
Params;

// Collection of parameters that do change during the calculation and need to be passed frequently.
// Some of these parameters will be used in later typedef struct dynamic_vars.

typedef struct dynamic_vars
{
    double Eb; // Absolute value of the bound energy
    double kb; // \sqrt(2*mass*Eb)
    double rc; // Turning point radius
    double Et; // Eb/Ea
    double xc; // ka*rc
    double rf; // Last point
    double xf; // ka*rf

    // New additions to bracket the energy eigenvalue search
    double Et_min;
    double Et_max;
}
DynamicVars;

extern Params PARAM_DATA; // Run specific data

// "extern" specifies that this variable can be accessed 
// by any function as long as "params.h" is included.

// Any variable declared as an extern in 
// a header file must be declared once in a .c file

#endif // end if

