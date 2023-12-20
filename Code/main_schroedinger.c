// File main_schroedinger.c
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "params.h"
#include "numerov_params.h"
#include "init.h"
#include "schroedinger.h"
#include "vector_mtx.h"

// Functions needed only in this file
// Read in Params parameters
void ReadIn_Params(char *input_file, DynamicVars *Dyn_Vars);

// Read in NumerovParams parameters
// This only reads in the number of data points for u_I and u_II
// The rest of the NumerovParams are read in by
// Schroedinger_InitializeCommonVar
void ReadIn_Num_Params(char *input_file_name, 
NumerovParams *Num_Params_f, NumerovParams *Num_Params_b);

// Record parameters
void Record_Params(NumerovParams Num_Params_f, NumerovParams Num_Params_b);

// Record results
void Record_Results(DynamicVars Dyn_Vars,
NumerovParams Num_Params_f, NumerovParams Num_Params_b,
double *yf, double *yb);

// This is declared as an extern in params.h
// As such, it needs to be declared outside params.h, but only once.
Params PARAM_DATA;

// This program is to be invoked as
// executable input_file1 input_file2
// the word "input_file1" and "input_file2" are then put into "argv" below
// argv[0] is the name of the executable (e.g. schroedinger)
// argv[1] is the name of the first input file (e.g. input_coulomb)
// argv[2] is the name of the second input file (e.g. input_n_params)

int main(int argc, char **argv)
{
    DynamicVars Dyn_Vars; // These parameters are calculated
    NumerovParams Num_Params_f; // For the forward evolution of u_I
    NumerovParams Num_Params_b; // For the backward evolution of u_II

    double *yf, *yb; // yf contains u_I, yb contains u_II

    // For PARAM_DATA
    // Reads in the initial data from the first input file
    ReadIn_Params(argv[1], &Dyn_Vars);
    
    // For Num_Params_f and Num_Params_b
    // Reads in the initial data from the second input file
    ReadIn_Num_Params(argv[2], &Num_Params_f, &Num_Params_b);

    // Get the energy and length scales
    Init_CalcScales();

    // Record the parameters
    Record_Params(Num_Params_f, Num_Params_b);

    // Allocate memory for the forward wavefunction yf
    yf = vector_malloc(Num_Params_f.nmax+1);

    // Allocate memory for the backward wavefunction yb
    yb = vector_malloc(Num_Params_b.nmax+1);

    Schroedinger_GetBoundState(&Dyn_Vars, &Num_Params_f, &Num_Params_b, yf, yb);
    Record_Results(Dyn_Vars, Num_Params_f, Num_Params_b, yf, yb);

    return 0;
}
// End of main

void ReadIn_Params(char *input_file, DynamicVars *Dyn_Vars)
{
FILE *input;
double x;
int ix;
char *mass_unit;

input = fopen(input_file, "r"); // Open the input file to "r"ead

// Read in the mass and put its value in x
fscanf(input, "%le", &x);
PARAM_DATA.mass = x / hbarc;

// Allocate enough memory to hold the mass unit
// 10 should be enough for any mass units
mass_unit = (char *) malloc(sizeof(char)*10);

/// Read in the mass unit
fscanf(input, "%s", mass_unit);
PARAM_DATA.mass_unit = mass_unit;

// Determine the length unit according to the mass unit
// strcmp is defined in string.h. It means "string compare."
// If the two strings are the same, then it returns 0.
// If the two strings are not the same, it returns non-zero.
if(strcmp(mass_unit, "eV")==0) 
{
    PARAM_DATA.length_unit = "nm";
}
else if(strcmp(mass_unit, "MeV")==0)
{
    PARAM_DATA.length_unit = "fm";
    }
else 
{
fprintf(stderr, "ReadIn_Params: %s is an unknown unit.\n", mass_unit);
fprintf(stderr, "Known units are eV and MeV.\n");
fprintf(stderr, "Exiting.\n");
exit(0);
}

// Read in the orbital angular momentum and put its value in ix
fscanf(input, "%d", &ix);
PARAM_DATA.ell = ix;

/* Note to self:
fscanf reads the file from the current position of the file pointer
until the end of the file or until the requested data is read.
Any subsequent calls to fscanf will start reading from 
where the previous call left off.
*/

// Read in the atomic mass A and put its value in x
fscanf(input, "%le", &x);
PARAM_DATA.nucA = x;

// Read in the atomic mass Z and put its value in x
fscanf(input, "%le", &x);
PARAM_DATA.nucZ = x;

// Read in the minimum value of Et for the eigenenergy search and put its value in x
fscanf(input, "%le", &x);
Dyn_Vars->Et_min = x;

// Read in the maximum value of Et for the eigenenergy search and put its value in x
fscanf(input, "%le", &x);
Dyn_Vars->Et_max = x;

fclose(input); // Always close an opened file
fprintf(stderr, "Done reading in.\n");

return;
}
// End of ReadIn_Params

void Record_Params(NumerovParams Num_Params_f, NumerovParams Num_Params_b)
{
    double x;
    int i;
    FILE *output;

    output = fopen("schroed_params.dat", "w"); // Open a file for writing

    // Note we multiply Ea by hbarc to recover units of energy
    fprintf(output, "mass = %e %s\n", PARAM_DATA.mass * hbarc, PARAM_DATA.mass_unit);
    // Note we multiply Ea by hbarc to recover units of energy
    fprintf(output, "Ea = %e %s\n", PARAM_DATA.Ea * hbarc, PARAM_DATA.mass_unit);
    // Note we multiply ka by hbarc to recover units of momentum
    fprintf(output, "ka = %e\n", PARAM_DATA.ka * hbarc);
    fprintf(output, "r0 = %e %s\n", PARAM_DATA.r0, PARAM_DATA.length_unit);
    fprintf(output, "x0 = %e\n", PARAM_DATA.x0);
    fprintf(output, "ell = %d\n", PARAM_DATA.ell);
    fprintf(output, "nucA = %e\n", PARAM_DATA.nucA);
    fprintf(output, "nucZ = %e\n", PARAM_DATA.nucZ);
    
    fprintf(output, "nmax_f = %d\n", Num_Params_f.nmax);
    fprintf(output, "nmax_b = %d\n", Num_Params_b.nmax);

    fclose(output); // Always close an open file

    return;
}
// End of Record_Params

// Read in NumerovParams data
void ReadIn_Num_Params(char *input_file_name, NumerovParams *Num_Params_f,
NumerovParams *Num_Params_b)
{
    FILE *input_file;
    double x;
    int ix;

    input_file = fopen(input_file_name, "r"); // Open the input file to "r"ead

    // Read in nmax for the forward evolution and put its value in x
    fscanf(input_file, "%d", &ix);
    Num_Params_f->nmax = ix;

    // Read in nmax for the backwards evolution and put its value in x
    fscanf(input_file, "%d", &ix);
    Num_Params_b->nmax = ix;

    fclose(input_file); // Always close an opened file

    return;
}
// End of ReadIn_Num_Params 

void Record_Results(DynamicVars Dyn_Vars, NumerovParams Num_Params_f, 
NumerovParams Num_Params_b, double *yf, double *yb)
{
    // Record the final results
    double Et, Eb;
    FILE *output;
    int n;
    double x;

    output = fopen("schroedinger.dat", "w");

    Et = Dyn_Vars.Et;
    fprintf(output, "Et = %e\n", Et);

    Eb = Et * PARAM_DATA.Ea;
    // Note we multiply Eb by hbarc to recover units of momentum
    fprintf(output, "Eb = %e\n", Eb * hbarc); 

    fclose(output);

    // Record the forward going solution
    output = fopen("yf.dat","w");
    for(n=0; n<=Num_Params_f.nmax; n++)
    {
        x = Num_Params_f.x_i + n * Num_Params_f.h; // Dimensionless
        x /= PARAM_DATA.ka; // Dimensions of 1/energy
        x *= hbarc; // Dimensions of length

        fprintf(output, "%e %e\n", x, yf[n] / yf[Num_Params_f.nmax]);
    }
    
    fclose(output);

    // Record the backward going solution
    output = fopen("yb.dat","w");
    for(n=0; n<=Num_Params_b.nmax; n++)
    {
        // 0 < x’ < x_f - x_c
        x = Dyn_Vars.xf - (Num_Params_b.x_i + n*Num_Params_b.h); // Dimensionless
        x /= PARAM_DATA.ka; // Dimensions of 1/energy
        x *= hbarc; // Dimensions of length

        fprintf(output, "%e %e\n", x, yb[n] / yb[Num_Params_b.nmax]);
    }
    fclose(output);
    
    return;
}
// End of Record_Results