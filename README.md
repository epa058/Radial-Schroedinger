# Radial-Schrodinger

## Objectives

In this project, we aim to numerically solve the radial component of the Schrödinger equation
$$\left(-\frac{1}{2\mu} \frac{d^2}{dr^2} + \frac{\ell(\ell + 1)}{2\mu r^2} + V(r)\right) u_{n, \ell}(r) = E_{n, \ell} u_{n, \ell}(r)$$
for any given confining potential V(r). 

The primary focus is on determining the bound state wavefunctions $u_{n, \ell}(r)$ and the corresponding energy eigenvalues $E_{n, \ell}$. We employ a range of numerical analysis techniques, including the Numerov method, bisection search, and Newton's method, to achieve our computational objectives. 

The algorithm was tested on both the Coulomb potential and the Woods–Saxon potential; the results are presented in their respective folders.

## Description

We break down the main file functions and interdependencies:
**main_schroedinger.c**: This is the primary file users should open first. It acts as the entry point of the program, orchestrating the overall process of solving the Schrödinger equation. It assembles the pieces of the algorithm and records the results. Mainly, it calls on:
- **schroedinger.c**: This file contains the central processes in solving the Schrödinger equation, with specific emphasis on wavefunction evolution. This file calls on:
  - **numerov.c**: This file implements the Numerov method.
  - **radial_eq_functions.c**: This file contains functions defining the potentials as well as the forward and backward evolution functions for some ODEs.

**init.c**: This file sets up the necessary initial parameters and variables. It prepares the necessary environment and values for the solver to function accurately. It is called by **schroedinger.c** and **main_schroedinger.c**. This file calls on:
- **extremum.c**: This file contains functions to find extremum points and the value of the second derivative at those points.
- **solve.c**: This file contains various search methods like Newton's method or the bisection search.
- **derivatives.c**: This file is responsible for computing the first and second derivatives.

Auxiliary files include:
**vector_mtx.c**: This file is involved in allocating memory space for 1D and 2D arrays. It is called by **numerov.c** and **main_schroedinger.c**.
**params.h** and **numerov_params.h**, collections of unchanging parameters from quantum mechanics and numerical analysis that assist in the calculations.

In summary, main_schroedinger.c is the starting point, initializing the environment with init.c, and then leveraging the specific functionalities provided in the other files like numerov.c, schroedinger.c, and solve.c to carry out the numerical solution of the Schrödinger equation.
