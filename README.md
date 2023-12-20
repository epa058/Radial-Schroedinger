# Radial-Schroedinger

## Objectives

In this project, we aim to numerically solve the radial component of the Schroedinger equation
$$\left(-\frac{1}{2\mu} \frac{d^2}{dr^2} + \frac{\ell(\ell + 1)}{2\mu r^2} + V(r)\right) u_{n, \ell}(r) = E_{n, \ell} u_{n, \ell}(r)$$
for any confining potential $V(r)$. 

The primary focus was on determining the bound state wavefunctions $u_{n, \ell}(r)$ and their corresponding energy eigenvalues $E_{n, \ell}$. We employ a diverse range of numerical analysis techniques, including the Numerov method, bisection search, and Newton's method, to achieve this objective. 

The algorithm was tested on both the Coulomb potential and the Woods–Saxon potential; the results are presented in their respective folders.

## Description

We break down the core files and interdependencies: 

- **main_schroedinger.c**: This is the entry point of the program, coordinates the overall process. This file assembles various pieces of the algorithm and records the results. It mainly calls on:
  - **schroedinger.c**: This file contains the central processes in solving the Schroedinger equation, with a specific emphasis on wavefunction evolution. This file calls on:
    - **numerov.c**: This file implements the Numerov method from numerical analysis.
    - **radial_eq_functions.c**: This file defines the potentials and sets up the forward and backward evolution of certain factors in the radial equation.

- **init.c**: This file establishes the initial parameters and the computation environment. Both **schroedinger.c** and **main_schroedinger.c** depend on it. It calls on:
  - **extremum.c**: This file locates extremum points and evaluates their second derivatives.
  - **solve.c**: This file contains various search techniques such as Newton's method and bisection search.
  - **derivatives.c**: This file defines the first and second derivatives of functions.

Auxiliary files include:

- **vector_mtx.c**: This file manages memory allocation for 1D and 2D arrays.
- **params.h** and **numerov_params.h**: These headers contain collections of static parameters from quantum mechanics and numerical analysis that assist in the above calculations.

In summary, the program begins with **main_schroedinger.c**, which sets up the initial environment through **init.c**. It then leverages functions from files such as **schroedinger.c**, **numerov.c**, and **solve.c** to numerically solve the radial Schroedinger equation.
