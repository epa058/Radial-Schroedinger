# Radial-Schrodinger

## Objectives

In this project, we aim to numerically solve the radial component of the Schroedinger equation
$$\left(-\frac{1}{2\mu} \frac{d^2}{dr^2} + \frac{\ell(\ell + 1)}{2\mu r^2} + V(r)\right) u_{n, \ell}(r) = E_{n, \ell} u_{n, \ell}(r)$$
for any given confining potential V(r). 

The primary focus is on determining the bound state wavefunctions $u_{n, \ell}(r)$ and the corresponding energy eigenvalues $E_{n, \ell}$. We employ a range of numerical analysis techniques, including the Numerov method, bisection search, and Newton's method, to achieve our computational objectives. 

The algorithm was tested on both the Coulomb potential and the Woods–Saxon potential; the results are presented in their respective folders.

## Description

We break down the core files and interdependencies: 

- **main_schroedinger.c**: Entry point of the program, coordinates the overall process of solving the Schroedinger equation. This file assembles various pieces of the algorithm and records the results. It mainly calls on:
  - **schroedinger.c**: Contains the central processes in solving the Schroedinger equation, with specific emphasis on wavefunction evolution. This file calls on:
    - **numerov.c**: Implements the Numerov method for numerical analysis.
    - **radial_eq_functions.c**: Defines the potentials and manages the forward and backward evolution of certain factors in the radial equation.

- **init.c**: Establishes the initial parameters and sets up the computation environment. Both **schroedinger.c** and **main_schroedinger.c** depend on it. This file calls on:
  - **extremum.c**: Contains routines for locating extremum points and their second derivative values.
  - **solve.c**: Contains various search techniques such as Newton's method and bisection search.
  - **derivatives.c**: Computes first and second derivatives.

Auxiliary files include:

- **vector_mtx.c**: Manages memory allocation for 1D and 2D arrays.
- **params.h** and **numerov_params.h**: Collections of static parameters from quantum mechanics and numerical analysis that assist in the calculations.

In summary, the program begins with **main_schroedinger.c**, which sets up the initial environment through **init.c**. It then leverages functions from files such as **schroedinger.c**, **numerov.c**, and **solve.c** to solve the radial Schroedinger equation numerically.
