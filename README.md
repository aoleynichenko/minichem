# minichem
A very simple quantum chemistry program.

Minichem is a tiny educational project, created in order to better understand the computational methods in quantum chemistry. I will improve it while studying new method.

Features:
 - Single-point energy calculation
 - Restricted Hartee-Fock method
 - DIIS convergence acceleration algorithm
 - OpenMP parallelization
 - Available elements: H-Ne
 - Available basis sets: any basis sets with only s and p functions
 - NWChem file format
 
Possibilities to be implemented in the future:
 - More scalable OpenMP parallelization
 - Analytic gradients
 - Unrestricted Hartree-Fock method
 
How to compile:
$ cd src
$ make
$ make install
To compile project, you must have any MPI implementation (I use openmpi on my laptop) and gcc compiler.
In addition, you must have any implementation of LAPACK on your computer, see src/Makefile for details.

How to uninstall:
$ cd src
$ make uninstall
