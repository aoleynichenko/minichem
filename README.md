# minichem
A very simple quantum chemistry program.

Minichem is a tiny educational project, created in order to better understand the computational methods of quantum chemistry. I will improve it while studying new method.
For details, see A. Szabo, N. Ostlund, "Modern Quantim Chemistry", and T. Helgaker, P. Jorgensen, J. Olsen, "Molecular Electronic-Structure Theory".

Features:
 - Single-point energy calculation
 - Restricted Hartee-Fock method
 - DIIS convergence acceleration algorithm
 - OpenMP parallelization
 - Available elements: H-Ne
 - Available basis sets: any basis sets with only s and p functions
 - NWChem input file format
 
Possibilities to be implemented in the future:
 - More scalable OpenMP parallelization
 - Analytic gradients and geometry optimization
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

For help and questions, ao2310@yandex.ru
