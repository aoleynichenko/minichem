# minichem
A very simple quantum chemistry program.

Minichem is a tiny educational project, created in order to better understand the computational methods of quantum chemistry. I will improve it while studying new methods.
For details, see A. Szabo, N. Ostlund, "Modern Quantim Chemistry", and T. Helgaker, P. Jorgensen, J. Olsen, "Molecular Electronic-Structure Theory". All example calculations were validated using NWChem.

Features:
 - Single-point energy calculation;
 - RHF and UHF methods;
 - DIIS convergence acceleration algorithm;
 - OpenMP parallelization (not yet implemented in C++ version);
 - Available elements: all periodic table (but only non-relativistic treatment is available!);
 - Available basis sets: cartesian and spherical with L up to g-functions;
 - Tiny basis set library (adopted from EMSL [1,2]).

Possibilities to be implemented in the future:
 - DIIS implementation for UHF;
 - More scalable OpenMP parallelization;
 - Analytic gradients and geometry optimization;
 - Projection population analysis;
 - Molecular symmetry;
 - I want to try to implement also four-index transformation and MP2

How to compile: <br>
 $ cd src <br>
 $ make <br>
 $ make install <br>
To compile minichem, you need the LAPACK library installed on your machine and any C++ compiler (I use g++ on my laptop).
For more details about installation of external libraries, see src/Makefile.

How to uninstall: <br>
 $ cd src <br>
 $ make uninstall <br>

For suggestions and questions, alexvoleynichenko@gmail.com, I look forward to any comments!

# Citations
[1] The Role of Databases in Support of Computational Chemistry Calculations
Feller, D., J. Comp. Chem., 17(13), 1571-1586, 1996.<br>
[2] Basis Set Exchange: A Community Database for Computational Sciences
Schuchardt, K.L., Didier, B.T., Elsethagen, T., Sun, L., Gurumoorthi, V., Chase, J., Li, J., and Windus, T.L.
J. Chem. Inf. Model., 47(3), 1045-1052, 2007, doi:10.1021/ci600510j. <br>
