/***********************************************************************
 * aoints.h
 * ========
 * 
 * Interface to minichem's molecular integrals evaluation module.
 * 
 * This module consists of:
 *  - aoints.h    this header file;
 *  - boys.c      provides implementation of the Boys function;
 *  - 1e.c        functions for normalization of primitive Gaussians,
 *                functions for evaluation of overlap, kinetic and
 *                potential energy integrals;
 *  - 2e.c        functions for evaluation electron-repulsion integrals
 *                (ERI's) and some helper functions for dealing with
 *                them.
 * Currently, only the Obara-Saika scheme is implemented.
 * 
 * For more details, see, for instance,
 *  T.Helgaker, P. Jorgensen, J. Olsen, "Molecular Electronic-Structure
 *  Theory".
 * 
 * 2016-2018 Alexander Oleynichenko
 **********************************************************************/

#ifndef AOINTS_H_INCLUDED
#define AOINTS_H_INCLUDED

#include "basis.h"
#include "chem.h"

void compute_aoints(Molecule_t *mol);

double boys(int n, double x);

double dist2(double *A, double *B);

double aoint_multipole(struct basis_function *fi,
					   struct basis_function *fj,
					   int *e);
double aoint_overlap(struct basis_function *fi,
					 struct basis_function *fj);
double aoint_kinetic(struct basis_function *fi,
					 struct basis_function *fj);
double aoint_potential(struct basis_function *fi,
					   struct basis_function *fj);
double aoint_eri(struct basis_function *fi,
				 struct basis_function *fj,
				 struct basis_function *fk,
				 struct basis_function *fl);
void print_ints(struct basis_function *funcs, int N);


#endif /* AOINTS_H_INCLUDED */







