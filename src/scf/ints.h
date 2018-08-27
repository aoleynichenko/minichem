/***********************************************************************
 * ints.h
 * 
 * Interface to minichem's molecular integrals evaluation module.
 * 
 * This module consists of:
 *  - ints.h      this header file;
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
 * (c) Alexander Oleynichenko, 2016
 ***********************************************************************/

#pragma once

#include "../input/basis.h"
#include "../input/chem.h"

double boys(int n, double x);

double cS00(double *A, double *B, double a, double b);
double sS00(double *A, double *B, double a, double b);
double cS10(double *A, double *B, double a, double b, int pm);
double cS11(double *A, double *B, double a, double b, int m1, int m2);
double cT00(double *A, double *B, double a, double b);
double cT01(double *A, double *B, double a, double b, int pm);
double cT11(double *A, double *B, double a, double b, int m1, int m2);

double N00(double a, double b);
double N01(double a, double b);
double N10(double a, double b);
double N11(double a, double b);

double dist2(double *A, double *B);

double N0(double a);
double N1(double a);

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










