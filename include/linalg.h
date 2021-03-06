/***********************************************************************
 * linalg.h
 * ========
 * 
 * Wrappers for BLAS/LAPACK routines and some helper functions.
 * 
 * 2018 Alexander Oleynichenko
 **********************************************************************/

#ifndef LINALG_H_INCLUDED
#define LINALG_H_INCLUDED

#include <cblas.h>
#include <lapacke.h>

void print_matrix(char *annot, double *A, int dim);
void linalg_square_dgemm(double *A, char ta, double *B, char tb, double *C, int dim);
void linalg_dsyev(double *A, double *eigval, double *V, int dim);

#endif /* LINALG_H_INCLUDED */
