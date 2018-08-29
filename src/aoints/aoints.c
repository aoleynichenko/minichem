/***********************************************************************
 * aoints.c
 * ========
 * 
 * Integral evaluation module -- driver routine.
 *
 * 2018 Alexander Oleynichenko
 **********************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "aoints.h"
#include "basis.h"
#include "chem.h"
#include "input.h"
#include "util.h"
#include "sys.h"

// locally used routines
void compute_1e_aoints(BasisFunc_t *bfns, int nbas, char *filename);


void compute_aoints()
{
	printf("\n");
	printf("AO integrals evaluation module\n");
	printf("Output files:\n");
	printf("  AOINTS1   one-electron integrals\n");
	printf("  AOINTS2   two-electron integrals\n");
	printf("One-electron integrals evaluation algorithm: Obara-Saika\n");
	printf("Two-electron integrals evaluation algorithm: Obara-Saika\n");
	
	printf("generating atom-centered basis set...\n");
	form_atom_centered_bfns(&calc_info.molecule, &bfns, &shells, &nbfns, &nshells);
	printf("  # bfns   = %d\n", nbfns);
	printf("  # shells = %d\n", nshells);
	
	compute_1e_aoints(bfns, nbfns, "AOINTS1");
	
	printf("\n");
	line_separator();
	printf("\n");
}


/***********************************************************************
 * compute_1e_aoints
 * 
 * Computes all one-electron integrals and write them to the binary
 * file 'filename'.
 * All matrix elements are stored (NOT only one triangle!).
 * 
 * Format of the file:
 * <nbf>   number of basis functions
 * <S>     overlap matrix (nbf x nbf)
 * <T>     kinetic energy operator integrals
 * <NA>    potential energy -- nuclear attraction
 * <X>     dipole moment operator matrices (operators are X, Y, Z)
 * <Y>
 * <Z>
 **********************************************************************/
void compute_1e_aoints(BasisFunc_t *bfns, int nbas, char *filename)
{
	int i, j, fd;
	double *A;   // temporary matrix
	int err;
	
	// prepare
	A = (double *) malloc(nbas*nbas * sizeof(double));
	fd = fastio_open(filename, "w");
	
	printf("begin 1e integrals...\n");
	fastio_write_int(fd, nbas);
	
	printf("  overlap\n");
	for (i = 0; i < nbas; i++)
		for (j = 0; j < nbas; j++) {
			struct basis_function *fi = &bfns[i];
			struct basis_function *fj = &bfns[j];
			A[nbas*i+j] = aoint_overlap(fi, fj);
		}
	fastio_write_doubles(fd, A, nbas*nbas);
	
	printf("  kinetic energy\n");
	for (i = 0; i < nbas; i++)
		for (j = 0; j < nbas; j++) {
			struct basis_function *fi = &bfns[i];
			struct basis_function *fj = &bfns[j];
			A[nbas*i+j] = aoint_kinetic(fi, fj);
		}
	fastio_write_doubles(fd, A, nbas*nbas);
	
	printf("  nuclear attraction\n");
	for (i = 0; i < nbas; i++)
		for (j = 0; j < nbas; j++) {
			struct basis_function *fi = &bfns[i];
			struct basis_function *fj = &bfns[j];
			A[nbas*i+j] = aoint_potential(fi, fj);
		}
	fastio_write_doubles(fd, A, nbas*nbas);
	
	fastio_close(fd);
	free(A);
	printf("done\n");
}


