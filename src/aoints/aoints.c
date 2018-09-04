/***********************************************************************
 * aoints.c
 * ========
 * 
 * Integral evaluation module -- driver routine.
 *
 * 2018 Alexander Oleynichenko
 **********************************************************************/

#include <math.h>
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
void compute_2e_aoints(BasisFunc_t *bfns, int nbas, char *filename);


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
	compute_2e_aoints(bfns, nbfns, "AOINTS2");
	
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
	int i, j, k, fd;
	double *A;   // temporary matrix
	int err;
	int xyz_pow[] = {0, 0, 0};
	
	timer_new_entry("1e", "One-electron integrals");
	timer_start("1e");
	
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

	printf("  X Y Z integrals\n");
	for (k = 0; k < 3; k++) {
		xyz_pow[0] = xyz_pow[1] = xyz_pow[2] = 0;
		xyz_pow[k] = 1;
		for (i = 0; i < nbas; i++)
			for (j = 0; j < nbas; j++) {
				struct basis_function *fi = &bfns[i];
				struct basis_function *fj = &bfns[j];
				A[nbas*i+j] = aoint_multipole(fi, fj, xyz_pow);
			}
		fastio_write_doubles(fd, A, nbas*nbas);
	}

	fastio_close(fd);
	free(A);
	
	timer_stop("1e");
	printf("done\n");
}


/***********************************************************************
 * compute_2e_aoints
 * 
 * Computes all permutationally-unique two-electron integrals and write
 * them to the binary file 'filename'.
 * Only non-zero integrals will are stored.
 * Format: <integral> <i1> <i2> <i3> <i4> (for each integral)
 **********************************************************************/
void compute_2e_aoints(BasisFunc_t *bfns, int nbas, char *filename)
{
	int m, n, p, q;
	int M = nbas;
	double V_mnpq;
	int i;
	int fd;
	typedef struct {
		double val;
		int i1, i2, i3, i4;
	} integral_t;
	const int BATCH_SIZE = 4096;
	integral_t buf[BATCH_SIZE];
	integral_t tmp;
	int n_uniq = 0, n_nonzero = 0;
	
	timer_new_entry("2e", "Two-electron integrals");
	timer_start("2e");
	
	printf("begin 2e integrals...\n");
	printf("  sizeof buf (bytes) = %d\n", sizeof(buf));
	printf("  sizeof integral (bytes) = %d\n", sizeof(integral_t));
	fd = fastio_open(filename, "w");
	
	i = 0;
	for (m = 0; m < M; m++)
	for (n = m; n < M; n++)
	for (p = m; p < M; p++)
	for (q = (p == m) ? n : p; q < M; q++) {
		n_uniq++;
		V_mnpq = aoint_eri(&bfns[m], &bfns[n], &bfns[p], &bfns[q]);
		// flush buffer to the disk
		if (i == BATCH_SIZE) {
			fastio_write(fd, buf, sizeof(buf));
			i = 0;
		}
		// store integral in the buffer if nonzero
		if (fabs(V_mnpq) < 1e-14) continue;

		n_nonzero++;
		tmp.val = V_mnpq;
		tmp.i1 = m;
		tmp.i2 = n;
		tmp.i3 = p;
		tmp.i4 = q;
		buf[i] = tmp;
		i++;
	}
	
	// store remaining integrals
	if (i > 0) {
		fastio_write(fd, buf, i*sizeof(integral_t));
	}
	
	printf("  # unique ERIs   = %d\n", n_uniq);
	printf("  # non-zero ERIs = %d\n", n_nonzero);
	
	fastio_close(fd);
	
	timer_stop("2e");
	printf("done\n");
}

