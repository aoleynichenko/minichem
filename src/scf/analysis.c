/***********************************************************************
 * analysis.c
 * ==========
 * 
 * Very simple analysis of SCF function:
 *  - Mulliken and Loewdin population analysis
 *  - dx, dy, dz dipole moments
 * 
 * 2016-2018 Alexander Oleynichenko
 **********************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "aoints.h"
#include "basis.h"
#include "chem.h"
#include "linalg.h"
#include "scf.h"
#include "util.h"
#include "sys.h"


// Loewdin and Mulliken algorithms only contains different matrices,
// PS for Mulliken's and S12PS12 for Loewdin's. But the formula is the
// same for both cases.
void general_pop_analysis(struct cart_mol *geom, struct basis_function *bfns, double *matrix, int dim)
{
	int i, j, count = 1;
	double nelec_i, qA;
	
	j = 0;
	for (i = 0; i < geom->size; i++) {
		nelec_i = 0.0;
		while (atoms_are_equal(bfns[j].a, &geom->atoms[i])) {
			nelec_i += matrix[j*dim+j];
			if (++j == dim)
				break;
		}
		
		qA = geom->atoms[i].Z - nelec_i;
		printf("      %-3d%-3s%13.8f%13.8f\n", count++, searchByZ(geom->atoms[i].Z)->sym, nelec_i, qA);
	}
	printf("\n");
}

void mulliken(struct cart_mol *geom, struct basis_function *bfns, double *P, double *S, int dim)
{
	int bytes = dim * dim * sizeof(double);
	double alpha = 1.0, beta  = 0.0;
	double *PS = (double *) qalloc(bytes);
	
	linalg_square_dgemm(P, 'N', S, 'N', PS, dim);
	
	printf("\n");
	printf("        Mulliken population analysis\n");
	printf("        ----------------------------\n\n");
	printf("      Atom     Population     Charge\n");
	
	general_pop_analysis(geom, bfns, PS, dim);
	
	qfree(PS, bytes);
}

/* tr(PS) = tr(S^1/2*P*S^1/2)
 * X - transforms S to diagonal matrix
 * Implementation of sqrt(S) is rather complicated in C!
 */
void loewdin(struct cart_mol *geom, struct basis_function *bfns, double *P, double *S, int dim)
{
	int i, j;
	int bytes = dim * dim * sizeof(double);
	double *val, *X, *S12, *T;
	double *S12PS12;
	
	// first, compute eigenvalues of S
	val = (double *) qalloc(sizeof(double) * dim);  // eigenvalues
	X   = (double *) qalloc(bytes);                 // eigenvectors
	linalg_dsyev(S, val, X, dim);
	
	// now, val contains overlap matrix' eigenvalues
	// form s^1/2 [diagonal]
	S12 = (double *) qalloc(bytes);
	for (i = 0; i < dim; i++)
		for (j = 0; j < dim; j++)
			if (i == j)
				S12[i*dim+i] = sqrt(val[i]);
			else
				S12[i*dim+j] = 0.0;
	qfree(val, dim * sizeof(double));
	
	// undiagonalization of s12 = s^1/2
	// X is row-stored, so
	// S12 = X*s12*X'  --->  S12 = X'*s12*X
	// T is temporary matrix for intermediate result  -  s12*X
	T = (double *) qalloc(bytes);
	linalg_square_dgemm(S12, 'N', X, 'N', T,   dim);
	linalg_square_dgemm(X,   'T', T, 'N', S12, dim);
	
	// now S12 contains S^1/2, we can calculate S12*P*S12
	// T <- P*S12
	// S12PS12 <- S12*T
	S12PS12 = X;
	linalg_square_dgemm(P,   'N', S12, 'N', T,       dim);
	linalg_square_dgemm(S12, 'N', T,   'N', S12PS12, dim);
	
	printf("\n");
	printf("        Loewdin population analysis\n");
	printf("        ---------------------------\n\n");
	printf("      Atom     Population     Charge\n");

	general_pop_analysis(geom, bfns, S12PS12, dim);
	
	qfree(T,       bytes);
	qfree(S12,     bytes);
	qfree(S12PS12, bytes);
}


/***********************************************************************
 * scf_properties
 * 
 * properties at SCF level.
 * we don't use P is Hermitian matrix and permut symmetry of integrals.
 * TODO: XX, XY ... (quadrupole) moments
 **********************************************************************/
void scf_properties(struct cart_mol *geom, double *Pa, double *Pb, int dim)
{
	double *Ptot;
	int nbytes = dim * dim * sizeof(double);
	int i, j, k;
	int fd;
	double D;
	double d[] = {0, 0, 0};
	double *coord[3];  // array of three coord matrices (X, Y, Z)
	
	printf("\n");
	printf("\t\tSCF properties\n");
	printf("\t\t--------------\n");
	printf("\n");
	
	// alloc coord (X,Y,Z) matrices
	coord[0] = (double *) qalloc(nbytes);
	coord[1] = (double *) qalloc(nbytes);
	coord[2] = (double *) qalloc(nbytes);
	
	// read integrals from disk (as square matrices)
	fd = fastio_open("AOINTS1", "r");
	fastio_read_int(fd, &dim);
	// skip S, T, V
	fastio_read_doubles(fd, coord[0], dim*dim);
	fastio_read_doubles(fd, coord[0], dim*dim);
	fastio_read_doubles(fd, coord[0], dim*dim);	
	// read X, Y, Z
	fastio_read_doubles(fd, coord[0], dim*dim);
	fastio_read_doubles(fd, coord[1], dim*dim);
	fastio_read_doubles(fd, coord[2], dim*dim);
	fastio_close(fd);

	// construct "total" density matrix (alpha + beta)
	Ptot = (double *) qalloc(nbytes);
	for (i = 0; i < dim; i++) {
		for (j = 0; j < dim; j++) {
			Ptot[i*dim+j] = Pa[i*dim+j] + Pb[i*dim+j];
		}
	}
	
	/* DIPOLE MOMENT */
	// electronic contribution:
	for (i = 0; i < dim; i++)
		for (j = 0; j < dim; j++) {
			for (k = 0; k < 3; k++) {  // loop over X, Y, Z
				d[k] -= Ptot[i*dim+j] * coord[k][j*dim+i];
			}
		}
	// nuclear contribution:
	for (i = 0; i < geom->size; i++)
		for (k = 0; k < 3; k++)
			d[k] += geom->atoms[i].Z * geom->atoms[i].r[k];
	D = sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);
	// print
	printf("  dipole moment:\n");
	printf("  D       x            y            z\n");
	printf("  %13.8f%13.8f%13.8f\n", d[0], d[1], d[2]);
	printf("  |D| = %.8f a.u. = %.8f Debye\n", D, D*2.541746230211);
	printf("\n");

	// cleanup
	qfree(coord[0], nbytes);
	qfree(coord[1], nbytes);
	qfree(coord[2], nbytes);
	qfree(Ptot, nbytes);
}




