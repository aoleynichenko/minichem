#include <cblas.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "scf.h"
#include "ints.h"
#include "../util/util.h"
#include "../input/chem.h"
#include "../input/basis.h"

int double_eq(double a, double b)
{
	return fabs(a - b) < 1e-14;
}

int atoms_are_equal(struct atom *a, struct atom *b)
{
	return (a->Z == b->Z &&
			double_eq(a->r[0], b->r[0]) &&
			double_eq(a->r[1], b->r[1]) &&
			double_eq(a->r[2], b->r[2]));
}

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
	
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dim, dim, dim, alpha, P, dim, S, dim, beta, PS, dim);
	
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
	int info, lwork, lda = dim, alpha = 1.0, beta = 0.0;
	double *val, *X, *S12, *work, wkopt, *T;
	double *S12PS12;
	
	// first, compute eigenvalues of S
	val = (double *) qalloc(sizeof(double) * dim);
	X   = (double *) qalloc(bytes);
	memcpy(X, S, bytes);
	lwork = -1;
	dsyev_("V", "U", &dim, X, &lda, val, &wkopt, &lwork, &info);
	lwork = (int) wkopt;
    work = (double*) qalloc(lwork * sizeof(double));
    dsyev_("V", "U", &dim, X, &lda, val, work, &lwork, &info);
    
    if(info > 0)   // Check for convergence
		errquit("in Loewdin population analysis: LAPACK failed to orthogonalize overlap matrix");
	qfree(work, lwork * sizeof(double));
	
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
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dim, dim, dim, alpha, S12, dim, X, dim, beta, T, dim);
	cblas_dgemm(CblasRowMajor, CblasTrans,   CblasNoTrans, dim, dim, dim, alpha, X, dim, T, dim, beta, S12, dim);
	
	// now S12 contains S^1/2, we can calculate S12*P*S12
	// T <- P*S12
	// S12PS12 <- S12*T
	S12PS12 = X;
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dim, dim, dim, alpha, P, dim, S12, dim, beta, T, dim);
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dim, dim, dim, alpha, S12, dim, T, dim, beta, S12PS12, dim);
	
	printf("\n");
	printf("        Loewdin population analysis\n");
	printf("        ---------------------------\n\n");
	printf("      Atom     Population     Charge\n");
	general_pop_analysis(geom, bfns, S12PS12, dim);
	
	qfree(T,       bytes);
	qfree(S12,     bytes);
	qfree(S12PS12, bytes);
}

void multipole_moments(struct cart_mol *geom, struct basis_function *bfns, double *P, int dim)
{
	int i, j, k;
	double D;
	double d[] = {0, 0, 0};
	//double q[] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
	int e[] = {0, 0, 0};
	
	// the simplest implementation, we don't use P is Hermitian matrix and
	// symmetry of dipole moments integrals
	for (i = 0; i < dim; i++)
		for (j = 0; j < dim; j++) {
			for (k = 0; k < 3; k++) {
				e[0] = e[1] = e[2] = 0;
				e[k] = 1;
				d[k] -= P[i*dim+j] * aoint_multipole(&bfns[j], &bfns[i], e);
			}
		}
	for (i = 0; i < geom->size; i++)
		for (k = 0; k < 3; k++)
			d[k] += geom->atoms[i].Z * geom->atoms[i].r[k];
	D = sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);
	
	// TODO: quadrupole
	
	printf("               Multipole moments\n");
	printf("               -----------------\n");
	printf("  D       x            y            z\n");
	printf("  %13.8f%13.8f%13.8f\n", d[0], d[1], d[2]);
	printf("\n");
	printf("  |D| = %.8f a.u. = %.8f Debye\n", D, D*2.541746230211);
	printf("\n");
}

void effconf(struct cart_mol *geom, struct basis_function *bfns, double *P, double *S, int dim)
{
	int i, j;
	double tr = 0.0;
	int bytes = dim * dim * sizeof(double);
	double alpha = 1.0, beta  = 0.0;
	double *PS = NULL;
	
	if (!scf_options.effconf)
		return;
	
	PS = (double *) qalloc(bytes);
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, dim, dim, dim, alpha, P, dim, S, dim, beta, PS, dim);
	
	printf("\n");
	printf("        Effective atomic configurations\n");
	printf("        -------------------------------\n\n");
	/*printf("PS matrix:\n");
	for (i = 0; i < dim; i++) {
		for (j = 0; j < dim; j++) {
			printf("%8.2f", PS[i*dim+j]);
			if (i == j)
				tr += PS[i*dim+j];
		}
		printf("\n");
	}
	printf("\ntrace = %.2f\n", tr);*/
	
	j = 0;
	int count = 0;
	for (i = 0; i < geom->size; i++) {
		double pop, pop_s = 0.0, pop_p = 0.0;
		int Z = geom->atoms[i].Z;
		while (atoms_are_equal(bfns[j].a, &geom->atoms[i])) {
			pop = PS[j*dim+j];
			if (bfns[j].f->L == 0)
				pop_s += pop;
			else if (bfns[j].f->L == 1)
				pop_p += pop;
			if (++j == dim)
				break;
		}
		
		printf("        %-3d%-3s", count++, searchByZ(Z)->sym);
		if (Z > 2) {
			printf("  [He]");
			pop_s -= 2;
		}
		printf("  s{%.2f}", pop_s);
		if (pop_p > 1e-14)
			printf("  p{%.2f}", pop_p);
		printf("\n");
	}
	
	qfree(PS, bytes);
}









