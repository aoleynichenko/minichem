/***********************************************************************
 * linalg.c
 * ========
 * 
 * Wrappers for BLAS/LAPACK routines and some helper functions.
 * 
 * 2018 Alexander Oleynichenko
 **********************************************************************/

#include <stdio.h>
#include <string.h>

#include "linalg.h"
//#include "error.h"
#include "util.h"

/***********************************************************************
 * linalg_square_dgemm
 * 
 * Multiplication of square real general matrices.
 * Both matrices are stored in a row-major style.
 * 
 * Arguments:
 * A, B    input matrices
 * C       resulting matrix
 * ta, tb  transpose matrices A and B or not ('N' or 'T')
 * dim     matrix dimension
 **********************************************************************/
void linalg_square_dgemm(double *A, char ta, double *B, char tb, double *C, int dim)
{
	double alpha = 1.0, beta = 0.0;
	
	/*
	int dgemm_(char *transa, char *transb, int *m, int *n, int *k,
		double *alpha, double *a, int *lda, double *b, int *ldb,
		double *beta, double *c, int *ldc);
	*/
	
	// C = 1.0 * A * B + 0.0 * C 
	cblas_dgemm(CblasRowMajor,
	            ta == 'N' ? CblasNoTrans : CblasTrans,
	            tb == 'N' ? CblasNoTrans : CblasTrans,
	            dim, dim, dim, alpha, A, dim, B, dim, beta, C, dim);
}


/***********************************************************************
 * linalg_dsyev
 * 
 * Solver for the eigenvalue problem.
 * Matrix is stored in a row-major style.
 * Eigenvectors will be stored in V row-wise.
 * A is not changed.
 * len(eigval) = dim
 **********************************************************************/
void linalg_dsyev(double *A, double *eigval, double *V, int dim)
{

	int info, lwork, lda = dim;
	double *work, wkopt;
	
	// save A
	memcpy(V, A, dim*dim*sizeof(double));
	
	lwork = -1;
	dsyev_("V", "U", &dim, V, &lda, eigval, &wkopt, &lwork, &info);
	lwork = (int) wkopt;
    work = (double*) qalloc(lwork * sizeof(double));
    dsyev_("V", "U", &dim, V, &lda, eigval, work, &lwork, &info);
    
    if(info > 0)   // Check for convergence
		errquit("LAPACK failed to orthogonalize overlap matrix");
	qfree(work, lwork * sizeof(double));
}


/***********************************************************************
 * print_matrix
 * 
 * Prints square matrix.
 **********************************************************************/
void print_matrix(char *annot, double *A, int dim)
{
	int i, j;
	
	if (annot)
		printf("%s\n", annot);
	for (i = 0; i < dim; i++) {
		for (j = 0; j < dim; j++)
			printf("%13.8f", A[i*dim+j]);
		printf("\n");
	}
}
