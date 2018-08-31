/* DIIS.
 * Идея: будем сохранять матрицы ошибок E[i] и Фока F[i] для каждого шага i
 * в односвязном списке. Новые матрицы будем добавлять в начало. Если
 * список достиг заданной длины (5-10 результатов итераций), то последний
 * его элемент освобождается. То есть, чем "старее" матрица, тем дальше она
 * от начала в этом списке.
 * 
 * Лучше всего эту часть алгоритма переписать на C++!
 */

#include <cblas.h>
#include <stdio.h>
#include <string.h>

#include "lapacke.h"
#include "scf.h"
#include "util.h"

DIISList_t *newDIISList(double *errm, double *fock, int dim)
{
	DIISList_t *leaf = (DIISList_t *) qalloc(sizeof(DIISList_t));
	leaf->next = NULL;
	leaf->dim = dim;
	leaf->E = errm;
	leaf->F = (double *) qalloc(dim*dim*sizeof(double));
	memcpy(leaf->F, fock, dim*dim*sizeof(double));
	return leaf;
}

/* Удаляет список данных DIIS начиная с элемента p включительно.
 * Работает рекурсивно, сначала уничтожает все элементы после данного,
 * потом уже данный.
 */
void removeDIISList(DIISList_t *p)
{
	if (!p)
		return;
	else {
		removeDIISList(p->next);
		qfree(p->E, sizeof(double) * p->dim * p->dim);
		qfree(p->F, sizeof(double) * p->dim * p->dim);
		qfree(p, sizeof(DIISList_t));
	}
}

DIISList_t *diis_store(DIISList_t *head, double *errm, double *fock, int dim, int diisbas)
{
	DIISList_t *leaf = (DIISList_t *) qalloc(sizeof(DIISList_t));
	DIISList_t *p = head, *last;
	int i = 0;
	
	*leaf = *head;  // copy all fields
	head->next = leaf;
	head->dim = dim;
	head->E = errm;
	head->F = (double *) qalloc(dim*dim*sizeof(double));
	memcpy(head->F, fock, dim*dim*sizeof(double));
	
	p = head;
	last = head;
	// обрезаем лишние элементы так, чтобы их было ровно diisbas в этом списке
	while (p && i < diisbas) {
		last = p;
		p = p->next;
		i++;
	}
	// теперь p указывает на первый из элементов, который должен быть уничтожен
	// идем по этому списку, уничтожая элементы и освобождая ресурсы
	// на деле, при равномерном росте цепи будет уничтожаться только
	// 1 элемент - последний, "хвост"
	last->next = NULL; // обрезаем список
	while (p) {
		last = p;
		p = p->next;
		qfree(last->E, sizeof(double) * last->dim * last->dim);
		qfree(last->F, sizeof(double) * last->dim * last->dim);
		qfree(last, sizeof(DIISList_t));
	}
	
	return leaf;
}

int diis_length(DIISList_t *p)
{
	int len = 0;
	while (p) {
		len++;
		p = p->next;
	}
	return len;
}

double matscalprod(double *A, double *B, int n)
{
	int i, j;
	double prod = 0.0;
	
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			prod += A[i*n+j]*B[i*n+j];
	return prod;
}

// ei = Fi*Di*S - S*Di*Fi
void make_error_matrix(double *e, double *F, double *D, double *S, int n)
{
	int bytes = sizeof(double) * n * n;
	double *T1; // temporary array
	double alpha = 1.0, // for dgemm
		   beta = 0.0;
	
	T1 = (double *) qalloc(bytes);
	
	// T1 <- D*F
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, alpha, D, n, F, n, beta, T1, n);
	// e <- S*T1, where 'e' is used as intermediate buffer to reduce memory usage
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, alpha, S, n, T1, n, beta, e, n);
	// T1 <- D*S
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, alpha, D, n, S, n, beta, T1, n);
	// e <- F*T1 - e
	beta = -1.0;
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, alpha, F, n, T1, n, beta, e, n);
	
	qfree(T1, bytes);
}


void diis_extrapolate(double *F, DIISList_t *diislist, int diisbas)
{
	int dim = diislist->dim, i, j, k;
	int diislen = diis_length(diislist);  // размерность итеративного подпространства
	int bdim = diislen + 1;   // B matrix size
	double *B = (double *) qalloc(bdim * bdim * sizeof(double));
	
	// compute Bij
	// максимально неэкономная схема
	// собираем сначала указатели на все матрицы ошибок
	double **errmats = (double **) qalloc(diislen * sizeof(double*));
	DIISList_t *p = diislist;
	i = 0;
	while (p) {
		errmats[i++] = p->E;
		p = p->next;
	}
	
	for (i = 0; i < diislen; i++) {
		for (j = 0; j < diislen; j++)
			B[i*bdim+j] = matscalprod(errmats[i], errmats[j], dim);
		B[i*bdim+j] = -1.0;
	}
	// -1 -1 -1 ... 0
	for (i = 0; i < bdim-1; i++)
		B[(bdim-1)*bdim+i] = -1.0;
	B[bdim*bdim-1] = 0.0;
	
	// A*x = B  --->  B*c = (0,0,0...-1) = r = right
	// B имеет размер bdim x bdim
	// r - это столбец высотой bdim
	// ld(B) = bdim
	// ld(r) = 1
	
	int info;
	int *ipiv = (int *) qalloc(sizeof(int)*bdim);
	double *right = (double *) qalloc(sizeof(double)*bdim);
	for (i = 0; i < bdim; i++)
		right[i] = 0.0;
	right[bdim-1] = -1;
	
        info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, bdim, 1, B, bdim, ipiv, right, 1);
        // Check for the exact singularity
        if( info > 0 ) {
                printf( "The diagonal element of the triangular factor of A,\n" );
                printf( "U(%i,%i) is zero, so that A is singular;\n", info, info );
                printf( "the solution could not be computed.\n" );
                errquit("error occured while DIIS extrapolation");
        }
	
	qfree(ipiv, sizeof(int)*bdim);
	qfree(B, bdim * bdim * sizeof(double));
	
	// и вот теперь создаем в F линейную комбинацию всех матриц Фока
	double **focks = errmats;  // to avoid reallocation
	p = diislist;
	i = 0;
	while (p) {
		focks[i++] = p->F;
		p = p->next;
	}
	for (i = 0; i < dim; i++)
		for (j = 0; j < dim; j++) {
			F[i*dim+j] = 0.0;
			for (k = 0; k < diislen; k++)
				F[i*dim+j] += focks[k][i*dim+j] * right[k];
		}
	
	qfree(focks, diislen * sizeof(double *));
	qfree(right, sizeof(double) * bdim);
}

