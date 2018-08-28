/***********************************************************************
 * basis.h
 * =======
 * 
 * 2016-2018 Alexander Oleynichenko
 **********************************************************************/

#ifndef BASIS_H_INCLUDED
#define BASIS_H_INCLUDED

#define MAX_CNTRCT_LEN   30
#define MAX_BASIS_SIZE   30

#define BFN_S  0
#define BFN_P  1
#define BFN_D  2
#define BFN_SP -1

#define SPHERICAL 0    // not used
#define CARTESIAN 1


struct cgtf {  /* Contracted Gaussian-type function */
	int n;     /* principal quantum number */
	int L;     /* angular momentum */
	int nprim;
	double c[MAX_CNTRCT_LEN];
	double exp[MAX_CNTRCT_LEN];
};

struct basis_set {
	char name[30];
	int type;      /* spherical or cartesian */
	int size;      /* number of contracted GTO's */
	struct cgtf cgtfs[MAX_BASIS_SIZE];
};

/* atom-centered basis-function */
struct basis_function {
	int m;
	struct atom *a;
	struct cgtf *f;
};

void directive_basis();
void print_basis_summary();

#endif /* BASIS_H_INCLUDED */
