#pragma once

#define BFN_S  0
#define BFN_P  1
#define BFN_D  2
#define BFN_SP -1

#define SPHERICAL 0
#define CARTESIAN 1

struct gto {  /* Gaussian-type orbital */
	
};

struct cgtf {  /* Contracted Gaussian-type function */
	int n;     /* principal quantum number */
	int L;     /* angular momentum */
	int nprim;
	double c[10];
	double exp[10];
};

struct basis_set {
	char name[30];
	int type;      /* spherical or cartesian */
	int size;      /* number of contracted GTO's */
	struct cgtf cgtfs[10];
};

/* atom-centered basis-function */
struct basis_function {
	int m;
	struct atom *a;
	struct cgtf *f;
};

void directive_basis();
void print_basis_summary();
