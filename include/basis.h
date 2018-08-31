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
	int ijk[3];  // for gaussian function: x^i y^j z^k * exp(...)
	struct atom *a;
	struct cgtf *f;
	double norm[MAX_BASIS_SIZE];  // precomputed normalization factors
};
typedef struct basis_function BasisFunc_t;


struct shell {
	int size;
	struct basis_function *start;
};

/* explicit set of basis functions */
/* all basis functions are centered on atom at (x,y,z) */
extern struct basis_function *bfns;
extern struct shell *shells;

/* number of shells */
extern int nshells;

/* task size -- number of basis functions */
extern int nbfns;

void directive_basis();
void print_basis_summary();

struct cart_mol;   // chem.h for details
void form_atom_centered_bfns(struct cart_mol *molecule, struct basis_function **bfns, struct shell **shs, int *M, int *nsh);

#endif /* BASIS_H_INCLUDED */
