/***********************************************************************
 * scf.h
 * =====
 * 
 * 2016-2018 Alexander Oleynichenko
 **********************************************************************/

#ifndef SCF_H_INCLUDED
#define SCF_H_INCLUDED

#include "chem.h"
#include "basis.h"

void scf_energy(struct cart_mol *);

void directive_scf();
void scf_init();
void uhf_loop(Molecule_t *molecule, BasisFunc_t *bfns, int M);
void rhf_loop(Molecule_t *molecule, BasisFunc_t *bfns, int M);
void read_1e_integrals(double *Hcore, double *S, struct basis_function *bfns, int dim);
void diag_fock(double *F, double *X, double *C, double *en, int M);
void orthobasis(double *S, double *X, int dim);
double enuc(Molecule_t *geom);
double rmsdens(double *P1, double *P2, int M);
void print_matrix(char *annot, double *A, int dim);
void scf_properties(struct cart_mol *geom, double *Pa, double *Pb, int dim);

#define SCF_RHF  0
#define SCF_UHF  1
#define SCF_ROHF 2

#define GUESS_BARE 0
#define GUESS_EHT  1


// for initial guess
void guess_F_eht(double *F, double *S, struct basis_function *bfn, int n);

// properties
void mulliken(struct cart_mol *geom, struct basis_function *bfns, double *P, double *S, int dim);
void loewdin(struct cart_mol *geom, struct basis_function *bfns, double *P, double *S, int dim);
void multipole_moments(struct cart_mol *geom, double *P, int dim);


// DIIS API
typedef struct diislist_t {
	struct diislist_t *next;
	double *E;
	double *F;
	double dim;
} DIISList_t;

void diis_extrapolate(double *F, DIISList_t *diislist, int diisbas);
void make_error_matrix(double *e, double *F, double *D, double *S, int n);
double matscalprod(double *A, double *B, int n);
int diis_length(DIISList_t *p);
DIISList_t *diis_store(DIISList_t *head, double *errm, double *fock, int dim, int diisbas);
DIISList_t *newDIISList(double *errm, double *fock, int dim);
void removeDIISList(DIISList_t *prev);
double maxerr(double *errmatrix, int n);

#endif /* SCF_H_INCLUDED */











