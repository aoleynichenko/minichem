#pragma once

#include "../input/chem.h"
#include "../input/basis.h"

void scf_energy(struct cart_mol *);
void scalapack();

extern void   pdlawrite_( char **filenam, int *m, int *n, double *A, int *ia, int *ja, int *descA, int *irwrit, int *icwrit, double *work);
extern void   pdelset_( double *A, int *ia, int *ja, int *desca, double *alpha);
extern double pdlamch_( int *ictxt, char *cmach);
extern int    indxg2p_( int *indxglob, int *nb, int *iproc, int *isrcproc, int *nprocs);
extern int    indxg2l_( int *indxglob, int *nb, int *iproc, int *isrcproc, int *nprocs);
extern int    numroc_( int *n, int *nb, int *iproc, int *isrcproc, int *nprocs);
extern void   descinit_( int *desc, int *m, int *n, int *mb, int *nb, int *irsrc, int *icsrc,
                                int *ictxt, int *lld, int *info);
extern void   pdlaset_( char *uplo, int *m, int *n, double *alpha, double *beta, double *A, int *ia, int *ja, int *descA );
extern double pdlange_( char *norm, int *m, int *n, double *A, int *ia, int *ja, int *desca, double *work);
extern void   pdlacpy_( char *uplo, int *m, int *n, double *a, int *ia, int *ja, int *desca,
                                double *b, int *ib, int *jb, int *descb);
extern void   pdgesv_( int *n, int *nrhs, double *A, int *ia, int *ja, int *desca, int* ipiv,
                                double *B, int *ib, int *jb, int *descb, int *info);
extern void   pdgesvd_( char *jobu, char *jobvt, int *m, int *n, double *a, int *ia, int *ja, int *desca,
                                double *s, double *u, int *iu, int *ju, int *descu,
                                double *vt, int *ivt, int *jvt, int *descvt, double *work, int *lwork, int *info);
extern void   pdgemm_( char *TRANSA, char *TRANSB, int * M, int * N, int * K, double * ALPHA,
                                double * A, int * IA, int * JA, int * DESCA, double * B, int * IB, int * JB, int * DESCB,
                                double * BETA, double * C, int * IC, int * JC, int * DESCC );
extern int    indxg2p_( int *indxglob, int *nb, int *iproc, int *isrcproc, int *nprocs);
extern void pdsyev_( char *jobz, char *uplo, int *n,
                double *a, int *ia, int *ja, int *desca, double *w,
                double *z, int *iz, int *jz, int *descz,
                double *work, int *lwork, int *info );

#ifdef F77_WITH_NO_UNDERSCORE
#define   numroc_      numroc
#define   descinit_    descinit
#define   pdlamch_     pdlamch
#define   pdlange_     pdlange
#define   pdlacpy_     pdlacpy
#define   pdgesv_      pdgesv
#define   pdgemm_      pdgemm
#define   indxg2p_     indxg2p
#endif

extern void   Cblacs_pinfo( int* mypnum, int* nprocs);
extern void   Cblacs_get( int context, int request, int* value);
extern int    Cblacs_gridinit( int* context, char * order, int np_row, int np_col);
extern void   Cblacs_gridinfo( int context, int*  np_row, int* np_col, int*  my_row, int*  my_col);
extern void   Cblacs_gridexit( int context);
extern void   Cblacs_exit( int error_code);

extern void   dscal_( int *n, double *da, double *dx, int *incx);
extern void   dsyev_( char* jobz, char* uplo, int* n, double* a, int* lda,
                      double* w, double* work, int* lwork, int* info );
extern void dgesv_( int* n, int* nrhs, double* a, int* lda, int* ipiv,
                double* b, int* ldb, int* info );

struct scf_opt {
	int wavefuntype;
	int guess;
	int diis;
	// print options
	int print_1eov;
	int print_1eke;
	int print_1epe;
	int print_2eri;
	int print_final_vectors;
	// convergence options
	int maxiter;
	int diisbas;
	double conv_dens;
	double conv_en;
};

struct scf_timing_t {
	double time_diag;
	double time_fock;
	double time_ortho;
	double time_guess;
	double time_dens;
};

extern struct scf_timing_t scf_timing;
extern struct scf_opt scf_options;

typedef struct basis_function BasisFunc_t;
typedef struct cart_mol       Molecule_t;

void directive_scf();
void scf_init();
void uhf_loop(Molecule_t *molecule, BasisFunc_t *bfns, int M);
void rhf_loop(Molecule_t *molecule, BasisFunc_t *bfns, int M);
void compute_1e(double *Hcore, double *S, struct basis_function *bfns, int dim);
void diag_fock(double *F, double *X, double *C, double *en, int M);
void orthobasis(double *S, double *X, int dim);
double enuc(Molecule_t *geom);
double rmsdens(double *P1, double *P2, int M);

#define SCF_RHF  0
#define SCF_UHF  1
#define SCF_ROHF 2

#define GUESS_BARE 0
#define GUESS_EHT  1



void print_scf_options(struct scf_opt *opt);

// for initial guess
void guess_F_eht(double *F, double *S, struct basis_function *bfn, int n);

// properties
void mulliken(struct cart_mol *geom, struct basis_function *bfns, double *P, double *S, int dim);
void loewdin(struct cart_mol *geom, struct basis_function *bfns, double *P, double *S, int dim);
void multipole_moments(struct cart_mol *geom, struct basis_function *bfns, double *P, int dim);


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













