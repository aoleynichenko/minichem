#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>

#include "input.h"
#include "basis.h"
#include "util.h"
#include "lapacke.h"

#include <omp.h>
#include <mpi.h>
#include <cblas.h>

#include "scf.h"
#include "ints.h"
#include "molden.h"

struct shell {
	int size;
	struct basis_function *start;
};

void print_blacs_grid_info();
void config_grid(int nprocs, int *nprow, int *npcol);
void form_atom_centered_bfns(struct cart_mol *molecule, struct basis_function **bfns, struct shell **shs, int *M, int *nshells);
void print_matrix(char *annot, double *A, int dim);

static int mpi_rank = -1;
static int mpi_size = -1;

/* scf parameters */
struct scf_opt scf_options;

/* task size - number of basis functions */
int M;

/* number of shells */
int nshells;

/* general BLACS variables */
int ictxt;   /* BLACS context */
int iam;     /* my number in BLACS grid */
int nprocs;  /* number of processes, invoked in BLACS grid */

/* process grid parameters */
int nprow;   /* number of rows */
int npcol;   /* number of columns */
int myrow;   /* my position in process grid - "y" */
int mycol;   /* my position in process grid - "x" */

/* size of submatrices */
int vb = 64; /* preferred vertical size of the block */
int hb = 64; /* horizontal size */
int vsize;   /* actual vertical size of submatrices */
int hsize;   /* actual horizontal size */

/* parameters for MPI scattering/gathering operations */
int *counts = 0;
int *displs = 0;

/* molecule */
struct cart_mol *geom;
int Nelecs;
int Nalpha;
int Nbeta;

/* matrices */
double *S;      /* overlap matrix */
double *X;		/* transformation matrix */
double *P;      /* density matrix */
double *Hcore;  /* core Hamiltonian */
double *F;      /* Fock matrix */
double *C;      /* vectors */
double *E;		/* energies of orbitals, F's eigenvalues */
double *OLDP;

/* explicit set of basis functions */
/* all basis functions are centered on atom at (x,y,z) */
struct basis_function *bfns;
struct shell *shells;

struct scf_timing_t scf_timing;


void scf_energy(struct cart_mol *molecule)
{
	int izero = 0;
	int displ, sqsize, startx, starty;
	int i, j;

	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

	if (mpi_size != 1)
		errquit("Sorry! Parallel SCF module hasn't implemented yet. Please, run minichem on one node!");

	printf("          ********************************\n");
	printf("          *      Parallel SCF Module     *\n");
	printf("          ********************************\n\n");
	
	scf_timing.time_dens  = 0.0;
	scf_timing.time_diag  = 0.0;
	scf_timing.time_fock  = 0.0;
	scf_timing.time_guess = 0.0;
	scf_timing.time_ortho = 0.0;
	omp_set_num_threads(calc_info.nproc);
		
	Nelecs = nalphabeta(molecule, &Nalpha, &Nbeta);
	// automatically set wavefunction type
	if (Nalpha == Nbeta)
		scf_options.wavefuntype = SCF_RHF;
	else
		scf_options.wavefuntype = SCF_UHF;
	
	print_scf_options(&scf_options);
	mol_summary(&calc_info.molecule);
	geom = molecule;
	form_atom_centered_bfns(molecule, &bfns, &shells, &M, &nshells);  // create atom-centered basis set
	
	if (scf_options.wavefuntype == SCF_RHF)
		rhf_loop(geom, bfns, M);
	else
		uhf_loop(geom, bfns, M);
}


void print_scf_options(struct scf_opt *opt)
{
	printf("                    SCF Parameters\n");
	printf("                    --------------\n");
	printf("               Wavefunction : %s\n", opt->wavefuntype == SCF_RHF ? "RHF" : "UHF");
	printf("             Max iterations : %d\n", opt->maxiter);
	printf("               Density conv : %g\n", opt->conv_dens);
	printf("                Energy conv : %g\n", opt->conv_en);
	printf("                       DIIS : %s\n", opt->diis ? "on" : "off");
	if (opt->diis)
	printf("                    diisbas : %d\n", opt->diisbas);
	printf("\n");
}

// Ортогонализация базиса
// Непараллельная реализация с использованием LAPACK
void orthobasis(double *S, double *X, int dim)
{
	int i, j;
	int mbytes = dim * dim * sizeof(double); // N bytes in dim x dim matrix
	int vbytes = dim * sizeof(double);       // N bytes in vector
	int info, lwork, lda = dim;
	double thresh = 1e-8;
	double t0 = MPI_Wtime();
	double *val, *work, wkopt;
	
	printf("Basis orthogonalization algorithm: canonical\n");
	printf("Basis functions elimination threshold: %g\n", thresh);
	
	// Solve eigenproblem for overlap matrix S
	val = (double *) qalloc(vbytes);
	memcpy(X, S, mbytes);
	lwork = -1;
	dsyev_("V", "U", &dim, X, &lda, val, &wkopt, &lwork, &info);
	lwork = (int) wkopt;
    work = (double*) qalloc(lwork * sizeof(double));
    dsyev_("V", "U", &dim, X, &lda, val, work, &lwork, &info);
    
    if(info > 0)   // Check for convergence
		errquit("LAPACK failed to orthogonalize overlap matrix");
	qfree(work, lwork * sizeof(double));
	
	printf("Overlap matrix lowest eigenvalue = %.8f\n", val[0]);
	
	if (val[0] < thresh)
		errquit("negative defined overlap matrix or zero eigenvalues!");
	
	for (i = 0; i < dim; i++)
		for (j = 0; j < dim; j++) {
			X[i*dim+j] = X[i*dim+j] / sqrt(val[i]);
		}
	
	qfree(val, vbytes);
	
	scf_timing.time_ortho = MPI_Wtime() - t0;
	printf("AO basis orthogonalization done in %.6f sec\n", MPI_Wtime()-t0);
}

void scf_init()
{
	scf_options.wavefuntype = SCF_RHF;
	scf_options.guess = GUESS_EHT;
	scf_options.diis = 1;
	// print options
	scf_options.print_1eov = 0;
	scf_options.print_1eke = 0;
	scf_options.print_1epe = 0;
	scf_options.print_2eri = 0;
	// convergence options
	scf_options.maxiter = 50;
	scf_options.diisbas = 5;
	scf_options.conv_dens = 1e-6;
	scf_options.conv_en = 1e-7;
}


void compute_1e(double *Hcore, double *S, struct basis_function *bfns, int dim)
{
	int i, j;
	double t1;
	
	if (mpi_rank == 0) {
		printf("One-electron integrals evaluation algorithm: Obara-Saika\n");
		printf("Two-electron integrals evaluation algorithm: Obara-Saika\n");
	}
	
	t1 = MPI_Wtime();
	
	// Overlap, Kinetic & Potential matrices
	for (i = 0; i < dim; i++)
		for (j = 0; j < dim; j++) {
			struct basis_function *fi = &bfns[i];
			struct basis_function *fj = &bfns[j];
			S[dim*i+j] = aoint_overlap(fi, fj);
			Hcore[dim*i+j] = aoint_kinetic(fi, fj) + aoint_potential(fi, fj);
		}
	
	MPI_Barrier(MPI_COMM_WORLD);
	if (mpi_rank == 0)
		printf("One-electron integrals done in %.6f sec\n", (MPI_Wtime()-t1)/1000);
	
	print_ints(bfns, dim);
}


double distance(double *A, double *B)
{
	return sqrt((A[0] - B[0])*(A[0] - B[0]) +
				(A[1] - B[1])*(A[1] - B[1]) +
				(A[2] - B[2])*(A[2] - B[2]));
}

double enuc(Molecule_t *geom)
{
	double e = 0.0;
	int i, j;
	
	for (i = 0; i < geom->size; i++)
		for (j = i+1; j < geom->size; j++)
			e += geom->atoms[i].Z*geom->atoms[j].Z/distance(geom->atoms[i].r, geom->atoms[j].r);
	
	return e;
}

double rmsdens(double *P1, double *P2, int M)
{
	int i, j;
	double s = 0.0, d;
	
	#pragma omp parallel for private(j,d) reduction(+:s)
	for (i = 0; i < M; i++)
		for (j = 0; j < M; j++) {
			d = P1[i*M+j]-P2[i*M+j];
			s += d*d;
		}
	return sqrt(s/4.0);
}

double maxerr(double *errmatrix, int n)
{
	int i, j;
	double max = DBL_MIN;
	
	for (i = 0; i < n; i++)
		for (j = i; j < n; j++)
			if (errmatrix[i*n+j] > max)
				max = errmatrix[i*n+j];
	return max;
}

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

#include "cblas.h"
int dgemm_(char *transa, char *transb, int *m, int *n, int *k,
	double *alpha, double *a, int *lda, double *b, int *ldb,
	double *beta, double *c, int *ldc);

void diag_fock(double *F, double *X, double *C, double *en, int M)
{
	int i, j;
	int info, lwork, lda = M, ldb = M, ldc = M;
	double alpha = 1.0;
	double beta = 0.0;
	double *work, wkopt, t0 = MPI_Wtime();
	
	double *Temp = (double *) qalloc(M*M*sizeof(double));
	double *TempF = (double *) qalloc(M*M*sizeof(double));
	memcpy(TempF, F, M*M*sizeof(double));
	
	// F = X*F*X', X stored in transposed form
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,   M, M, M, alpha, TempF, M, X, M, beta, Temp, M);
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, M, M, alpha, X, M, Temp, M, beta, TempF, M);
	
	memcpy(C, TempF, sizeof(double) * M * M);
	lwork = -1;
	dsyev_("V", "U", &M, C, &lda, en, &wkopt, &lwork, &info);
	lwork = (int) wkopt;
    work = (double*) qalloc(lwork * sizeof(double));
    dsyev_("V", "U", &M, C, &lda, en, work, &lwork, &info);
    
    if(info > 0)   // Check for convergence
		errquit("LAPACK failed to orthogonalize Fock matrix");
	qfree(work, lwork * sizeof(double));
	
	// transform C:  C = X*C
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, M, M, alpha, C, M, X, M, beta, Temp, M);
	memcpy(C, Temp, M*M*sizeof(double));
	// Now, transformed vectors are 'lying' in the C matrix rows
	
	qfree(Temp, M*M*sizeof(double));
	qfree(TempF, M*M*sizeof(double));
	scf_timing.time_diag += MPI_Wtime() - t0;
}

/* Creates atom-centered basis functions from molecular data.
 * Returns: vector of atom-centered functions bfns with length M.
 * This function should be executed by all processes.
 * */
void form_atom_centered_bfns(struct cart_mol *molecule, struct basis_function **bfns, struct shell **shs, int *M, int *nsh)
{
	int i, j;
	int K = 0;
	int shn = 0;
	struct basis_function *p;
	struct shell *s;

	for (i = 0; i < molecule->size; i++) {
		struct elem_info *e = searchByZ(molecule->atoms[i].Z);
		if (e) {
			for (j = 0; j < e->bas->size; j++) {
				K += 2*e->bas->cgtfs[j].L + 1;  // пока и так сойдет
				shn++;
			}
		}
	}
	p = (struct basis_function *) malloc(K * sizeof(struct basis_function));
	s = (struct shell *) malloc(shn * sizeof(struct shell));
	*bfns = p;
	*shs = s;
	*M = K;
	*nsh = shn;
	for (i = 0; i < molecule->size; i++) {
		struct elem_info *e = searchByZ(molecule->atoms[i].Z);
		if (!e)
			continue;
		for (j = 0; j < e->bas->size; j++) {  // add atom-centered bfn
			int L = e->bas->cgtfs[j].L;
			int m;
			
			(*s).size = 2*L+1;
			(*s).start = p;
			for (m = 0; m < 2*L+1; m++) {
				(*p).m = m;
				(*p).f = &e->bas->cgtfs[j];
				(*p).a = &molecule->atoms[i];
				p++;
			}
			s++;
		}
	}
}


/* This function prints information about BLACS process grid and
 * submatrices to stdout. This function is important for debugging. */
void print_blacs_grid_info()
{
	// от каждого потока нам надо:
	//  x - положение в сетке процессов по горизонтали
	//  y - положение в сетке процессов по вертикали
	//  h - вертикальная   размерность подматрицы
	//  v - горизонтальная размерность подматрицы

	if (mpi_rank == 0) {
		int i, j;
		int *info = (int *) malloc(4 * mpi_size * sizeof(int));
		int *map =  (int *) malloc(mpi_size * sizeof(int));
		MPI_Request *reqs = (MPI_Request *) malloc((mpi_size-1)*sizeof(MPI_Request));
		MPI_Status *stats = (MPI_Status *)  malloc((mpi_size-1)*sizeof(MPI_Status));

		info[0] = myrow;
		info[1] = mycol;
		info[2] = vsize;
		info[3] = hsize;

		for (i = 0; i < mpi_size-1; i++)  // receive data from all other processes
			MPI_Irecv(&info[4*(i+1)], 4, MPI_INT, i+1, 0, MPI_COMM_WORLD, &reqs[i]);
		MPI_Waitall(mpi_size-1, reqs, stats);

		printf("\nBLACS grid information:\n");
		printf("  nprocs = %d\n", nprocs);
		printf("  nprow x npcol: %d x %d\n", nprow, npcol);
		printf("  preferred submatrix size (vert x hor): %d x %d\n", vb, hb);
		//for (i = 0; i < mpi_size; i++)
			//printf("%d: (%d,%d), %dx%d\n", i, info[4*i], info[4*i+1], info[4*i+2], info[4*i+3]);
		printf("\n                    PROCESS GRID\n");

		for (i = 0; i < mpi_size; i++) { /* for more beautiful table */
			int y = info[4*i];
			int x = info[4*i+1];
			map[y*npcol+x] = i;
		}
		printf("        ");
		for (i = 0; i < npcol; i++)
			printf("   %-3d   ", i);
		printf("\n");
		printf("       -");
		for (i = 0; i < npcol; i++)
			printf("---------");
		printf("-\n");
		for (i = 0; i < nprow; i++) {
			printf("   %-3d |", i);
			for (j = 0; j < npcol; j++)
				printf("   %-3d   ", map[i*npcol+j]); //print process number
			printf("|\n");
			printf("       |");
			for (j = 0; j < npcol; j++) {
				int n = map[i*npcol+j];
				printf(" %3dx%-3d ", info[4*n+2], info[4*n+3]);  //print submatrices size, v x h
			}
			printf("|\n");
		}
		printf("       -");
		for (i = 0; i < npcol; i++)
			printf("---------");
		printf("-\n\n");  /* end of this decorative table */


		free(info);
		free(map);
		free(reqs);
		free(stats);
	}
	else {  // send 4 integers to master-process
		int info[4];

		info[0] = myrow;  // proc's y
		info[1] = mycol;  // x
		info[2] = vsize;  // actual vertical size of submatrix
		info[3] = hsize;  // horizontal size
		MPI_Send(info, 4, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}
}

void config_grid(int nprocs, int *nprow, int *npcol)
{
	int _nprow = (int) sqrt(mpi_size);
	int _npcol = mpi_size / _nprow;
	while (_npcol * _nprow != mpi_size) {
		_nprow++;
		_npcol = mpi_size / _nprow;
	}
	*nprow = _nprow;
	*npcol = _npcol;
}







