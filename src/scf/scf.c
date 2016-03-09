#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>

#include "../input/input.h"
#include "../input/basis.h"
#include "../util/util.h"
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

static int max(int a, int b);
static int min(int a, int b);
void print_blacs_grid_info();
void config_grid(int nprocs, int *nprow, int *npcol);
void form_atom_centered_bfns(struct cart_mol *molecule, struct basis_function **bfns, struct shell **shs, int *M, int *nshells);
void compute_integrals();
void alloc_matrices();
void free_matrices();
void build_hcore();
void init_guess();
void orthobasis();
void makefock();
void scf_loop();
void diag_fock(double *en);
void makedensity(double *P, double *C, int nelec);
double hf_energy(double *P, double *F, double *H);
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

struct {
	double time_diag;
	double time_fock;
	double time_ortho;
	double time_guess;
	double time_dens;
} scf_timing;


void scf_energy(struct cart_mol *molecule)
{
	int izero = 0;
	int displ, sqsize, startx, starty;
	int i, j;

	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

	if (mpi_size != 1)
		errquit("Sorry! Parallel SCF module hasn't implemented yet. Please, run minichem on one node!");
	if (mpi_rank == 0) {
		printf("          ********************************\n");
		printf("          *      Parallel SCF Module     *\n");
		printf("          ********************************\n\n");
		
		mol_summary(&calc_info.molecule);
		line_separator();
	}
	
	scf_timing.time_dens  = 0.0;
	scf_timing.time_diag  = 0.0;
	scf_timing.time_fock  = 0.0;
	scf_timing.time_guess = 0.0;
	scf_timing.time_ortho = 0.0;
	
	// теперь мы должны понять, как распределены по узлам матрицы:
	//  - перекрывания S
	//  - плотности P
	//  - остовного гамильтониана Hcore
	//  - Фока F
	//  - матрица собственных векторов C
	// если M - размер базиса, то все используемые матрицы - M x M.
	Nelecs = nelec(molecule);
	if (molecule->mult != 1)
		errquit("only RHF calculations can be performed");
	if ((Nelecs % 2) != 0)
		errquit("odd number of electrons! Only RHF calculations can be performed");
	geom = molecule;
	form_atom_centered_bfns(molecule, &bfns, &shells, &M, &nshells);  // create atom-centered basis set

	alloc_matrices();
	compute_integrals();
	
	orthobasis();
	init_guess();
	
	scf_loop();
	
	// освобождаем ресурсы
	free_matrices();
}

double hf_energy(double *P, double *F, double *H)
{
	int i, j;
	double E0 = 0.0;
	
	for (i = 0; i < M; i++)
		for (j = 0; j < M; j++)
			E0 += 0.5*P[i*M+j]*(H[i*M+j] + F[i*M+j]);
	return E0;
}

// Ортогонализация базиса
// Непараллельная реализация с использованием LAPACK
void orthobasis()
{
	int i, j;
	int info, lwork, lda = M;
	double thresh = 1e-8;
	double t0 = MPI_Wtime();
	FILE *f = NULL;
	char fn[150];
	double *val, *work, wkopt;
	
	printf("\n****************** BASIS ORTHOGONALIZATION *******************\n");
	printf("Algorithm: canonical\n");
	printf("Basis functions elimination threshold: %g\n", thresh);
	
	val = (double *) qalloc(sizeof(double) * M);
	memcpy(X, S, sizeof(double) * M * M);
	lwork = -1;
	dsyev_("V", "U", &M, X, &lda, val, &wkopt, &lwork, &info);
	lwork = (int) wkopt;
    work = (double*) qalloc(lwork * sizeof(double));
    dsyev_("V", "U", &M, X, &lda, val, work, &lwork, &info);
    
    if(info > 0)   // Check for convergence
		errquit("LAPACK failed to orthogonalize overlap matrix");
	qfree(work, lwork * sizeof(double));
	
	printf("Overlap matrix lowest eigenvalue = %.8f\n", val[0]);
	
	if (val[0] < thresh)
		errquit("negative defined overlap matrix or zero eigenvalues!");
	
	for (i = 0; i < M; i++)
		for (j = 0; j < M; j++) {
			X[i*M+j] = X[i*M+j] / sqrt(val[i]);
		}
	
	qfree(val, M * sizeof(double));
	
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

void init_guess()
{
	int i, j;
	double t = MPI_Wtime();
	
	printf("\nInitial guess: %s\n", scf_options.guess == GUESS_EHT ?
		"extended Huckel theory" : "bare nuclei");
	
	if (scf_options.guess == GUESS_EHT) {
		guess_F_eht(F, S, bfns, M);
	}
	else {
		memcpy(F, Hcore, M*M*sizeof(double));
	}

	diag_fock(E);

	makedensity(P, C, Nelecs);  // generate initial density guess
	
	printf("Initial energy = %.8f\n", hf_energy(P, F, Hcore));
	
    if (mpi_rank == 0)
		printf("Initial guess done in %.6f sec\n\n", MPI_Wtime()-t);
	scf_timing.time_guess += MPI_Wtime() - t;
}

void compute_integrals()
{
	int i, j, k, h, lap, nuc;
	int startx, starty;
	double t1, t2;

	starty = myrow * vb;
	startx = mycol * hb;
	
	if (mpi_rank == 0) {
		printf("One-electron integrals evaluation algorithm: Obara-Saika\n");
		printf("Two-electron integrals evaluation algorithm: Obara-Saika\n");
	}
	
	t1 = MPI_Wtime();
	
	// Overlap, Kinetic & Potential matrices
	for (i = 0; i < M; i++)
		for (j = 0; j < M; j++) {
			struct basis_function *fi = &bfns[starty+i];
			struct basis_function *fj = &bfns[startx+j];
			S[M*i+j] = aoint_overlap(fi, fj);
			Hcore[M*i+j] = aoint_kinetic(fi, fj) + aoint_potential(fi, fj);
		}
	
	MPI_Barrier(MPI_COMM_WORLD);
	if (mpi_rank == 0)
		printf("One-electron integrals done in %.6f sec\n", (MPI_Wtime()-t1)/1000);
	
	print_ints(bfns, M);
}

double distance(double *A, double *B)
{
	return sqrt((A[0] - B[0])*(A[0] - B[0]) +
				(A[1] - B[1])*(A[1] - B[1]) +
				(A[2] - B[2])*(A[2] - B[2]));
}

double enuc()
{
	double e = 0.0;
	int i, j;
	
	for (i = 0; i < geom->size; i++)
		for (j = i+1; j < geom->size; j++)
			e += geom->atoms[i].Z*geom->atoms[j].Z/distance(geom->atoms[i].r, geom->atoms[j].r);
	
	return e;
}

double rmsdens(double *P1, double *P2)
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

// текущую матрицу плотности записывает в OLDP
void save_P()
{
	memcpy(OLDP, P, M*M*sizeof(double));
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

void scf_loop()
{
	int n = 1, i, j, diisbas = scf_options.diisbas;
	double Enuc = enuc();
	double olde = hf_energy(P, F, Hcore) + Enuc;
	double hfe, deltap;
	double t0 = MPI_Wtime();
	DIISList_t *diislist = 0;
	
	printf("#bfns = %d\n", M);
	printf("#eris = %d\n\n", (M*M*M*M+2*M*M*M+3*M*M+2*M)/8);
	//printf("Enuc = %.8f\n", Enuc);
	//printf("#shells     = %d\n", nshells);
	//printf("#quartets   = %d\n", (nshells*nshells*nshells*nshells +
	//	2*nshells*nshells*nshells+3*nshells*nshells+2*nshells)/8);
	
	printf(" iter.       Energy         Delta E       RMS-Dens       DIIS-Err     time\n");
	printf("----------------------------------------------------------------------------\n");
	while (1) {
		if (n > scf_options.maxiter) {
			printf("----------------------------------------------------------------------------\n");
			printf("      not converged!\n");
			errquit("no convergence of SCF equations! Try to increase scf:maxiter\n");
		}
			
		save_P();
		
		makefock();
		
		// in fact, now we have Fi and Di, used to contruct this Fi
		double *errm = (double *) qalloc(sizeof(double)*M*M);
		double sum = 0;
		make_error_matrix(errm, F, P, S, M);
		double diiserror = maxerr(errm, M);
		
		if (scf_options.diis && diisbas != 0) {
			if (!diislist)
				diislist = newDIISList(errm, F, M);
			else
				diis_store(diislist, errm, F, M, diisbas);
			// extrapolate new F:
			diis_extrapolate(F, diislist, diisbas);
		}
		
		diag_fock(E);
		makedensity(P, C, Nelecs);
		
		deltap = rmsdens(P, OLDP);
		hfe = hf_energy(P, F, Hcore) + Enuc;
		printf("%4d%17.8f%15.8f%15.8f%15.8f%8.2f\n", n, hfe, hfe-olde, deltap, diiserror, MPI_Wtime()-t0);
		if (deltap < 1e-6)
			break;
		olde = hfe;
		n++;
	}
	printf("----------------------------------------------------------------------------\n");
	printf("          Total SCF energy =%15.8f\n", hfe);
	printf("  Nuclear repulsion energy =%15.8f\n", Enuc);
	
	printf("\n\n");
	printf("      Time for:           sec\n");
	printf("---------------------------------\n");
	printf("  Orthogonalization  %9.3f\n", scf_timing.time_ortho);
	printf("  Initial guess      %9.3f\n", scf_timing.time_guess);
	printf("  Density matrix     %9.3f\n", scf_timing.time_dens);
	printf("  Fock matrix        %9.3f\n", scf_timing.time_fock);
	printf("  Diagonalization    %9.3f\n", scf_timing.time_diag);
	printf("---------------------------------\n\n");
	
	// Вывод результатов
	// Энергии орбиталей
	printf("\n");
	printf("         Molecular Orbitals Summary\n");
	printf("  +-----+-----+----------------+----------+\n");
	printf("  | No  | Occ |     Energy     | Symmetry |\n");
	printf("  +-----+-----+----------------+----------+\n");
	for (i = 0; i < M; i++)
		printf("  | %-3d |  %1d  | %14.8f |     ?    |\n", i+1, (i < Nelecs/2) ? 2 : 0, E[i]);
	printf("  +-----+-----+----------------+----------+\n");
	
	// Вывод векторов в файл в формате Molden
	if (calc_info.out_molden_vectors) {
		int i;
		int *occ = (int) qalloc(sizeof(int) * M);
		memset(occ, 0, sizeof(int) * M);
		for (i = 0; i < Nelecs/2; i++)
			occ[i] = 2;
		vectors_molden(geom, C, E, occ, M, calc_info.name);
		qfree(occ, sizeof(int) * M);
	}
	
	// Анализ заселенностей
	mulliken(geom, bfns, P, S, M);
	loewdin (geom, bfns, P, S, M);
	// Расчет мультипольных моментов
	multipole_moments(geom, bfns, P, M);
}

#include "cblas.h"
int dgemm_(char *transa, char *transb, int *m, int *n, int *k,
	double *alpha, double *a, int *lda, double *b, int *ldb,
	double *beta, double *c, int *ldc);

void diag_fock(double *en)
{
	int i, j;
	int info, lwork, lda = M, ldb = M, ldc = M;
	double alpha = 1.0;
	double beta = 0.0;
	double *work, wkopt, t0 = MPI_Wtime();
	
	double *Temp = (double *) qalloc(M*M*sizeof(double));  // BAD!!! Alloc in critical block
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

void makedensity(double *P, double *C, int nelec)
{
	int i, j, a;
	double p, t0 = MPI_Wtime();
	
	#pragma omp parallel for private(p,j,a)
	for (i = 0; i < M; i++)
		for (j = 0; j < M; j++) {
			p = 0.0;
			for (a = 0; a < nelec/2; a++)
				p += C[a*M+i] * C[a*M+j];
			P[i*M+j] = p * 2.0;
		}
	
	scf_timing.time_dens += MPI_Wtime() - t0;
}

void makefock_naive()
{
	int i, j, k, l;
	
	for (i = 0; i < M; i++)
	for (j = 0; j < M; j++) {
		double Gij = 0.0;
		for (k = 0; k < M; k++)
		for (l = 0; l < M; l++) {
			struct basis_function *fi = &bfns[i];
			struct basis_function *fj = &bfns[j];
			struct basis_function *fk = &bfns[k];
			struct basis_function *fl = &bfns[l];
			Gij += P[k*M+l]*(aoint_eri(fi, fj, fk, fl) - 0.5*aoint_eri(fi, fk, fj, fl));
		}
		F[i*M+j] = Hcore[i*M+j] + Gij;
	}
}

/*
F(m, n) += D(p, q) * I(m, n, p, q)
F(n, m) += D(p, q) * I(n, m, p, q)
F(m, n) += D(q, p) * I(m, n, q, p)
F(n, m) += D(q, p) * I(n, m, q, p)
F(p, q) += D(m, n) * I(p, q, m, n)
F(p, q) += D(n, m) * I(p, q, n, m)
F(q, p) += D(m, n) * I(q, p, m, n)
F(q, p) += D(n, m) * I(q, p, n, m)

F(m, p) -= 0.5 * D(n, q) * I(m, n, p, q)
F(p, m) -= 0.5 * D(q, n) * I(p, q, m, n)
F(n, p) -= 0.5 * D(m, q) * I(n, m, p, q)
F(p, n) -= 0.5 * D(q, m) * I(p, q, n, m)
F(m, q) -= 0.5 * D(n, p) * I(m, n, q, p)
F(q, m) -= 0.5 * D(p, n) * I(q, p, m, n)
F(n, q) -= 0.5 * D(m, p) * I(n, m, q, p)
F(q, n) -= 0.5 * D(p, m) * I(q, p, n, m)
*/
void makefock()
{
	//int m, n, p, q;
	int m;
	int neri = 0;
	double t0 = MPI_Wtime();
	int work[] = {0,0,0,0,0,0,0,0};
	
	memcpy(F, Hcore, M*M*sizeof(double));
	
	#pragma omp parallel for schedule(dynamic,1)
	for (m = 0; m < M; m++) {
		int n, p, q;
	for (n = m; n < M; n++)
	for (p = m; p < M; p++)
	for (q = (p == m) ? n : p; q < M; q++) {
		struct basis_function *fm, *fn, *fp, *fq;
		double Int;
		
		double Dmm = P[m*M+m];
		double Dnn = P[n*M+n];
		double Dpp = P[p*M+p];
		double Dqq = P[q*M+q];
		
		double Dpq = P[p*M+q];
		double Dqp = P[q*M+p];
		double Dmn = P[m*M+n];
		double Dnm = P[n*M+m];
		
		double Dnq = P[n*M+q];
		double Dqn = P[q*M+n];
		double Dmq = P[m*M+q];
		double Dqm = P[q*M+m];
		double Dnp = P[n*M+p];
		double Dpn = P[p*M+n];
		double Dmp = P[m*M+p];
		double Dpm = P[p*M+m];
		
		fm = &bfns[m];
		fn = &bfns[n];
		fp = &bfns[p];
		fq = &bfns[q];
		
		Int = aoint_eri(fm, fn, fp, fq);
		
		//#pragma omp critical
		{
		if (m == n && m == p && m == q) {  // (mm|mm) - 1
			/*F[m*M+m] += Dmm * Int;
			F[m*M+m] -= 0.5 * Dmm * Int;*/
			#pragma omp critical
			{
			F[m*M+m] += 0.5 * Dmm * Int;
			}
		}
		else if (n == m && p == q) {  // (mm|pp) - 2
			#pragma omp critical
			{
			F[m*M+m] += Dpp * Int;
			F[p*M+p] += Dmm * Int;
				
			F[m*M+p] -= 0.5 * Dmp * Int;
			F[p*M+m] -= 0.5 * Dpm * Int;
			}
		}
		else if (m == p && n == q) { // (mn|mn) - 4
			/*
			F[m*M+n] += Dnm * Int;
			F[n*M+m] += Dnm * Int;
			F[m*M+n] += Dmn * Int;
			F[n*M+m] += Dmn * Int;
			
			F[m*M+n] -= 0.5 * Dnm * Int;
			F[n*M+m] -= 0.5 * Dmn * Int;
			F[n*M+n] -= 0.5 * Dmm * Int;
			F[m*M+m] -= 0.5 * Dnn * Int;*/
			#pragma omp critical
			{
			F[m*M+n] += Int * (Dmn + 0.5*Dnm);
			F[n*M+m] += Int * (Dnm + 0.5*Dmn);
			
			F[n*M+n] -= 0.5 * Dmm * Int;
			F[m*M+m] -= 0.5 * Dnn * Int;
			}
		}
		else if (n == m) { // (mm|pq) - 4
			#pragma omp critical
			{
			F[m*M+m] += Dpq * Int;
			F[m*M+m] += Dqp * Int;
			F[p*M+q] += Dmm * Int;
			F[q*M+p] += Dmm * Int;
			
			F[m*M+p] -= 0.5 * Dmq * Int;
			F[p*M+m] -= 0.5 * Dqm * Int;
			F[m*M+q] -= 0.5 * Dmp * Int;
			F[q*M+m] -= 0.5 * Dpm * Int;
			}
		}
		else if (p == q) {  // (mn|pp) - 4
			#pragma omp critical
			{
			F[m*M+n] += Dpp * Int;
			F[n*M+m] += Dpp * Int;
			F[p*M+p] += Dmn * Int;
			F[p*M+p] += Dnm * Int;
			
			F[m*M+p] -= 0.5 * Dnp * Int;
			F[p*M+m] -= 0.5 * Dpn * Int;
			F[n*M+p] -= 0.5 * Dmp * Int;
			F[p*M+n] -= 0.5 * Dpm * Int;
			}
		}
		else {  // (mn|pq) - 8
			#pragma omp critical
			{
			F[m*M+n] += Dpq * Int;
			F[n*M+m] += Dpq * Int;
			F[m*M+n] += Dqp * Int;
			F[n*M+m] += Dqp * Int;
			F[p*M+q] += Dmn * Int;
			F[p*M+q] += Dnm * Int;
			F[q*M+p] += Dmn * Int;
			F[q*M+p] += Dnm * Int;
				
			F[m*M+p] -= 0.5 * Dnq * Int;
			F[p*M+m] -= 0.5 * Dqn * Int;
			F[n*M+p] -= 0.5 * Dmq * Int;
			F[p*M+n] -= 0.5 * Dqm * Int;
			F[m*M+q] -= 0.5 * Dnp * Int;
			F[q*M+m] -= 0.5 * Dpn * Int;
			F[n*M+q] -= 0.5 * Dmp * Int;
			F[q*M+n] -= 0.5 * Dpm * Int;
			}
		}
	}
	}
	} // omp parallel
	
	scf_timing.time_fock += MPI_Wtime() - t0;
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

void alloc_matrices()
{
	int bytes = M * M * sizeof(double);
	E = (double *) qalloc(M * sizeof(double));
	S = (double *) qalloc(bytes);
	Hcore = (double *) qalloc(bytes);
	F = (double *) qalloc(bytes);
	X = (double *) qalloc(bytes);
	P = (double *) qalloc(bytes);
	OLDP = (double *) qalloc(bytes);
	C = (double *) qalloc(bytes);
}

void free_matrices()
{
	int bytes = M * M * sizeof(double);
	qfree(E, M * sizeof(double));
	qfree(S, bytes);
	qfree(Hcore, bytes);
	qfree(F, bytes);
	qfree(X, bytes);
	qfree(P, bytes);
	qfree(OLDP, bytes);
	qfree(C, bytes);
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







