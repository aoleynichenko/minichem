#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>

#include "../input/input.h"
#include "../input/basis.h"
#include "../util/util.h"
#include "../scalapack/lapacke.h"

#include <mpi.h>
#include <cblas.h>

#include "scf.h"
#include "ints.h"

static int max(int a, int b);
static int min(int a, int b);
void print_blacs_grid_info();
void config_grid(int nprocs, int *nprow, int *npcol);
void form_atom_centered_bfns(struct cart_mol *molecule, struct basis_function **bfns, int *M);
void compute_integrals();
void alloc_matrices();
void free_matrices();
void build_hcore();
void init_guess();
void orthobasis();

static int mpi_rank = -1;
static int mpi_size = -1;

/* scf parameters */
struct scf_opt scf_options;

/* task size - number of basis functions */
int M;

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

/* matrices */
double *S;      /* overlap matrix */
double *T;      /* kinetic energy matrix */
double *V;      /* potential energy matrix */
double *P;      /* density matrix */
double *Hcore;  /* core Hamiltonian */
double *F;      /* Fock matrix */
double *C;      /* vectors */
double *E;		/* energies of orbitals, F's eigenvalues */

/* explicit set of basis functions */
/* all basis functions are centered on atom at (x,y,z) */
struct basis_function *bfns;

void scf_energy(struct cart_mol *molecule)
{
	int izero = 0;
	int displ, sqsize, startx, starty;
	int i, j;

	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

	if (mpi_rank == 0) {
		printf("          ********************************\n");
		printf("          *      Parallel SCF Module     *\n");
		printf("          ********************************\n\n");
		
		mol_summary(&calc_info.molecule);
		line_separator();
	}

	// теперь мы должны понять, как распределены по узлам матрицы:
	//  - перекрывания S
	//  - плотности P
	//  - остовного гамильтониана Hcore
	//  - Фока F
	//  - матрица собственных векторов C
	// если M - размер базиса, то все используемые матрицы - M x M.
	form_atom_centered_bfns(molecule, &bfns, &M);  // create atom-centered basis set
	geom = molecule;

	// процедуру ССП начинаем с инициализации BLACS
	Cblacs_pinfo(&iam, &nprocs);
	if (nprocs != mpi_size)
		errquit("BLACS nprocs != mpi_size, what to do? (from scf/scf.c)");
    Cblacs_get(-1, 0, &ictxt);
    config_grid(nprocs, &nprow, &npcol);
    Cblacs_gridinit(&ictxt, "Row", nprow, npcol);  //номера процессов будут расти слева направо и снизу вверх
    Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);
    // расчет размеров подматриц. vb - размер блока по вертикали, hb - по горизонтали.
    // каждый поток владеет соответственными блоками F, S, Hcore, P, C
    vb = M / nprow;
    if (vb*nprow < M)
		vb++;
    hb = M / npcol;
    if (hb*npcol < M)
		hb++;
	vsize = numroc_(&M, &vb, &myrow, &izero, &nprow);
    hsize = numroc_(&M, &hb, &mycol, &izero, &npcol);
    
    // обмениваемся информацией о распределении кусков матриц MxM (например, P)
    sqsize = vsize*hsize;
    displ = 0;
    counts = (int *) malloc(sizeof(int) * mpi_size);
    displs = (int *) malloc(sizeof(int) * mpi_size);
    counts[mpi_rank] = sqsize;
    MPI_Allgather(&sqsize, 1, MPI_INT, counts, 1, MPI_INT, MPI_COMM_WORLD);
    for (i = 0; i < mpi_rank; i++)
		displ += counts[i];
    displs[mpi_rank] = displ;
    MPI_Allgather(&displ, 1, MPI_INT, displs, 1, MPI_INT, MPI_COMM_WORLD);
    
    // вывод краткой информации о распределении работы по процессам
	print_blacs_grid_info();

	alloc_matrices();
	compute_integrals();
	
	build_hcore();
	
	init_guess();
	
	orthobasis();
	
	// освобождаем ресурсы
	free_matrices();

	// финализуем подсистему BLACS
	Cblacs_gridexit(0);
}

// Ортогонализация базиса
void orthobasis()
{
	int i, j;
	char jobz = 'V', uplo = 'U';
    int info, itemp, lwork;
    int izero = 0, ione = 1;
    double *A = NULL, *Z = NULL, *W = NULL, *work = NULL;
    int descA[9], descZ[9];
    int starty = myrow * vb;
	int startx = mycol * hb;
	int ia = 1;
	int ja = 1;
	printf("%d: {%d; %d}\n", mpi_rank, ia, ja);
	
	printf("\nOverlap matrix at %d:\n", mpi_rank);
	for (i = 0; i < vsize; i++) {
		for (j = 0; j < hsize; j++)
			printf("%13.8f", S[i*hsize+j]);
		printf("\n");
	}
	printf("\n");
	MPI_Barrier(MPI_COMM_WORLD);
	

	// выделяем матрицы и заполняем их
	A = S;
	Z = (double *) malloc(hsize*vsize*sizeof(double));
	W = (double *) malloc(M*sizeof(double));
	
	// инициализация дескрипторов матриц
	itemp = max(1, vsize);
	descinit_(descA, &M, &M, &vb, &hb, &izero, &izero, &ictxt, &itemp, &info);
	descinit_(descZ, &M, &M, &vb, &hb, &izero, &izero, &ictxt, &itemp, &info);

	//printf("info = %d\n", info);
	lwork = -1;
	work = (double *) malloc(2*sizeof(double));

	MPI_Barrier(MPI_COMM_WORLD);
	pdsyev_(&jobz, &uplo, &M, A, &ia, &ja, descA, W, Z, &ia, &ja, descZ, work, &lwork, &info);
	lwork = (int) work[0];
	free(work);
	work = (double *) malloc(lwork * sizeof(double));
	if (!work)
		errquit("basis orthogonalization: work == NULL");

	pdsyev_(&jobz, &uplo, &M, A, &ia, &ja, descA, W, Z, &ia, &ja, descZ, work, &lwork, &info);
	free(work);

	MPI_Barrier(MPI_COMM_WORLD);
	if (mpi_rank == 0) {
		printf("\n");
		for (i = 0; i < M; i++)
			printf("%13.8f", W[i]);
	}
}

void eigen()
{
	char jobz = 'V', uplo = 'U';
	int i, j;
    int iam, nprocs;
    int ictxt, info, itemp, lwork;    // context?
    int n = 4;
    // разбиение процессов в прямоугольную сетку
    // по горизонтали располагаются nprow процессов, по вертикали - npcol.
    int nprow, npcol, myrow, mycol;
    int izero = 0, ione = 1, nb = 64, mb = 64, min_mn;
    int hsize, vsize;
    double *A = NULL, *Z = NULL, *W = NULL, *work = NULL;
    int descA[9], descZ[9];

	MPI_Barrier(MPI_COMM_WORLD);

	// сконфигурируем сетку такой, чтобы она была максимально близка к квадратной
	// по хорошему, надо подгонять под максимально равное распределение
	// между процессами и nb - размер матричного блока
	// nb == 1 как ни странно, организует равное распределение работы!!!
	nprow = (int) sqrt(mpi_size);
	npcol = mpi_size / nprow;
	while (npcol * nprow != mpi_size) {
		nprow++;
		npcol = mpi_size / nprow;
	}
	if (mpi_rank == 0)
		printf("BLACS process grid: nprow = %d, npcol = %d\n", nprow, npcol);

	// init blacs grid
	Cblacs_pinfo( &iam, &nprocs ) ;
	printf("iam = %d, nprocs = %d\n", iam, nprocs);
    Cblacs_get(-1, 0, &ictxt);  // ictxt == 0
    Cblacs_gridinit(&ictxt, "Row", nprow, npcol);  //номера процессов будут расти слева направо и снизу вверх
    Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);
    printf("%d: nprow = %d\n  npcol = %d\n  myrow = %d\n  mycol = %d\n", mpi_rank, nprow, npcol, myrow, mycol);

    // compute the NUMber of Rows Or Columns of a distributed matrix owned by the process
    // блоки режутся максимально близкими к квадратным
    hsize = numroc_(&n, &nb, &myrow, &izero, &nprow);
    vsize = numroc_(&n, &mb, &mycol, &izero, &npcol);

    MPI_Barrier(MPI_COMM_WORLD);
    printf("%d: nb = %d, hsize = %d, vsize = %d, elements = %d\n", mpi_rank, nb, hsize, vsize, hsize*vsize);
    MPI_Barrier(MPI_COMM_WORLD);

    if (mpi_rank == 0)
		printf("**************************************************************\n");

	// выделяем матрицы и заполняем их
	A = (double *) malloc(hsize*vsize*sizeof(double));
	Z = (double *) malloc(hsize*vsize*sizeof(double));
	W = (double *) malloc(n*sizeof(double));
	for (i = 0; i < hsize; i++)
		for (j = 0; j < vsize; j++) {
			A[j * hsize + i] = mpi_rank;
		}


	// инициализация дескрипторов матриц
	itemp = max(1, hsize);
	descinit_(descA, &n, &n, &nb, &nb, &izero, &izero, &ictxt, &itemp, &info);
	descinit_(descZ, &n, &n, &nb, &nb, &izero, &izero, &ictxt, &itemp, &info);
	printf("%d: A [ ", mpi_rank);
	for (i = 0; i < 9; i++)
		printf("%d ", descA[i]);
	printf("] ");
	printf("Z [ ");
	for (i = 0; i < 9; i++)
		printf("%d ", descZ[i]);
	printf("] \n");

	lwork = -1;
	work = (double *) malloc(2*sizeof(double));

	MPI_Barrier(MPI_COMM_WORLD);
	if (mpi_rank == 0)
		printf("Invoke pdsyev_ first time\n");
	pdsyev_(&jobz, &uplo, &n, A, &ione, &ione, descA, W, Z, &ione, &ione, descZ, work, &lwork, &info);
	if (mpi_rank == 0)
		printf("After first pdsyev_. work = [%g, %g]\n", work[0], work[1]);
	lwork = (int) work[0];
	free(work);
	work = (double *) malloc(lwork * sizeof(double));
	if (!work)
		errquit("work == NULL");

	if (mpi_rank == 0)
		printf("Start second pdsyev_\n");
	pdsyev_(&jobz, &uplo, &n, A, &ione, &ione, descA, W, Z, &ione, &ione, descZ, work, &lwork, &info);

	if (mpi_rank == 0)
		for (i = 0; i < n; i++)
			printf("%g ", W[i]);
	if (mpi_rank == 1)
		for (i = 0; i < n; i++)
			printf("%g ", W[i]);

	Cblacs_gridexit(0);
}

void scf_init()
{
	scf_options.wavefuntype = SCF_RHF;
	// print options
	scf_options.print_1eov = 0;
	scf_options.print_1eke = 0;
	scf_options.print_1epe = 0;
	scf_options.print_2eri = 0;
}

void init_guess()
{
	int i, j;
	double t = MPI_Wtime();
	
	for (i = 0; i < vsize; i++)
		for (j = 0; j < hsize; j++)
			P[hsize*i + j] = 0.0;
	
    if (mpi_rank == 0)
		printf("\nInitial guess done in %.6f sec\n", MPI_Wtime()-t);
}

void compute_integrals()
{
	int i, j, k, h, lap, nuc;
	int startx, starty;
	double t1, t2;

	if (mpi_rank == 0) {
		printf("   Atom-centered basis functions\n\n");
		for (i = 0; i < M; i++) {
			struct atom *a = bfns[i].a;
			struct cgtf *f = bfns[i].f;
			printf("%d (%.8f, %8f, %8f)\n", a->Z, a->r[0], a->r[1], a->r[2]);
			printf("  %d%15.8f%15.8f\n", f->L, f->exp[0], f->c[0]);
			for (j = 1; j < f->nprim; j++)
				printf("   %15.8f%15.8f\n", f->exp[j], f->c[j]);
		}
		printf("\n");
	}

	starty = myrow * vb;
	startx = mycol * hb;
	
	if (mpi_rank == 0) {
		printf("One-electron integrals evaluation algorithm: Obara-Saika\n");
		printf("Two-electron integrals evaluation algorithm: Obara-Saika\n");
	}
	
	t1 = MPI_Wtime();
	
	// Overlap, Kinetic & Potential matrices
	for (i = 0; i < vsize; i++)
		for (j = 0; j < hsize; j++) {
			struct basis_function *fi = &bfns[starty+i];
			struct basis_function *fj = &bfns[startx+j];
			
			S[hsize*i+j] = aoint_overlap(fi, fj);
			T[hsize*i+j] = aoint_kinetic(fi, fj);
			V[hsize*i+j] = aoint_potential(fi, fj);
		}
	
	MPI_Barrier(MPI_COMM_WORLD);
	if (mpi_rank == 0)
		printf("One-electron integrals done in %.6f sec\n", (MPI_Wtime()-t1)/1000);
	
	print_ints(bfns, M);
}

void build_hcore()
{
	int i, j, n = hsize*vsize;
	for (i = 0; i < n; i++)
		Hcore[i] = T[i] + V[i];
}

void makefock()
{
	
}

/* Creates atom-centered basis functions from molecular data.
 * Returns: vector of atom-centered functions bfns with length M.
 * This function should be executed by all processes.
 * */
void form_atom_centered_bfns(struct cart_mol *molecule, struct basis_function **bfns, int *M)
{
	int i, j;
	int K = 0;
	struct basis_function *p;

	for (i = 0; i < molecule->size; i++) {
		struct elem_info *e = searchByZ(molecule->atoms[i].Z);
		if (e) {
			for (j = 0; j < e->bas->size; j++)
				K += 2*e->bas->cgtfs[j].L + 1;  // пока и так сойдет
		}
	}
	p = (struct basis_function *) malloc(K * sizeof(struct basis_function));
	*bfns = p;
	*M = K;
	for (i = 0; i < molecule->size; i++) {
		struct elem_info *e = searchByZ(molecule->atoms[i].Z);
		if (!e)
			continue;
		for (j = 0; j < e->bas->size; j++) {  // add atom-centered bfn
			int L = e->bas->cgtfs[j].L;
			int m;
			for (m = 0; m < 2*L+1; m++) {
				(*p).m = m;
				(*p).f = &e->bas->cgtfs[j];
				(*p).a = &molecule->atoms[i];
				p++;
			}
		}
	}
}

void alloc_matrices()
{
	int bytes = hsize * vsize * sizeof(double);
	S = (double *) qalloc(bytes);
	T = (double *) qalloc(bytes);
	V = (double *) qalloc(bytes);
	Hcore = (double *) qalloc(bytes);
	F = (double *) qalloc(bytes);
	P = (double *) qalloc(M * M * sizeof(double));  /* density matrix is fully stored */
	C = (double *) qalloc(bytes);
}

void free_matrices()
{
	int bytes = hsize * vsize * sizeof(double);
	qfree(S, bytes);
	qfree(T, bytes);
	qfree(V, bytes);
	qfree(Hcore, bytes);
	qfree(F, bytes);
	qfree(P, M * M * sizeof(double));
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

static int max( int a, int b )
{
	if (a > b)
		return a;
	else
		return b;
}

static int min( int a, int b )
{
	if (a < b)
		return a;
	else
		return b;
}







