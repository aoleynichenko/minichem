/***********************************************************************
 * scf.c
 * =====
 * 
 * Hartree-Fock method -- general routines for RHF and UHF.
 * 
 * 2016-2018 Alexander Oleynichenko
 **********************************************************************/

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>

#include <omp.h>
#include <mpi.h>

#include "input.h"
#include "basis.h"
#include "util.h"
#include "sys.h"
#include "linalg.h"
#include "scf.h"
#include "aoints.h"

static int mpi_rank = -1;
static int mpi_size = -1;

/* task size - number of basis functions */
int M;

/* molecular data */
struct cart_mol *geom;
int Nelecs;
int Nalpha;
int Nbeta;


/***********************************************************************
 * scf_energy
 * 
 * "main" function of the SCF module.
 * performs SCF calculation (RHF or UHF).
 **********************************************************************/
void scf_energy(struct cart_mol *molecule)
{
	int izero = 0;
	int displ, sqsize, startx, starty;
	int i, j;
	
	// scf options -- we extract the from the RTDB
	int scf_wf_type;
	int scf_maxiter;
	int scf_direct;
	int scf_diis;
	int scf_diisbas;
	double scf_conv_dens;
	double scf_conv_en;
	int nproc;

	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
	
	// evaluate AO integrals and write them to disk
	geom = molecule;
	compute_aoints(geom);

	if (mpi_size != 1)
		errquit("Sorry! Parallel SCF module hasn't implemented yet. Please, run minichem on one node!");

	printf("          ********************************\n");
	printf("          *      Parallel SCF Module     *\n");
	printf("          ********************************\n\n");

	rtdb_get("top:nproc", &nproc);
	omp_set_num_threads(nproc);
		
	Nelecs = nalphabeta(molecule, &Nalpha, &Nbeta);
	// automatically set wavefunction type
	if (Nalpha == Nbeta) {
		rtdb_set("scf:wf", "%i", SCF_RHF);
	}
	else {
		rtdb_set("scf:wf", "%i", SCF_UHF);
	}
	
	// extract SCF options from the RTDB
	rtdb_get("scf:wf",        &scf_wf_type);
	rtdb_get("scf:maxiter",   &scf_maxiter);
	rtdb_get("scf:direct",    &scf_direct);
	rtdb_get("scf:conv_dens", &scf_conv_dens);
	rtdb_get("scf:conv_en",   &scf_conv_en);
	rtdb_get("scf:diisbas",   &scf_diisbas);
	rtdb_get("scf:diis",      &scf_diis);
	
	// print SCF options
	printf("                    SCF Parameters\n");
	printf("                    --------------\n");
	printf("               Wavefunction : %s\n", scf_wf_type == SCF_RHF ? "RHF" : "UHF");
	printf("             Max iterations : %d\n", scf_maxiter);
	printf("               Density conv : %g\n", scf_conv_dens);
	printf("                Energy conv : %g\n", scf_conv_en);
	printf("            Integral direct : %s\n", scf_direct ? "on" : "off");
	printf("                       DIIS : %s\n", scf_diis ? "on" : "off");
	if (scf_diis) {
		printf("                    diisbas : %d\n", scf_diisbas);
	}
	printf("\n");

	mol_summary(molecule);
	
	//form_atom_centered_bfns(molecule, &bfns, &shells, &M, &nshells);  // create atom-centered basis set
	M = nbfns;
	
	if (scf_wf_type == SCF_RHF)
		rhf_loop(geom, bfns, M);
	else
		uhf_loop(geom, bfns, M);
	
	// print timing
	timer_stats();
	
	// all data stored in the rtdb
	rtdb_print_meta();
}


/***********************************************************************
 * orthobasis
 * 
 * AO basis set orthogonalization. Canonical algorithm is employed.
 *
 * Arguments:
 * S   -- [input]  AO overlap matrix (square)
 * X   -- [output] transformation matrix (square)
 * dim -- [input]  dimension of matrices
 **********************************************************************/
void orthobasis(double *S, double *X, int dim)
{
	int i, j;
	int mbytes = dim * dim * sizeof(double); // N bytes in dim x dim matrix
	int vbytes = dim * sizeof(double);       // N bytes in vector
	double thresh = 1e-8;
	double t0 = MPI_Wtime();
	double *val;
	
	timer_new_entry("ortho", "Basis set orthogonalization");
	timer_start("ortho");
	
	printf("Basis orthogonalization algorithm: canonical\n");
	printf("Basis functions elimination threshold: %g\n", thresh);
	
	// Solve eigenproblem for overlap matrix S
	val = (double *) qalloc(vbytes);  // eigenvalues
	linalg_dsyev(S, val, X, dim);
	
	printf("Overlap matrix lowest eigenvalue = %.8f\n", val[0]);
	
	if (val[0] < thresh)
		errquit("negative defined overlap matrix or zero eigenvalues!");
	
	for (i = 0; i < dim; i++)
		for (j = 0; j < dim; j++) {
			X[i*dim+j] = X[i*dim+j] / sqrt(val[i]);
		}
	
	qfree(val, vbytes);
	
	timer_stop("ortho");
	printf("AO basis orthogonalization done in %.6f sec\n", MPI_Wtime()-t0);
}


/***********************************************************************
 * scf_init
 * 
 * Sets default values for the SCF options.
 ***********************************************************************/
void scf_init()
{
	rtdb_set("scf:wf",    "%i", SCF_RHF);
	rtdb_set("scf:guess", "%i", GUESS_EHT);
	rtdb_set("scf:diis",  "%i", 1);
	// print options
	rtdb_set("aoints:print:overlap",   "%i", 0);
	rtdb_set("aoints:print:kinetic",   "%i", 0);
	rtdb_set("aoints:print:potential", "%i", 0);
	rtdb_set("aoints:print:eri",       "%i", 0);
	// convergence options
	rtdb_set("scf:maxiter", "%i", 50);
	rtdb_set("scf:diisbas", "%i", 5);
	rtdb_set("scf:conv_dens", "%d", 1e-6);
	rtdb_set("scf:conv_en",   "%d", 1e-7);
	// direct/conventional scf
	rtdb_set("scf:direct", "%i", 0);  // conventional
}


/***********************************************************************
 * read_1e_integrals 
 * 
 * Reads 1e integrals from disk (they are ALWAYS stored on disk)
 ***********************************************************************/
void read_1e_integrals(double *Hcore, double *S, struct basis_function *bfns, int dim)
{
	int i, j;
	double t1;
	double *T, *V;
	int fd;
	int nbytes;
	
	nbytes = sizeof(double) * dim * dim;
	T = (double *) qalloc(nbytes);
	V = (double *) qalloc(nbytes);
	
	t1 = MPI_Wtime();
	
	// read integrals as square matrices
	fd = fastio_open("AOINTS1", "r");
	fastio_read_int(fd, &dim);
	fastio_read_doubles(fd, S, dim*dim);
	fastio_read_doubles(fd, T, dim*dim);
	fastio_read_doubles(fd, V, dim*dim);	
	fastio_close(fd);
	
	// construct Hcore = T + V
	for (i = 0; i < dim; i++)
		for (j = 0; j < dim; j++) {
			Hcore[dim*i+j] = T[dim*i+j] + V[dim*i+j];
		}
	
	printf("One-electron integrals read in %.6f sec\n", (MPI_Wtime()-t1)/1000);
	
	// TODO: remove
	print_ints(bfns, dim);
	
	// cleanup
	qfree(T, nbytes);
	qfree(V, nbytes);
}


/***********************************************************************
 * enuc
 * 
 * Calculates nuclear-repulsion energy for the given molecular geometry.
 **********************************************************************/
double enuc(Molecule_t *geom)
{
	double e = 0.0;
	double r_ij;
	int i, j;
	int n_at = geom->size;
	Atom_t *atoms = geom->atoms;
	
	for (i = 0; i < n_at; i++)
		for (j = i+1; j < n_at; j++) {
			r_ij = distance(geom->atoms[i].r, geom->atoms[j].r);
			e += atoms[i].Z * atoms[j].Z / r_ij;
		}
	
	return e;
}


/***********************************************************************
 * rmsdens
 * 
 * "root mean square" of two [density] matrices
 **********************************************************************/
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


/***********************************************************************
 * maxerr
 * 
 * returns max element of the [error] matrix.
 * 
 * TODO: move to the linalg module as a general routine.
 * can be rewritten without some knowledge about matrix representation
 **********************************************************************/
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


/***********************************************************************
 * diag_fock
 * 
 * Fock matrix diagonalization routine.
 * 
 * Agruments:
 * [input]
 *  F  -- non-transformed Fock matrix in AO basis
 *  X  -- transformation matrix (was used to orthogonalize the AO basis)
 *  M  -- dimensions of all matrices
 * [output]
 *  C  -- (MxM matrix) MO expansion coefficients (row-wise?)
 *  en -- (len=M vector) eigenvalues (orbital energies)
 **********************************************************************/
void diag_fock(double *F, double *X, double *C, double *en, int M)
{
	int i, j;
	
	double *Temp = (double *) qalloc(M*M*sizeof(double));
	double *TempF = (double *) qalloc(M*M*sizeof(double));
	
	timer_new_entry("diagfock", "Fock matrix diagonalization");
	timer_start("diagfock");
	
	memcpy(TempF, F, M*M*sizeof(double));
	
	// F = X*F*X', X stored in transposed form
	linalg_square_dgemm(TempF, 'N', X,    'T', Temp,  M);
	linalg_square_dgemm(X,     'N', Temp, 'N', TempF, M);
	
	// diagonalize F'
	linalg_dsyev(TempF, en, C, M);
	
	// transform C:  C = X*C
	linalg_square_dgemm(C, 'N', X, 'N', Temp, M);
	memcpy(C, Temp, M*M*sizeof(double));
	// Now, transformed vectors are 'lying' in the C matrix rows
	
	qfree(Temp, M*M*sizeof(double));
	qfree(TempF, M*M*sizeof(double));
	
	timer_stop("diagfock");
}


