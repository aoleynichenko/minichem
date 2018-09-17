/***********************************************************************
 * rhf.c
 * =====
 * 
 * Restricted Hartree-Fock method.
 * 
 * 2016-2018 Alexander Oleynichenko
 **********************************************************************/

#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "basis.h"
#include "input.h"
#include "aoints.h"
#include "visual.h"
#include "scf.h"
#include "util.h"
#include "sys.h"

/* Prototypes of locally used functions */
void   rhf_loop(Molecule_t *molecule, BasisFunc_t *bfns, int M);
void   rhf_guess(double *F, double *H, double *P, double *S, double *X,
			   double *C, double *E, Molecule_t *molecule,
			   BasisFunc_t *bfns, int M);
void   rhf_makedensity(double *P, double *C, int nelec, int M);
double rhf_energy(double *P, double *F, double *H, int M);
void   rhf_makefock(double *F, double *H, double *P, int M);
void   rhf_makefock_direct(double *F, double *H, double *P, BasisFunc_t *bfns, int M);
void scf_properties_rhf(struct cart_mol *molecule, double *P, int M);

void rhf_loop(Molecule_t *molecule, BasisFunc_t *bfns, int M)
{
	int i, n, Nelecs;
	double t0;
	double hfe, olde, deltap;
	double diiserror;
	double *H;   // core Hamiltonian
	double *S;   // overlap
	double *X;   // transformation to orthonormal basis, stored by rows
	double *F;   // Fock matrix
	double *P;   // density
	double *P0;  // stored density from previous step
	double *C;   // AO coefficients in MO
	double *E;   // orbital energies
	double *ErrM;// error matrix, e=FDS-SDF
	int nbytes;  // bytes in matrix to allocate
	int maxiter; // max # scf iterations
	int direct;  // enable direct scf (1/0)
	int do_diis; // is diis enabled (0/1)
	int diisbas; // diis subspace dimension
	int molden_out; // print vectors to molden format or not
	double Enuc; // nuclei repulsion energy
	DIISList_t *diislist; // list with stored Fock and error matrices
	
	// read parameters from the rtdb
	rtdb_get("scf:maxiter", &maxiter);
	rtdb_get("scf:direct",  &direct );
	rtdb_get("scf:diis",    &do_diis);
	rtdb_get("scf:diisbas", &diisbas);
	rtdb_get("visual:molden", &molden_out);
	
	n = 1;  // iteration number 1
	t0 = MPI_Wtime();
	nbytes = M * M * sizeof(double);
	diislist = NULL;
	Enuc = enuc(molecule);
	Nelecs = nelec(molecule);
	
	H = (double *) qalloc(nbytes);
	S = (double *) qalloc(nbytes);
	X = (double *) qalloc(nbytes);
	F = (double *) qalloc(nbytes);
	P = (double *) qalloc(nbytes);
	P0= (double *) qalloc(nbytes);
	C = (double *) qalloc(nbytes);
	E = (double *) qalloc(M*sizeof(double));
	
	// compute core Hamiltonian and overlap matrices
	read_1e_integrals(H, S, bfns, M);
	
	// basis orthogonalization
	orthobasis(S, X, M);
	
	// guess
	rhf_guess(F, H, P, S, X, C, E, molecule, bfns, M);
	olde = rhf_energy(P, F, H, M) + Enuc;
	
	printf(" iter.       Energy         Delta E       RMS-Dens       DIIS-Err     time\n");
	printf("----------------------------------------------------------------------------\n");
	while (1) {
		if (n > maxiter) {
			printf("----------------------------------------------------------------------------\n");
			printf("      not converged!\n");
			errquit("no convergence of SCF equations! Try to increase scf:maxiter\n");
		}
			
		memcpy(P0, P, nbytes);  // store actual P
		
		if (direct) {
			// integral-direct
			rhf_makefock_direct(F, H, P, bfns, M);
		}
		else {
			// conventional scf (2e int-s from disk)
			rhf_makefock(F, H, P, M);
		}
		
		// in fact, now we have Fi and Di, used to contruct this Fi
		ErrM = (double *) qalloc(nbytes);
		make_error_matrix(ErrM, F, P, S, M);
		diiserror = maxerr(ErrM, M);
		
		// if DIIS is enabled
		if (do_diis && diisbas != 0) {
			if (!diislist)
				diislist = newDIISList(ErrM, F, M);
			else
				diis_store(diislist, ErrM, F, M, diisbas);
			// extrapolate new F:
			diis_extrapolate(F, diislist, diisbas);
		}
		
		diag_fock(F, X, C, E, M);
		rhf_makedensity(P, C, Nelecs, M);
		
		deltap = rmsdens(P, P0, M);
		hfe = rhf_energy(P, F, H, M) + Enuc;
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
	
	/* save results to rtdb */
	rtdb_set("scf:etot", "%d", hfe);
	rtdb_set("scf:enuc", "%d", Enuc);
	
	// Освобождение ресурсов, которые были заняты DIIS
	if (do_diis && diislist)
		removeDIISList(diislist);
	
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
	if (molden_out) {
		int i;
		int *occ = (int *) qalloc(sizeof(int) * M);
		memset(occ, 0, sizeof(int) * M);
		for (i = 0; i < Nelecs/2; i++)
			occ[i] = 2;
		vectors_molden(molecule, C, E, occ, M);
		qfree(occ, sizeof(int) * M);
	}
	
	// population analysis
	mulliken(molecule, bfns, P, S, M);
	loewdin (molecule, bfns, P, S, M);

	// properties
	scf_properties_rhf(molecule, P, M);
	
	qfree(H, nbytes);
	qfree(S, nbytes);
	qfree(X, nbytes);
	qfree(F, nbytes);
	qfree(P, nbytes);
	qfree(P0,nbytes);
	qfree(C, nbytes);
	qfree(E, M*sizeof(double));
}

void rhf_guess(double *F, double *H, double *P, double *S, double *X,
			   double *C, double *E, Molecule_t *molecule,
			   BasisFunc_t *bfns, int M)
{
	int Nelecs = nelec(molecule);
	int guess_type;
	
	rtdb_get("scf:guess", &guess_type);
	
	printf("\nInitial guess: %s\n", guess_type == GUESS_EHT ?
		"extended Huckel theory" : "bare nuclei");
	
	timer_new_entry("guess", "Initial guess");
	timer_start("guess");
	
	if (guess_type == GUESS_EHT)
		guess_F_eht(F, S, bfns, M);
	else
		memcpy(F, H, M*M*sizeof(double));

	diag_fock(F, X, C, E, M);

	rhf_makedensity(P, C, Nelecs, M);  // generate initial density guess
	
	timer_stop("guess");
	printf("Initial energy = %.8f\n\n", rhf_energy(P, F, H, M));
}

void rhf_makedensity(double *P, double *C, int nelec, int M)
{
	int i, j, a;
	double p;
	
	timer_new_entry("dens", "Density matrix construction");
	timer_start("dens");
	
	#pragma omp parallel for private(p,j,a)
	for (i = 0; i < M; i++)
		for (j = 0; j < M; j++) {
			p = 0.0;
			for (a = 0; a < nelec/2; a++)
				p += C[a*M+i] * C[a*M+j];
			P[i*M+j] = p * 2.0;
		}

	timer_stop("dens");
}

double rhf_energy(double *P, double *F, double *H, int M)
{
	int i, j;
	double E0 = 0.0;
	
	for (i = 0; i < M; i++)
		for (j = 0; j < M; j++)
			E0 += 0.5*P[i*M+j]*(H[i*M+j] + F[i*M+j]);
	return E0;
}


/***********************************************************************
 * scf_properties_rhf
 * 
 * wrapper for the 'scf_properties' function.
 * transforms density matrix into [equal] alpha and beta components.
 **********************************************************************/
void scf_properties_rhf(struct cart_mol *molecule, double *P, int M)
{
	double *Pa, *Pb;
	int nbytes = M * M * sizeof(double);
	int i;
	
	Pa = (double *) qalloc(nbytes);
	Pb = (double *) qalloc(nbytes);
	
	// RHF: Pa = Pb = P/2
	for (i = 0; i < M * M; i++) {
		Pa[i] = P[i] * 0.5;
		Pb[i] = P[i] * 0.5;
	}
	
	scf_properties(molecule, Pa, Pb, M);
	
	qfree(Pa, nbytes);
	qfree(Pb, nbytes);
}


/***********************************************************************
 * rhf_makefock
 * 
 * Constructs Fock matrix; 2e integrals are stored on disk.
 **********************************************************************/
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
void rhf_makefock(double *F, double *H, double *P, int M)
{
	int m, i, j;
	int n, p, q;
	int fd;
	typedef struct {
		double val;
		int i1, i2, i3, i4;
	} integral_t;
	integral_t *tmpi;
	const int BATCH_SIZE = 4096;
	integral_t buf[BATCH_SIZE];
	int n_read;  // bytes
	int n_int_read;
	double Int;
	
	timer_new_entry("makefock", "Fock matrix construction");
	timer_start("makefock");

	memcpy(F, H, M*M*sizeof(double));
	
	fd = fastio_open("AOINTS2", "r");
	while ((n_read = fastio_read(fd, buf, sizeof(integral_t)*BATCH_SIZE)) > 0) {
		n_int_read = n_read / sizeof(integral_t);
		
		for (i = 0; i < n_int_read; i++) {
			tmpi = &buf[i];
			Int = tmpi->val;
			m = tmpi->i1;
			n = tmpi->i2;
			p = tmpi->i3;
			q = tmpi->i4;
			
			if (fabs(Int) < 1e-14) continue;
			
			double Dmm = P[m*M+m];
			double Dnn = P[n*M+n];
			double Dpp = P[p*M+p];
			//double Dqq = P[q*M+q];
			
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
			
			if (m == n && m == p && m == q) {  // (mm|mm) - 1
				F[m*M+m] += 0.5 * Dmm * Int;
			}
			else if (n == m && p == q) {  // (mm|pp) - 2
				F[m*M+m] += Dpp * Int;
				F[p*M+p] += Dmm * Int;
					
				F[m*M+p] -= 0.5 * Dmp * Int;
				F[p*M+m] -= 0.5 * Dpm * Int;
			}
			else if (m == p && n == q) { // (mn|mn) - 4
				F[m*M+n] += Int * (Dmn + 0.5*Dnm);
				F[n*M+m] += Int * (Dnm + 0.5*Dmn);
				
				F[n*M+n] -= 0.5 * Dmm * Int;
				F[m*M+m] -= 0.5 * Dnn * Int;
			}
			else if (n == m) { // (mm|pq) - 4
				F[m*M+m] += Dpq * Int;
				F[m*M+m] += Dqp * Int;
				F[p*M+q] += Dmm * Int;
				F[q*M+p] += Dmm * Int;
				
				F[m*M+p] -= 0.5 * Dmq * Int;
				F[p*M+m] -= 0.5 * Dqm * Int;
				F[m*M+q] -= 0.5 * Dmp * Int;
				F[q*M+m] -= 0.5 * Dpm * Int;
			}
			else if (p == q) {  // (mn|pp) - 4
				F[m*M+n] += Dpp * Int;
				F[n*M+m] += Dpp * Int;
				F[p*M+p] += Dmn * Int;
				F[p*M+p] += Dnm * Int;
				
				F[m*M+p] -= 0.5 * Dnp * Int;
				F[p*M+m] -= 0.5 * Dpn * Int;
				F[n*M+p] -= 0.5 * Dmp * Int;
				F[p*M+n] -= 0.5 * Dpm * Int;
			}
			else {  // (mn|pq) - 8
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
	
	
	fastio_close(fd);
	
	timer_stop("makefock");
}


/***********************************************************************
 * rhf_makefock_direct
 * 
 * Constructs Fock matrix -- direct algorithm (2e integrals are
 * evaluated just-in-time).
 **********************************************************************/
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
void rhf_makefock_direct(double *F, double *H, double *P, BasisFunc_t *bfns, int M)
{
	int m, i, j;
	int maxthreads = omp_get_max_threads();
	double **Fi = (double **) qalloc(maxthreads * sizeof(double *));
	
	timer_new_entry("makefock", "Fock matrix construction");
	timer_start("makefock");
	
	for (i = 0; i < maxthreads; i++) {
		Fi[i] = (double *) qalloc(M*M*sizeof(double));
		memset(Fi[i], 0, M*M*sizeof(double));
	}
	
	#pragma omp parallel for schedule(dynamic,1)
	for (m = 0; m < M; m++) {
		int n, p, q;
		int tnum = omp_get_thread_num();
		double *f = Fi[tnum];
	for (n = m; n < M; n++)
	for (p = m; p < M; p++)
	for (q = (p == m) ? n : p; q < M; q++) {
		struct basis_function *fm, *fn, *fp, *fq;
		double Int;
		
		double Dmm = P[m*M+m];
		double Dnn = P[n*M+n];
		double Dpp = P[p*M+p];
		//double Dqq = P[q*M+q];
		
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
		
		if (m == n && m == p && m == q) {  // (mm|mm) - 1
			f[m*M+m] += 0.5 * Dmm * Int;
		}
		else if (n == m && p == q) {  // (mm|pp) - 2
			f[m*M+m] += Dpp * Int;
			f[p*M+p] += Dmm * Int;
				
			f[m*M+p] -= 0.5 * Dmp * Int;
			f[p*M+m] -= 0.5 * Dpm * Int;
		}
		else if (m == p && n == q) { // (mn|mn) - 4
			f[m*M+n] += Int * (Dmn + 0.5*Dnm);
			f[n*M+m] += Int * (Dnm + 0.5*Dmn);
			
			f[n*M+n] -= 0.5 * Dmm * Int;
			f[m*M+m] -= 0.5 * Dnn * Int;
		}
		else if (n == m) { // (mm|pq) - 4
			f[m*M+m] += Dpq * Int;
			f[m*M+m] += Dqp * Int;
			f[p*M+q] += Dmm * Int;
			f[q*M+p] += Dmm * Int;
			
			f[m*M+p] -= 0.5 * Dmq * Int;
			f[p*M+m] -= 0.5 * Dqm * Int;
			f[m*M+q] -= 0.5 * Dmp * Int;
			f[q*M+m] -= 0.5 * Dpm * Int;
		}
		else if (p == q) {  // (mn|pp) - 4
			f[m*M+n] += Dpp * Int;
			f[n*M+m] += Dpp * Int;
			f[p*M+p] += Dmn * Int;
			f[p*M+p] += Dnm * Int;
			
			f[m*M+p] -= 0.5 * Dnp * Int;
			f[p*M+m] -= 0.5 * Dpn * Int;
			f[n*M+p] -= 0.5 * Dmp * Int;
			f[p*M+n] -= 0.5 * Dpm * Int;
		}
		else {  // (mn|pq) - 8
			f[m*M+n] += Dpq * Int;
			f[n*M+m] += Dpq * Int;
			f[m*M+n] += Dqp * Int;
			f[n*M+m] += Dqp * Int;
			f[p*M+q] += Dmn * Int;
			f[p*M+q] += Dnm * Int;
			f[q*M+p] += Dmn * Int;
			f[q*M+p] += Dnm * Int;
				
			f[m*M+p] -= 0.5 * Dnq * Int;
			f[p*M+m] -= 0.5 * Dqn * Int;
			f[n*M+p] -= 0.5 * Dmq * Int;
			f[p*M+n] -= 0.5 * Dqm * Int;
			f[m*M+q] -= 0.5 * Dnp * Int;
			f[q*M+m] -= 0.5 * Dpn * Int;
			f[n*M+q] -= 0.5 * Dmp * Int;
			f[q*M+n] -= 0.5 * Dpm * Int;
		}
	}
	} // end loop
	
	memcpy(F, H, M*M*sizeof(double));
	
	int k;
	#pragma omp parallel for private(j,k)
	for (i = 0; i < M; i++)
		for (j = 0; j < M; j++)
			for (k = 0; k < maxthreads; k++)
				F[i*M+j] += Fi[k][i*M+j];
	
	for (i = 0; i < maxthreads; i++)
		qfree(Fi[i], M*M*sizeof(double));
	qfree(Fi, maxthreads * sizeof(double *));

	timer_stop("makefock");
}

