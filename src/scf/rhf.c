#include <cblas.h>
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

/* Prototypes of locally used functions */
void   rhf_loop(Molecule_t *molecule, BasisFunc_t *bfns, int M);
void   rhf_guess(double *F, double *H, double *P, double *S, double *X,
			   double *C, double *E, Molecule_t *molecule,
			   BasisFunc_t *bfns, int M);
void   rhf_makedensity(double *P, double *C, int nelec, int M);
double rhf_energy(double *P, double *F, double *H, int M);
void   rhf_makefock_direct(double *F, double *H, double *P, BasisFunc_t *bfns, int M);

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
	int diisbas; // diis subspace dimension
	double Enuc; // nuclei repulsion energy
	DIISList_t *diislist; // list with stored Fock and error matrices
	
	n = 1;  // iteration number 1
	t0 = MPI_Wtime();
	nbytes = M * M * sizeof(double);
	diisbas = scf_options.diisbas;
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
	compute_1e(H, S, bfns, M);
	
	// basis orthogonalization
	orthobasis(S, X, M);
	
	// guess
	rhf_guess(F, H, P, S, X, C, E, molecule, bfns, M);
	olde = rhf_energy(P, F, H, M) + Enuc;
	
	printf("#bfns = %d\n", M);
	printf("#eris = %d\n\n", (M*M*M*M+2*M*M*M+3*M*M+2*M)/8);
	printf(" iter.       Energy         Delta E       RMS-Dens       DIIS-Err     time\n");
	printf("----------------------------------------------------------------------------\n");
	while (1) {
		if (n > scf_options.maxiter) {
			printf("----------------------------------------------------------------------------\n");
			printf("      not converged!\n");
			errquit("no convergence of SCF equations! Try to increase scf:maxiter\n");
		}
			
		memcpy(P0, P, nbytes);  // store actual P
		
		rhf_makefock_direct(F, H, P, bfns, M);
		
		// in fact, now we have Fi and Di, used to contruct this Fi
		ErrM = (double *) qalloc(nbytes);
		make_error_matrix(ErrM, F, P, S, M);
		diiserror = maxerr(ErrM, M);
		
		// if DIIS is enabled
		if (scf_options.diis && diisbas != 0) {
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
	printf("      Time for:           sec\n");
	printf("---------------------------------\n");
	printf("  Orthogonalization  %9.3f\n", scf_timing.time_ortho);
	printf("  Initial guess      %9.3f\n", scf_timing.time_guess);
	printf("  Density matrix     %9.3f\n", scf_timing.time_dens);
	printf("  Fock matrix        %9.3f\n", scf_timing.time_fock);
	printf("  Diagonalization    %9.3f\n", scf_timing.time_diag);
	printf("---------------------------------\n\n");
	
	// Освобождение ресурсов, которые были заняты DIIS
	if (scf_options.diis && diislist)
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
	if (calc_info.out_molden_vectors) {
		int i;
		int *occ = (int *) qalloc(sizeof(int) * M);
		memset(occ, 0, sizeof(int) * M);
		for (i = 0; i < Nelecs/2; i++)
			occ[i] = 2;
		vectors_molden(molecule, C, E, occ, M, calc_info.name);
		qfree(occ, sizeof(int) * M);
	}
	
	// Анализ заселенностей
	mulliken(molecule, bfns, P, S, M);
	loewdin (molecule, bfns, P, S, M);
	// Расчет мультипольных моментов
	multipole_moments(molecule, bfns, P, M);
	
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
	double t = MPI_Wtime();
	
	printf("\nInitial guess: %s\n", scf_options.guess == GUESS_EHT ?
		"extended Huckel theory" : "bare nuclei");
	
	if (scf_options.guess == GUESS_EHT)
		guess_F_eht(F, S, bfns, M);
	else
		memcpy(F, H, M*M*sizeof(double));

	diag_fock(F, X, C, E, M);

	rhf_makedensity(P, C, Nelecs, M);  // generate initial density guess
	
	printf("Initial energy = %.8f\n", rhf_energy(P, F, H, M));
	printf("Initial guess done in %.6f sec\n\n", MPI_Wtime()-t);
	
	scf_timing.time_guess += MPI_Wtime() - t;
}

void rhf_makedensity(double *P, double *C, int nelec, int M)
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
	double t0 = MPI_Wtime();
	int maxthreads = omp_get_max_threads();
	double **Fi = (double **) qalloc(maxthreads * sizeof(double *));
	
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
	
	scf_timing.time_fock += MPI_Wtime() - t0;
}

