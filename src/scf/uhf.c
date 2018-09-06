/***********************************************************************
 * uhf.c
 * =====
 * 
 * Unrestricted Hartree-Fock method.
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

void uhf_loop(Molecule_t *molecule, BasisFunc_t *bfns, int M);
void uhf_guess(double *Fa, double *Fb, double *H, double *Pa, double *Pb,
			   double *S, double *X, double *Ca, double *Cb,
			   double *Ea, double *Eb, Molecule_t *molecule,
			   BasisFunc_t *bfns, int M);
double uhf_energy(double *Pa, double *Pb, double *Fa, double *Fb, double *H, int M);
void uhf_makefock(double *Fa, double *Fb, double *H, double *Pa, double *Pb, int M);
void uhf_makefock_direct(double *Fa, double *Fb, double *H, double *Pa, double *Pb, BasisFunc_t *bfns, int M);
void uhf_makedensity(double *P, double *C, int nelec, int dim);
double exact_S2(int Nalpha, int Nbeta);
double uhf_S2(int Nalpha, int Nbeta, double *Ca, double *Cb, double *S, int dim);

double trace(double *A, double *B, int M)
{
	int i, j;
	double trace = 0.0;
	for (i = 0; i < M; i++)
		for (j = 0; j < M; j++)
			trace += A[i*M+j]*B[i*M+j];
	return trace;
}

void uhf_loop(Molecule_t *molecule, BasisFunc_t *bfns, int M)
{
	int n, i;
	int Nalpha, Nbeta;
	double t0;
	double hfe, olde, deltap_a, deltap_b;
	int nbytes;   // bytes in matrix to allocate
	double Enuc;  // nuclei repulsion energy
	double s2;    // current UHF <S^2>
	int maxiter;
	int guess;
	int direct_scf;
	int diis;
	
	double *H;    // core Hamiltonian
	double *S;    // overlap
	double *X;    // transformation to orthonormal basis, stored by rows
	
	double *Fa, *Fb;   // Fock matrices
	double *Pa, *Pb;   // density matrices
	double *P0a,*P0b;  // stored density matrices from previous step
	double *Ca, *Cb;   // AO coefficients in MO with respect to spin part
	double *Ea, *Eb;   // orbital energies, alpha & beta

	// for DIIS
	double *ErrMa, *ErrMb;  // error matrix, e=FDS-SDF
	int diisbas;            // diis subspace dimension
	DIISList_t *diislist_a, *diislist_b; // list with stored Fock and error matrices
	double diiserror_a, diiserror_b;
	
	// read scf parameters from the rtdb
	rtdb_get("scf:maxiter", &maxiter);
	rtdb_get("scf:guess",   &guess);
	rtdb_get("scf:direct",    &direct_scf);
	rtdb_get("scf:diis",    &diis);
	rtdb_get("scf:diisbas", &diisbas);
	
	n = 1;  // iteration number 1
	t0 = MPI_Wtime();
	nbytes = M * M * sizeof(double);
	Enuc = enuc(molecule);
	nalphabeta(molecule, &Nalpha, &Nbeta);
	diislist_a = NULL;
	diislist_b = NULL;
	
	H =  (double *) qalloc(nbytes);
	S =  (double *) qalloc(nbytes);
	X =  (double *) qalloc(nbytes);
	Fa = (double *) qalloc(nbytes); Fb = (double *) qalloc(nbytes);
	Pa = (double *) qalloc(nbytes); Pb = (double *) qalloc(nbytes);
	P0a= (double *) qalloc(nbytes); P0b= (double *) qalloc(nbytes);
	Ca = (double *) qalloc(nbytes); Cb = (double *) qalloc(nbytes);
	Ea = (double *) qalloc(M*sizeof(double));
	Eb = (double *) qalloc(M*sizeof(double));
	
	// compute core Hamiltonian and overlap matrices
	read_1e_integrals(H, S, bfns, M);
	
	// basis orthogonalization
	orthobasis(S, X, M);
	
	// guess
	uhf_guess(Fa, Fb, H, Pa, Pb, S, X, Ca, Cb, Ea, Eb, molecule, bfns, M);
	olde = uhf_energy(Pa, Pb, Fa, Fb, H, M) + Enuc;
	
	// scf loop
	printf("Nalpha = %d\nNbeta = %d\n\n", Nalpha, Nbeta);
	//printf("#bfns = %d\n", M);
	//printf("#eris = %d\n\n", (M*M*M*M+2*M*M*M+3*M*M+2*M)/8);
	printf(" iter.       Energy         Delta E       RMS-Dens      <S^2>       time\n");
	printf("--------------------------------------------------------------------------\n");
	while (1) {
		if (n > maxiter) {
			printf("--------------------------------------------------------------------------\n");
			printf("      not converged!\n");
			errquit("no convergence of SCF equations! Try to increase scf:maxiter\n");
		}
			
		memcpy(P0a, Pa, nbytes);  // store actual P
		memcpy(P0b, Pb, nbytes);
		
		if (direct_scf) {
			// integral-direct
			uhf_makefock_direct(Fa, Fb, H, Pa, Pb, bfns, M);
		}
		else {
			// conventional scf (2e int-s from disk)
			uhf_makefock(Fa, Fb, H, Pa, Pb, M);
		}
		
		// in fact, now we have Fi and Di, used to contruct this Fi
		// alpha
		ErrMa = (double *) qalloc(nbytes);
		make_error_matrix(ErrMa, Fa, Pa, S, M);  // Pa or P = Pa + Pb?
		diiserror_a = maxerr(ErrMa, M);
		// beta
		ErrMb = (double *) qalloc(nbytes);
		make_error_matrix(ErrMb, Fb, Pb, S, M);
		diiserror_b = maxerr(ErrMb, M);
		
		// if DIIS is enabled
		if (diis && diisbas != 0) {
			if (!diislist_a) {
				diislist_a = newDIISList(ErrMa, Fa, M);
				diislist_b = newDIISList(ErrMb, Fb, M);
			}
			else {
				diis_store(diislist_a, ErrMa, Fa, M, diisbas);
				diis_store(diislist_b, ErrMb, Fb, M, diisbas);
			}
			// extrapolate new Fa and Fb:
			diis_extrapolate(Fa, diislist_a, diisbas);
			diis_extrapolate(Fb, diislist_b, diisbas);
		}
		
		diag_fock(Fa, X, Ca, Ea, M);
		diag_fock(Fb, X, Cb, Eb, M);
		uhf_makedensity(Pa, Ca, Nalpha, M);
		uhf_makedensity(Pb, Cb, Nbeta, M);
		
		deltap_a = rmsdens(Pa, P0a, M);
		deltap_b = rmsdens(Pb, P0b, M);
		hfe = uhf_energy(Pa, Pb, Fa, Fb, H, M) + Enuc;
		s2 = uhf_S2(Nalpha, Nbeta, Ca, Cb, S, M);
		printf("%4d%17.8f%15.8f%15.8f%10.4f%11.2f\n", n, hfe, hfe-olde, deltap_a, s2, MPI_Wtime()-t0);
		printf("                                    %15.8f\n", deltap_b);
		if (deltap_a < 1e-6 && deltap_b < 1e-6)
			break;
		olde = hfe;
		n++;
	}
	printf("--------------------------------------------------------------------------\n");
	printf("          Total SCF energy =%15.8f\n", hfe);
	printf("  Nuclear repulsion energy =%15.8f\n", Enuc);
	printf("                 UHF <S^2> =%11.4f\n", uhf_S2(Nalpha, Nbeta, Ca, Cb, S, M));
	printf("               Exact <S^2> =%11.4f\n", exact_S2(Nalpha, Nbeta));
	
	/* save results to rtdb */
	rtdb_set("scf:etot", "%d", hfe);
	rtdb_set("scf:enuc", "%d", Enuc);
	rtdb_set("scf:uhf_s2", "%d", s2);
	
	// MO analysis
	printf("\n");
	printf("      Alpha Molecular Orbitals Summary\n");
	printf("  +-----+-----+----------------+----------+\n");
	printf("  | No  | Occ |     Energy     | Symmetry |\n");
	printf("  +-----+-----+----------------+----------+\n");
	for (i = 0; i < M; i++)
		printf("  | %-3d |  %1d  | %14.8f |     ?    |\n", i+1, (i < Nalpha) ? 1 : 0, Ea[i]);
	printf("  +-----+-----+----------------+----------+\n");
	printf("\n");
	printf("       Beta Molecular Orbitals Summary\n");
	printf("  +-----+-----+----------------+----------+\n");
	printf("  | No  | Occ |     Energy     | Symmetry |\n");
	printf("  +-----+-----+----------------+----------+\n");
	for (i = 0; i < M; i++)
		printf("  | %-3d |  %1d  | %14.8f |     ?    |\n", i+1, (i < Nbeta) ? 1 : 0, Eb[i]);
	printf("  +-----+-----+----------------+----------+\n");
	
	
	// cleanup DIIS
	if (diis && diislist_a) {
		removeDIISList(diislist_a);
		removeDIISList(diislist_b);
	}
	
	// cleanup
	qfree(H, nbytes);
	qfree(S, nbytes);
	qfree(X, nbytes);
	qfree(Fa, nbytes);    qfree(Fb, nbytes);
	qfree(Pa, nbytes);    qfree(Pb, nbytes);
	qfree(P0a,nbytes);    qfree(P0b,nbytes);
	qfree(Ca, nbytes);    qfree(Cb, nbytes);
	qfree(Ea, M*sizeof(double)); qfree(Eb, M*sizeof(double));
}

void uhf_guess(double *Fa, double *Fb, double *H, double *Pa, double *Pb,
			   double *S, double *X, double *Ca, double *Cb,
			   double *Ea, double *Eb, Molecule_t *molecule,
			   BasisFunc_t *bfns, int M)
{
	int Na, Nb;
	int guess_type;
	
	rtdb_get("scf:guess", &guess_type);

	timer_new_entry("guess", "Initial guess");
	timer_start("guess");
	
	printf("\nInitial guess: %s\n", guess_type == GUESS_EHT ?
		"extended Huckel theory" : "bare nuclei");
	nalphabeta(molecule, &Na, &Nb);
	
	if (guess_type == GUESS_EHT)
		guess_F_eht(Fa, S, bfns, M);
	else
		memcpy(Fa, H, M*M*sizeof(double));

	diag_fock(Fa, X, Ca, Ea, M);
	// initial Fa, Ca and Fb, Cb are equal
	memcpy(Fb, Fa, M*M*sizeof(double));
	memcpy(Cb, Ca, M*M*sizeof(double));

	// initial density matrices ARE NOT EQUAL, respect to different
	// number of alpha and beta electrons
	uhf_makedensity(Pa, Ca, Na, M);
	uhf_makedensity(Pb, Cb, Nb, M);
	
	printf("Initial energy = %.8f\n\n", uhf_energy(Pa, Pb, Fa, Fb, H, M));
	
	timer_stop("guess");
}

double uhf_energy(double *Pa, double *Pb, double *Fa, double *Fb, double *H, int M)
{
	int m, n;
	double E0 = 0.0;
	
	for (m = 0; m < M; m++)
		for (n = 0; n < M; n++)
			E0 += (Pa[n*M+m] + Pb[n*M+m])*H[m*M+n] + Pa[n*M+m]*Fa[m*M+n] + Pb[n*M+m]*Fb[m*M+n];
	return 0.5*E0;
}


/***********************************************************************
 * uhf_makefock
 * 
 * Constructs alpha and beta Fock matrices -- conventional algorithm
 * (2e integrals from disk).
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
void uhf_makefock(double *Fa, double *Fb, double *H, double *Pa, double *Pb, int M)
{
	int m, i, j, n, p, q;
	double t0 = MPI_Wtime();
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
	double *Pt = (double *) qalloc(M*M*sizeof(double));
	
	timer_new_entry("makefock", "Fock matrix construction");
	timer_start("makefock");
		
	for (m = 0; m < M*M; m++)
		Pt[m] = Pa[m] + Pb[m];
	
	memcpy(Fa, H, M*M*sizeof(double));
	memcpy(Fb, H, M*M*sizeof(double));
	
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
		
			if (m == n && m == p && m == q) {  // (mm|mm) - 1
				Fa[m*M+m] += Pt[m*M+m] * Int;
				Fb[m*M+m] += Pt[m*M+m] * Int;
				
				Fa[m*M+m] -= Pa[m*M+m] * Int;
				Fb[m*M+m] -= Pb[m*M+m] * Int;
			}
			else if (n == m && p == q) {  // (mm|pp) - 2
				Fa[m*M+m] += Pt[p*M+p] * Int;
				Fa[p*M+p] += Pt[m*M+m] * Int;
				Fb[m*M+m] += Pt[p*M+p] * Int;
				Fb[p*M+p] += Pt[m*M+m] * Int;

				Fa[m*M+p] -= Pa[m*M+p] * Int;
				Fa[p*M+m] -= Pa[p*M+m] * Int;
				Fb[m*M+p] -= Pb[m*M+p] * Int;
				Fb[p*M+m] -= Pb[p*M+m] * Int;
			}
			else if (m == p && n == q) { // (mn|mn) - 4
				Fa[m*M+n] += Pt[m*M+n] * Int;
				Fa[n*M+m] += Pt[m*M+n] * Int;
				Fa[m*M+n] += Pt[n*M+m] * Int;
				Fa[n*M+m] += Pt[n*M+m] * Int;
				
				Fb[m*M+n] += Pt[m*M+n] * Int;
				Fb[n*M+m] += Pt[m*M+n] * Int;
				Fb[m*M+n] += Pt[n*M+m] * Int;
				Fb[n*M+m] += Pt[n*M+m] * Int;
				
				// exchange
				Fa[m*M+m] -= Pa[n*M+n] * Int;
				Fa[n*M+m] -= Pa[m*M+n] * Int;
				Fa[m*M+n] -= Pa[n*M+m] * Int;
				Fa[n*M+n] -= Pa[m*M+m] * Int;
				
				Fb[m*M+m] -= Pb[n*M+n] * Int;
				Fb[n*M+m] -= Pb[m*M+n] * Int;
				Fb[m*M+n] -= Pb[n*M+m] * Int;
				Fb[n*M+n] -= Pb[m*M+m] * Int;
			}
			else if (n == m) { // (mm|pq) - 4
				Fa[m*M+m] += Pt[p*M+q] * Int;
				Fa[m*M+m] += Pt[q*M+p] * Int;
				Fa[p*M+q] += Pt[m*M+m] * Int;
				Fa[q*M+p] += Pt[m*M+m] * Int;
				
				Fb[m*M+m] += Pt[p*M+q] * Int;
				Fb[m*M+m] += Pt[q*M+p] * Int;
				Fb[p*M+q] += Pt[m*M+m] * Int;
				Fb[q*M+p] += Pt[m*M+m] * Int;
				
				Fa[m*M+p] -= Pa[m*M+q] * Int;
				Fa[p*M+m] -= Pa[q*M+m] * Int;
				Fa[m*M+q] -= Pa[m*M+p] * Int;
				Fa[q*M+m] -= Pa[p*M+m] * Int;
				
				Fb[m*M+p] -= Pb[m*M+q] * Int;
				Fb[p*M+m] -= Pb[q*M+m] * Int;
				Fb[m*M+q] -= Pb[m*M+p] * Int;
				Fb[q*M+m] -= Pb[p*M+m] * Int;
			}
			else if (p == q) {  // (mn|pp) - 4
				Fa[m*M+n] += Pt[p*M+p] * Int;
				Fa[n*M+m] += Pt[p*M+p] * Int;
				Fa[p*M+p] += Pt[m*M+n] * Int;
				Fa[p*M+p] += Pt[n*M+m] * Int;
				
				Fb[m*M+n] += Pt[p*M+p] * Int;
				Fb[n*M+m] += Pt[p*M+p] * Int;
				Fb[p*M+p] += Pt[m*M+n] * Int;
				Fb[p*M+p] += Pt[n*M+m] * Int;
				
				Fa[m*M+p] -= Pa[n*M+p] * Int;
				Fa[p*M+m] -= Pa[p*M+n] * Int;
				Fa[n*M+p] -= Pa[m*M+p] * Int;
				Fa[p*M+n] -= Pa[p*M+m] * Int;
				
				Fb[m*M+p] -= Pb[n*M+p] * Int;
				Fb[p*M+m] -= Pb[p*M+n] * Int;
				Fb[n*M+p] -= Pb[m*M+p] * Int;
				Fb[p*M+n] -= Pb[p*M+m] * Int;
			}
			else {  // (mn|pq) - 8
				// coulomb
				Fa[m*M+n] += Pt[p*M+q] * Int;
				Fa[n*M+m] += Pt[p*M+q] * Int;
				Fa[m*M+n] += Pt[q*M+p] * Int;
				Fa[n*M+m] += Pt[q*M+p] * Int;
				Fa[p*M+q] += Pt[m*M+n] * Int;
				Fa[p*M+q] += Pt[n*M+m] * Int;
				Fa[q*M+p] += Pt[m*M+n] * Int;
				Fa[q*M+p] += Pt[n*M+m] * Int;
				
				Fb[m*M+n] += Pt[p*M+q] * Int;
				Fb[n*M+m] += Pt[p*M+q] * Int;
				Fb[m*M+n] += Pt[q*M+p] * Int;
				Fb[n*M+m] += Pt[q*M+p] * Int;
				Fb[p*M+q] += Pt[m*M+n] * Int;
				Fb[p*M+q] += Pt[n*M+m] * Int;
				Fb[q*M+p] += Pt[m*M+n] * Int;
				Fb[q*M+p] += Pt[n*M+m] * Int;
				
				// exchange
				Fa[m*M+p] -= Pa[n*M+q] * Int;
				Fa[p*M+m] -= Pa[q*M+n] * Int;
				Fa[n*M+p] -= Pa[m*M+q] * Int;
				Fa[p*M+n] -= Pa[q*M+m] * Int;
				Fa[m*M+q] -= Pa[n*M+p] * Int;
				Fa[q*M+m] -= Pa[p*M+n] * Int;
				Fa[n*M+q] -= Pa[m*M+p] * Int;
				Fa[q*M+n] -= Pa[p*M+m] * Int;
				
				Fb[m*M+p] -= Pb[n*M+q] * Int;
				Fb[p*M+m] -= Pb[q*M+n] * Int;
				Fb[n*M+p] -= Pb[m*M+q] * Int;
				Fb[p*M+n] -= Pb[q*M+m] * Int;
				Fb[m*M+q] -= Pb[n*M+p] * Int;
				Fb[q*M+m] -= Pb[p*M+n] * Int;
				Fb[n*M+q] -= Pb[m*M+p] * Int;
				Fb[q*M+n] -= Pb[p*M+m] * Int;
			}
		}
	} // end loop
	
	fastio_close(fd);
	qfree(Pt, M*M*sizeof(double));
	
	timer_stop("makefock");
}


/***********************************************************************
 * uhf_makefock_direct
 * 
 * Constructs alpha and beta Fock matrices -- direct algorithm
 * (2e integrals are evaluated just-in-time).
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
void uhf_makefock_direct(double *Fa, double *Fb, double *H, double *Pa, double *Pb, BasisFunc_t *bfns, int M)
{
	int m, i, j;
	double t0 = MPI_Wtime();
	int maxthreads = omp_get_max_threads();
	double **Fai = (double **) qalloc(maxthreads * sizeof(double *));
	double **Fbi = (double **) qalloc(maxthreads * sizeof(double *));
	double *Pt = (double *) qalloc(M*M*sizeof(double));
	
	timer_new_entry("makefock", "Fock matrix construction");
	timer_start("makefock");
	
	for (i = 0; i < maxthreads; i++) {
		Fai[i] = (double *) qalloc(M*M*sizeof(double));
		Fbi[i] = (double *) qalloc(M*M*sizeof(double));
		memset(Fai[i], 0, M*M*sizeof(double));
		memset(Fbi[i], 0, M*M*sizeof(double));
	}
	
	#pragma omp parallel for
	for (m = 0; m < M*M; m++)
		Pt[m] = Pa[m] + Pb[m];
	
	#pragma omp parallel for schedule(dynamic,1)
	for (m = 0; m < M; m++) {
		int n, p, q;
		int tnum = omp_get_thread_num();
		double *fa = Fai[tnum];
		double *fb = Fbi[tnum];
	for (n = m; n < M; n++)
	for (p = m; p < M; p++)
	for (q = (p == m) ? n : p; q < M; q++) {
		struct basis_function *fm, *fn, *fp, *fq;
		double Int;
		
		fm = &bfns[m];
		fn = &bfns[n];
		fp = &bfns[p];
		fq = &bfns[q];
		Int = aoint_eri(fm, fn, fp, fq);
		
		if (fabs(Int) < 1e-14) continue;
		
		if (m == n && m == p && m == q) {  // (mm|mm) - 1
			fa[m*M+m] += Pt[m*M+m] * Int;
			fb[m*M+m] += Pt[m*M+m] * Int;
			
			fa[m*M+m] -= Pa[m*M+m] * Int;
			fb[m*M+m] -= Pb[m*M+m] * Int;
		}
		else if (n == m && p == q) {  // (mm|pp) - 2
			fa[m*M+m] += Pt[p*M+p] * Int;
			fa[p*M+p] += Pt[m*M+m] * Int;
			fb[m*M+m] += Pt[p*M+p] * Int;
			fb[p*M+p] += Pt[m*M+m] * Int;
			
			fa[m*M+p] -= Pa[m*M+p] * Int;
			fa[p*M+m] -= Pa[p*M+m] * Int;
			fb[m*M+p] -= Pb[m*M+p] * Int;
			fb[p*M+m] -= Pb[p*M+m] * Int;
		}
		else if (m == p && n == q) { // (mn|mn) - 4
			fa[m*M+n] += Pt[m*M+n] * Int;
			fa[n*M+m] += Pt[m*M+n] * Int;
			fa[m*M+n] += Pt[n*M+m] * Int;
			fa[n*M+m] += Pt[n*M+m] * Int;
			
			fb[m*M+n] += Pt[m*M+n] * Int;
			fb[n*M+m] += Pt[m*M+n] * Int;
			fb[m*M+n] += Pt[n*M+m] * Int;
			fb[n*M+m] += Pt[n*M+m] * Int;
			
			// exchange
			fa[m*M+m] -= Pa[n*M+n] * Int;
			fa[n*M+m] -= Pa[m*M+n] * Int;
			fa[m*M+n] -= Pa[n*M+m] * Int;
			fa[n*M+n] -= Pa[m*M+m] * Int;
			
			fb[m*M+m] -= Pb[n*M+n] * Int;
			fb[n*M+m] -= Pb[m*M+n] * Int;
			fb[m*M+n] -= Pb[n*M+m] * Int;
			fb[n*M+n] -= Pb[m*M+m] * Int;
		}
		else if (n == m) { // (mm|pq) - 4
			fa[m*M+m] += Pt[p*M+q] * Int;
			fa[m*M+m] += Pt[q*M+p] * Int;
			fa[p*M+q] += Pt[m*M+m] * Int;
			fa[q*M+p] += Pt[m*M+m] * Int;
			
			fb[m*M+m] += Pt[p*M+q] * Int;
			fb[m*M+m] += Pt[q*M+p] * Int;
			fb[p*M+q] += Pt[m*M+m] * Int;
			fb[q*M+p] += Pt[m*M+m] * Int;
			
			fa[m*M+p] -= Pa[m*M+q] * Int;
			fa[p*M+m] -= Pa[q*M+m] * Int;
			fa[m*M+q] -= Pa[m*M+p] * Int;
			fa[q*M+m] -= Pa[p*M+m] * Int;
			
			fb[m*M+p] -= Pb[m*M+q] * Int;
			fb[p*M+m] -= Pb[q*M+m] * Int;
			fb[m*M+q] -= Pb[m*M+p] * Int;
			fb[q*M+m] -= Pb[p*M+m] * Int;
		}
		else if (p == q) {  // (mn|pp) - 4
			fa[m*M+n] += Pt[p*M+p] * Int;
			fa[n*M+m] += Pt[p*M+p] * Int;
			fa[p*M+p] += Pt[m*M+n] * Int;
			fa[p*M+p] += Pt[n*M+m] * Int;
			
			fb[m*M+n] += Pt[p*M+p] * Int;
			fb[n*M+m] += Pt[p*M+p] * Int;
			fb[p*M+p] += Pt[m*M+n] * Int;
			fb[p*M+p] += Pt[n*M+m] * Int;
			
			fa[m*M+p] -= Pa[n*M+p] * Int;
			fa[p*M+m] -= Pa[p*M+n] * Int;
			fa[n*M+p] -= Pa[m*M+p] * Int;
			fa[p*M+n] -= Pa[p*M+m] * Int;
			
			fb[m*M+p] -= Pb[n*M+p] * Int;
			fb[p*M+m] -= Pb[p*M+n] * Int;
			fb[n*M+p] -= Pb[m*M+p] * Int;
			fb[p*M+n] -= Pb[p*M+m] * Int;
		}
		else {  // (mn|pq) - 8
			// coulomb
			fa[m*M+n] += Pt[p*M+q] * Int;
			fa[n*M+m] += Pt[p*M+q] * Int;
			fa[m*M+n] += Pt[q*M+p] * Int;
			fa[n*M+m] += Pt[q*M+p] * Int;
			fa[p*M+q] += Pt[m*M+n] * Int;
			fa[p*M+q] += Pt[n*M+m] * Int;
			fa[q*M+p] += Pt[m*M+n] * Int;
			fa[q*M+p] += Pt[n*M+m] * Int;
			
			fb[m*M+n] += Pt[p*M+q] * Int;
			fb[n*M+m] += Pt[p*M+q] * Int;
			fb[m*M+n] += Pt[q*M+p] * Int;
			fb[n*M+m] += Pt[q*M+p] * Int;
			fb[p*M+q] += Pt[m*M+n] * Int;
			fb[p*M+q] += Pt[n*M+m] * Int;
			fb[q*M+p] += Pt[m*M+n] * Int;
			fb[q*M+p] += Pt[n*M+m] * Int;
			
			// exchange
			fa[m*M+p] -= Pa[n*M+q] * Int;
			fa[p*M+m] -= Pa[q*M+n] * Int;
			fa[n*M+p] -= Pa[m*M+q] * Int;
			fa[p*M+n] -= Pa[q*M+m] * Int;
			fa[m*M+q] -= Pa[n*M+p] * Int;
			fa[q*M+m] -= Pa[p*M+n] * Int;
			fa[n*M+q] -= Pa[m*M+p] * Int;
			fa[q*M+n] -= Pa[p*M+m] * Int;
			
			fb[m*M+p] -= Pb[n*M+q] * Int;
			fb[p*M+m] -= Pb[q*M+n] * Int;
			fb[n*M+p] -= Pb[m*M+q] * Int;
			fb[p*M+n] -= Pb[q*M+m] * Int;
			fb[m*M+q] -= Pb[n*M+p] * Int;
			fb[q*M+m] -= Pb[p*M+n] * Int;
			fb[n*M+q] -= Pb[m*M+p] * Int;
			fb[q*M+n] -= Pb[p*M+m] * Int;
		}
	}
	} // end loop
	
	memcpy(Fa, H, M*M*sizeof(double));
	memcpy(Fb, H, M*M*sizeof(double));
	
	int k;
	#pragma omp parallel for private(j,k)
	for (i = 0; i < M; i++)
		for (j = 0; j < M; j++)
			for (k = 0; k < maxthreads; k++) {
				Fa[i*M+j] += Fai[k][i*M+j];
				Fb[i*M+j] += Fbi[k][i*M+j];
			}
	
	for (i = 0; i < maxthreads; i++) {
		qfree(Fai[i], M*M*sizeof(double));
		qfree(Fbi[i], M*M*sizeof(double));
	}
	qfree(Fai, maxthreads * sizeof(double *));
	qfree(Fbi, maxthreads * sizeof(double *));
	qfree(Pt, M*M*sizeof(double));
	
	timer_stop("makefock");
}

// Calpha ---> Palpha
// nelec == Nalpha
void uhf_makedensity(double *P, double *C, int nelec, int dim)
{
	int i, j, a;
	double p;
	
	timer_new_entry("dens", "Density matrix construction");
	timer_start("dens");
	
	#pragma omp parallel for private(p,j,a)
	for (i = 0; i < dim; i++)
		for (j = 0; j < dim; j++) {
			p = 0.0;
			for (a = 0; a < nelec; a++)
				p += C[a*dim+i] * C[a*dim+j];
			P[i*dim+j] = p;
		}
	
	timer_stop("dens");
}

double exact_S2(int Nalpha, int Nbeta)
{
	return 0.5*(Nalpha - Nbeta)*(0.5*(Nalpha - Nbeta) + 1);
}

double uhf_S2(int Nalpha, int Nbeta, double *Ca, double *Cb, double *S, int M)
{
	int i, j, m, n;
	double dS2 = 0.0;
	
	for (i = 0; i < Nalpha; i++)
		for (j = 0; j < Nbeta; j++) {
			double Sijab = 0.0;
			for (m = 0; m < M; m++)
				for (n = 0; n < M; n++)
					Sijab += Ca[i*M+m]*Cb[j*M+n]*S[m*M+n];
			dS2 += Sijab*Sijab;
		}
	
	return exact_S2(Nalpha, Nbeta) + Nbeta - dS2;
}





