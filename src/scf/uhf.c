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

void uhf_loop(Molecule_t *molecule, BasisFunc_t *bfns, int M);
void uhf_guess(double *Fa, double *Fb, double *H, double *Pa, double *Pb,
			   double *S, double *X, double *Ca, double *Cb,
			   double *Ea, double *Eb, Molecule_t *molecule,
			   BasisFunc_t *bfns, int M);
double uhf_energy(double *Pa, double *Pb, double *Fa, double *Fb, double *H, int M);
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
	
	double *H;    // core Hamiltonian
	double *S;    // overlap
	double *X;    // transformation to orthonormal basis, stored by rows
	
	double *Fa, *Fb;   // Fock matrices
	double *Pa, *Pb;   // density matrices
	double *P0a,*P0b;  // stored density matrices from previous step
	double *Ca, *Cb;   // AO coefficients in MO with respect to spin part
	double *Ea, *Eb;   // orbital energies, alpha & beta
	
	n = 1;  // iteration number 1
	t0 = MPI_Wtime();
	nbytes = M * M * sizeof(double);
	Enuc = enuc(molecule);
	nalphabeta(molecule, &Nalpha, &Nbeta);
	
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
	compute_1e(H, S, bfns, M);
	
	// basis orthogonalization
	orthobasis(S, X, M);
	
	// guess
	uhf_guess(Fa, Fb, H, Pa, Pb, S, X, Ca, Cb, Ea, Eb, molecule, bfns, M);
	olde = uhf_energy(Pa, Pb, Fa, Fb, H, M) + Enuc;
	
	// scf loop
	printf("Nalpha = %d\nNbeta = %d\n", Nalpha, Nbeta);
	printf("#bfns = %d\n", M);
	printf("#eris = %d\n\n", (M*M*M*M+2*M*M*M+3*M*M+2*M)/8);
	printf(" iter.       Energy         Delta E       RMS-Dens       time\n");
	printf("---------------------------------------------------------------\n");
	while (1) {
		if (n > scf_options.maxiter) {
			printf("---------------------------------------------------------------\n");
			printf("      not converged!\n");
			errquit("no convergence of SCF equations! Try to increase scf:maxiter\n");
		}
			
		memcpy(P0a, Pa, nbytes);  // store actual P
		memcpy(P0b, Pb, nbytes);
		
		uhf_makefock_direct(Fa, Fb, H, Pa, Pb, bfns, M);
		
		diag_fock(Fa, X, Ca, Ea, M);
		diag_fock(Fb, X, Cb, Eb, M);
		uhf_makedensity(Pa, Ca, Nalpha, M);
		uhf_makedensity(Pb, Cb, Nbeta, M);
		
		deltap_a = rmsdens(Pa, P0a, M);
		deltap_b = rmsdens(Pb, P0b, M);
		hfe = uhf_energy(Pa, Pb, Fa, Fb, H, M) + Enuc;
		printf("%4d%17.8f%15.8f%15.8f%8.2f\n", n, hfe, hfe-olde, deltap_a, MPI_Wtime()-t0);
		printf("                                    %15.8f\n", deltap_b);
		if (deltap_a < 1e-6 && deltap_b < 1e-6)
			break;
		olde = hfe;
		n++;
	}
	printf("---------------------------------------------------------------\n");
	printf("          Total SCF energy =%15.8f\n", hfe);
	printf("  Nuclear repulsion energy =%15.8f\n", Enuc);
	printf("                 UHF <S^2> =%11.4f\n", uhf_S2(Nalpha, Nbeta, Ca, Cb, S, M));
	printf("               Exact <S^2> =%11.4f\n", exact_S2(Nalpha, Nbeta));
	
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
	double t = MPI_Wtime();
	
	printf("\nInitial guess: %s\n", scf_options.guess == GUESS_EHT ?
		"extended Huckel theory" : "bare nuclei");
	nalphabeta(molecule, &Na, &Nb);
	
	if (scf_options.guess == GUESS_EHT)
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
	
	printf("Initial energy = %.8f\n", uhf_energy(Pa, Pb, Fa, Fb, H, M));
	printf("Initial guess done in %.6f sec\n\n", MPI_Wtime()-t);
	
	scf_timing.time_guess += MPI_Wtime() - t;
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
	
	scf_timing.time_fock += MPI_Wtime() - t0;
}

// Calpha ---> Palpha
// nelec == Nalpha
void uhf_makedensity(double *P, double *C, int nelec, int dim)
{
	int i, j, a;
	double p, t0 = MPI_Wtime();
	
	#pragma omp parallel for private(p,j,a)
	for (i = 0; i < dim; i++)
		for (j = 0; j < dim; j++) {
			p = 0.0;
			for (a = 0; a < nelec; a++)
				p += C[a*dim+i] * C[a*dim+j];
			P[i*dim+j] = p;
		}
	
	scf_timing.time_dens += MPI_Wtime() - t0;
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





