/***********************************************************************
 * 1e.c
 * ====
 * 
 * Molecular integral evaluation module. One-electron integrals.
 * NOTE: not thread-safe!
 * 
 * For details, see, for instance,
 * T. Helgaker, P. Jorgensen, J. Olsen, "Molecular Electronic-Structure Theory".
 *
 * 2016-2018 Alexander Oleynichenko
 **********************************************************************/

#include <math.h>

#include "aoints.h"
#include "basis.h"
#include "util.h"

#ifndef M_PI
  #define M_PI 3.14159265358979323846
#endif

/* Some help functions */
double dist(double *A, double *B)
{
	return sqrt((A[0] - B[0])*(A[0] - B[0]) +
				(A[1] - B[1])*(A[1] - B[1]) +
				(A[2] - B[2])*(A[2] - B[2]));
}

double dist2(double *A, double *B)
{
	return ((A[0] - B[0])*(A[0] - B[0]) +
			(A[1] - B[1])*(A[1] - B[1]) +
			(A[2] - B[2])*(A[2] - B[2]));
}

/* Normalization */

double N0(double a)
{
	return pow((2.0*a/M_PI), 0.75);
}

double N1(double a)
{
	return 2.0*sqrt(a)*pow((2.0*a/M_PI), 0.75);
}

double N00(double a, double b)
{
	return pow(4.0*a*b/(M_PI*M_PI), 0.75);
}

double N01(double a, double b)
{
	return 2.0*sqrt(b)*pow(4.0*a*b/(M_PI*M_PI), 0.75);
}

double N10(double a, double b)
{
	return 2.0*sqrt(a)*pow(4.0*a*b/(M_PI*M_PI), 0.75);
}

double N11(double a, double b)
{
	return 4.0*sqrt(a*b)*pow(4.0*a*b/(M_PI*M_PI), 0.75);
}

/***********************************************************************
 *                              GLOBALS
 *   we use it to reduce number of parameters
 ***********************************************************************/
double alpha;
double beta;
double p, q;
double mu;
double Rpc2;
double Kab;

double Pat[3], Qat[3];
double *Aat, *Bat, *Cat, *Dat;

// geometry
extern struct cart_mol *geom;

/***********************************************************************
 ****                     OVERLAP INTEGRALS                         ****
 ***********************************************************************/
double Sij(int i, int j, double xA, double xB)
{
	if (i + j == 0)
		return sqrt(M_PI/p)*exp(-mu*(xA-xB)*(xA-xB));
	
	else if (i >= j) {
		i -= 1;
		return -beta/p*(xA-xB)*Sij(i, j, xA, xB) + 0.5/p *
			((i == 0 ? 0.0 : i*Sij(i-1, j, xA, xB)) + (j == 0 ? 0.0 : j*Sij(i, j-1, xA, xB)));
	}
	else {  // i < j
		j -= 1;
		return alpha/p*(xA-xB)*Sij(i, j, xA, xB) + 0.5/p *
			((i == 0 ? 0.0 : i*Sij(i-1, j, xA, xB)) + (j == 0 ? 0.0 : j*Sij(i, j-1, xA, xB)));
	}
}

// returns overlap integral on two primitive gaussians with exponents
// a and b, centered on atoms A and B, respectively
double Sab(int *ijklmn, double *A, double *B, double a, double b)
{
	p = a + b;
	mu = a*b/p;
	alpha = a;
	beta = b;
	
	return Sij(ijklmn[0], ijklmn[1], A[0], B[0]) *
		   Sij(ijklmn[2], ijklmn[3], A[1], B[1]) *
		   Sij(ijklmn[4], ijklmn[5], A[2], B[2]);
}

double aoint_overlap(struct basis_function *fi, struct basis_function *fj)
{
	int k, h;
	int L1 = fi->f->L;
	int L2 = fj->f->L;
	double s = 0.0;
	int set[] = {0, 0, 0, 0, 0, 0};
	double (*Norm)(double a, double b);
	
	if (L1 == BFN_S && L2 == BFN_S) {
		Norm = N00;
	}
	else if (L1 == BFN_S && L2 == BFN_P) {
		Norm = N01;
		set[2*fj->m+1] = 1;
	}
	else if (L1 == BFN_P && L2 == BFN_S) {
		Norm = N10;
		set[2*fi->m] = 1;
	}
	else if (L1 == BFN_P && L2 == BFN_P) {
		Norm = N11;
		set[2*fi->m] = 1;
		set[2*fj->m+1] = 1;
	}
	else
		errquit("angular momentum of basis function not equal to 0,1 hasn't implemented yet");
			
	for (k = 0; k < fi->f->nprim; k++)
		for (h = 0; h < fj->f->nprim; h++) {
			double a = fi->f->exp[k];
			double b = fj->f->exp[h];
			double c = Norm(a, b) * fi->f->c[k] * fj->f->c[h];
			s += c * Sab(set, fi->a->r, fj->a->r, a, b);
		}
	return s;
}

/***********************************************************************
 ****                  KINETIC-ENERGY INTEGRALS                     ****
 ***********************************************************************/

double Tij(int i, int j, double xA, double xB)
{
	double result;
	
	if (i + j == 0)
		return (alpha - 2.0*alpha*alpha*(beta*beta/(p*p)*(xA-xB)*(xA-xB) + 1.0/(2.0*p))) * Sij(i, j, xA, xB);
	else if (i >= j) {
		i -= 1;
		result = -beta/p*(xA-xB)*Tij(i, j, xA, xB) + 2.0*mu*Sij(i+1, j, xA, xB);
		if (i != 0)
			result += i*(Tij(i-1, j, xA, xB)/(2.0*p) - beta/p*Sij(i-1, j, xA, xB));
		if (j != 0)
			result += j*0.5*Tij(i, j-1, xA, xB)/p;
	}
	else {
		j -= 1;
		result = alpha/p*(xA-xB)*Tij(i, j, xA, xB) + 2.0*mu*Sij(i, j+1, xA, xB);
		if (j != 0)
			result += j*(Tij(i, j-1, xA, xB)/(2.0*p) - alpha/p*Sij(i, j-1, xA, xB));
		if (i != 0)
			result += i*0.5*Tij(i-1, j, xA, xB)/p;
	}
	return result;
}

double Tab(int *ijklmn, double *A, double *B, double a, double b)
{
	double Tabx, Taby, Tabz;
	p = a + b;
	mu = a*b/p;
	alpha = a;
	beta = b;
	
	Tabx =
		Tij(ijklmn[0], ijklmn[1], A[0], B[0]) *
		Sij(ijklmn[2], ijklmn[3], A[1], B[1]) *
		Sij(ijklmn[4], ijklmn[5], A[2], B[2]);
	Taby =
		Sij(ijklmn[0], ijklmn[1], A[0], B[0]) *
		Tij(ijklmn[2], ijklmn[3], A[1], B[1]) *
		Sij(ijklmn[4], ijklmn[5], A[2], B[2]);
	Tabz =
		Sij(ijklmn[0], ijklmn[1], A[0], B[0]) *
		Sij(ijklmn[2], ijklmn[3], A[1], B[1]) *
		Tij(ijklmn[4], ijklmn[5], A[2], B[2]);
	return Tabx + Taby + Tabz;
}

double aoint_kinetic(struct basis_function *fi, struct basis_function *fj)
{
	int k, h;
	int L1 = fi->f->L;
	int L2 = fj->f->L;
	double t = 0.0;
	int set[] = {0, 0, 0, 0, 0, 0};
	double (*Norm)(double a, double b);
	
	if (L1 == BFN_S && L2 == BFN_S) {
		Norm = N00;
	}
	else if (L1 == BFN_S && L2 == BFN_P) {
		Norm = N01;
		set[2*fj->m+1] = 1;
	}
	else if (L1 == BFN_P && L2 == BFN_S) {
		Norm = N10;
		set[2*fi->m] = 1;
	}
	else if (L1 == BFN_P && L2 == BFN_P) {
		Norm = N11;
		set[2*fi->m] = 1;
		set[2*fj->m+1] = 1;
	}
	else
		errquit("angular momentum of basis function not equal to 0,1 hasn't implemented yet");
			
	for (k = 0; k < fi->f->nprim; k++)
		for (h = 0; h < fj->f->nprim; h++) {
			double a = fi->f->exp[k];
			double b = fj->f->exp[h];
			double c = Norm(a, b) * fi->f->c[k] * fj->f->c[h];
			t += c * Tab(set, fi->a->r, fj->a->r, a, b);
		}
	return t;
}

/***********************************************************************
 ****                 POTENTIAL-ENERGY INTEGRALS                    ****
 ***********************************************************************/
double Sigma(int N, int *ijklmn)
{
	int i, j, t;
	double result = 0.0, Xpa, Xpc, Xpb;
	
	if (ijklmn[0] + ijklmn[1] + ijklmn[2] + ijklmn[3] + ijklmn[4] + ijklmn[5] == 0) {
		return 2.0*M_PI/p*Kab*boys(N, p*Rpc2);
	}
	
	// recursion
	if (ijklmn[0] != 0 || ijklmn[1] != 0) {
		i = ijklmn[0];
		j = ijklmn[1];
		t = 0;
	}	
	else if (ijklmn[2] != 0 || ijklmn[3] != 0) {
		i = ijklmn[2];
		j = ijklmn[3];
		t = 1;
	}
	else {
		i = ijklmn[4];
		j = ijklmn[5];
		t = 2;
	}
	Xpa = Pat[t] - Aat[t];
	Xpc = Pat[t] - Cat[t];
	Xpb = Pat[t] - Bat[t];
	
	if (i >= j) { // downward step by i
		i -= 1;
		ijklmn[2*t] -= 1;
		result += Xpa*Sigma(N, ijklmn);
		result -= Xpc*Sigma(N+1, ijklmn);
		if (i != 0) {
			ijklmn[2*t] -= 1;
			result += 0.5/p*i*Sigma(N, ijklmn);
			result -= 0.5/p*i*Sigma(N+1, ijklmn);
			ijklmn[2*t] += 1;
		}
		if (j != 0) {
			ijklmn[2*t+1] -= 1;
			result += 0.5/p*j*Sigma(N, ijklmn);
			result -= 0.5/p*j*Sigma(N+1, ijklmn);
			ijklmn[2*t+1] += 1;
		}
		ijklmn[2*t] += 1;
	}
	else { // downward step by j
		j -= 1;
		ijklmn[2*t+1] -= 1;
		result += Xpb*Sigma(N, ijklmn);
		result -= Xpc*Sigma(N+1, ijklmn);
		if (i != 0) {
			ijklmn[2*t] -= 1;
			result += 0.5/p*i*Sigma(N, ijklmn);
			result -= 0.5/p*i*Sigma(N+1, ijklmn);
			ijklmn[2*t] += 1;
		}
		if (j != 0) {
			ijklmn[2*t+1] -= 1;
			result += 0.5/p*j*Sigma(N, ijklmn);
			result -= 0.5/p*j*Sigma(N+1, ijklmn);
			ijklmn[2*t+1] += 1;
		}
		ijklmn[2*t+1] += 1;
	}
	return result;
}

double Vab(int *ijklmn, double *A, double *B, double *C, double a, double b)
{
	int i;
	
	p = a + b;
	mu = a*b/p;
	alpha = a;
	beta = b;
	Aat = A;
	Bat = B;
	Cat = C;
	
	for (i = 0; i < 3; i++)
		Pat[i] = (a*A[i]+b*B[i])/p;
	
	Rpc2 = dist2(Pat, C);
	Kab = exp(-mu*dist2(A, B));
	
	return Sigma(0, ijklmn);
}

double aoint_potential(struct basis_function *fi, struct basis_function *fj)
{
	int k, h, nuc;
	int L1 = fi->f->L;
	int L2 = fj->f->L;
	double v = 0.0;
	int set[] = {0, 0, 0, 0, 0, 0};
	double (*Norm)(double a, double b);
	
	if (L1 == BFN_S && L2 == BFN_S) {
		Norm = N00;
	}
	else if (L1 == BFN_S && L2 == BFN_P) {
		Norm = N01;
		set[2*fj->m+1] = 1;
	}
	else if (L1 == BFN_P && L2 == BFN_S) {
		Norm = N10;
		set[2*fi->m] = 1;
	}
	else if (L1 == BFN_P && L2 == BFN_P) {
		Norm = N11;
		set[2*fi->m] = 1;
		set[2*fj->m+1] = 1;
	}
	else
		errquit("angular momentum of basis function not equal to 0,1 hasn't implemented yet");
			
	for (k = 0; k < fi->f->nprim; k++)
		for (h = 0; h < fj->f->nprim; h++) {
			double a = fi->f->exp[k];
			double b = fj->f->exp[h];
			double c = Norm(a, b) * fi->f->c[k] * fj->f->c[h];
			for (nuc = 0; nuc < geom->size; nuc++) {
				double vab;
				vab = Vab(set, fi->a->r, fj->a->r, geom->atoms[nuc].r, a, b);
				v += -((double) geom->atoms[nuc].Z) * c * vab;
			}
		}
	return v;
}


/***********************************************************************
 ****                 MULTIPOLE-MOMENT INTEGRALS                    ****
 ***********************************************************************/
double Seij(int i, int j, int e, double xA, double xB)
{
	if (i + j + e == 0)
		return sqrt(M_PI/p)*exp(-mu*(xA-xB)*(xA-xB));
	
	if (i + j == 0) { // e != 0
		return Seij(i, j+1, e-1, xA, xB) + xB * Seij(i, j, e-1, xA, xB);
	}
	else if (i >= j) {
		i -= 1;
		return -beta/p*(xA-xB)*Seij(i, j, e, xA, xB) + 0.5/p *
			((i == 0 ? 0.0 : i*Seij(i-1, j, e, xA, xB)) +
			 (j == 0 ? 0.0 : j*Seij(i, j-1, e, xA, xB)) +
			 (e == 0 ? 0.0 : e*Seij(i, j, e-1, xA, xB)));
	}
	else {  // i < j
		j -= 1;
		return alpha/p*(xA-xB)*Seij(i, j, e, xA, xB) + 0.5/p *
			((i == 0 ? 0.0 : i*Seij(i-1, j, e, xA, xB)) +
			 (j == 0 ? 0.0 : j*Seij(i, j-1, e, xA, xB)) +
			 (e == 0 ? 0.0 : e*Seij(i, j, e-1, xA, xB)));
	}
}

double Seab(int *ijklmn, int *e, double *A, double *B, double a, double b)
{
	p = a + b;
	mu = a*b/p;
	alpha = a;
	beta = b;
	
	return Seij(ijklmn[0], ijklmn[1], e[0], A[0], B[0]) *
		   Seij(ijklmn[2], ijklmn[3], e[1], A[1], B[1]) *
		   Seij(ijklmn[4], ijklmn[5], e[2], A[2], B[2]);
}

double aoint_multipole(struct basis_function *fi, struct basis_function *fj, int *e)
{
	int k, h;
	int L1 = fi->f->L;
	int L2 = fj->f->L;
	double s = 0.0;
	int set[] = {0, 0, 0, 0, 0, 0};
	double (*Norm)(double a, double b);
	
	if (L1 == BFN_S && L2 == BFN_S) {
		Norm = N00;
	}
	else if (L1 == BFN_S && L2 == BFN_P) {
		Norm = N01;
		set[2*fj->m+1] = 1;
	}
	else if (L1 == BFN_P && L2 == BFN_S) {
		Norm = N10;
		set[2*fi->m] = 1;
	}
	else if (L1 == BFN_P && L2 == BFN_P) {
		Norm = N11;
		set[2*fi->m] = 1;
		set[2*fj->m+1] = 1;
	}
	else
		errquit("aoint_multipole: angular momentum of basis function not equal to 0,1 hasn't implemented yet");
			
	for (k = 0; k < fi->f->nprim; k++)
		for (h = 0; h < fj->f->nprim; h++) {
			double a = fi->f->exp[k];
			double b = fj->f->exp[h];
			double c = Norm(a, b) * fi->f->c[k] * fj->f->c[h];
			s += c * Seab(set, e, fi->a->r, fj->a->r, a, b);
		}
	return s;
}


