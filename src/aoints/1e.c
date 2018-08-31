/***********************************************************************
 * 1e.c
 * ====
 * 
 * Molecular integral evaluation module. One-electron integrals.
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

/* Some helper functions */
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


/* for thread safety */
struct oneint_data {
	double alpha;
	double beta;
	double p, q;
	double mu;
	double Rpc2;
	double Kab;
	double Pat[3];
	double Qat[3];
	double *Aat;
	double *Bat;
	double *Cat;
	double *Dat;
};

// geometry
extern struct cart_mol *geom;

/***********************************************************************
 ****                     OVERLAP INTEGRALS                         ****
 ***********************************************************************/
double Sij(struct oneint_data *data, int i, int j, double xA, double xB)
{
	if (i + j == 0)
		return sqrt(M_PI/data->p)*exp(- data->mu*(xA-xB)*(xA-xB));
	
	else if (i >= j) {
		i -= 1;
		return - data->beta/data->p * (xA-xB) * Sij(data, i, j, xA, xB) + 0.5/data->p *
			((i == 0 ? 0.0 : i*Sij(data, i-1, j, xA, xB)) + (j == 0 ? 0.0 : j*Sij(data, i, j-1, xA, xB)));
	}
	else {  // i < j
		j -= 1;
		return data->alpha/data->p * (xA-xB) * Sij(data, i, j, xA, xB) + 0.5/data->p *
			((i == 0 ? 0.0 : i*Sij(data, i-1, j, xA, xB)) + (j == 0 ? 0.0 : j*Sij(data, i, j-1, xA, xB)));
	}
}

// returns overlap integral on two primitive gaussians with exponents
// a and b, centered on atoms A and B, respectively
double Sab(int *ijklmn, double *A, double *B, double a, double b)
{
	struct oneint_data data;
	data.p = a + b;
	data.mu = a*b/data.p;
	data.alpha = a;
	data.beta = b;
	
	return Sij(&data, ijklmn[0], ijklmn[1], A[0], B[0]) *
		   Sij(&data, ijklmn[2], ijklmn[3], A[1], B[1]) *
		   Sij(&data, ijklmn[4], ijklmn[5], A[2], B[2]);
}

double aoint_overlap(struct basis_function *fi, struct basis_function *fj)
{
	int k, h;
	int L1 = fi->f->L;
	int L2 = fj->f->L;
	double s = 0.0;
	int set[] = {0, 0, 0, 0, 0, 0};
	
	// X
	set[0] = fi->ijk[0];
	set[1] = fj->ijk[0];
	// Y
	set[2] = fi->ijk[1];
	set[3] = fj->ijk[1];
	// Z
	set[4] = fi->ijk[2];
	set[5] = fj->ijk[2];

	for (k = 0; k < fi->f->nprim; k++)
		for (h = 0; h < fj->f->nprim; h++) {
			double a = fi->f->exp[k];
			double b = fj->f->exp[h];
			double N = fi->norm[k] * fj->norm[h];
			double c = N * fi->f->c[k] * fj->f->c[h];
			s += c * Sab(set, fi->a->r, fj->a->r, a, b);
		}
	return s;
}

/***********************************************************************
 ****                  KINETIC-ENERGY INTEGRALS                     ****
 ***********************************************************************/

double Tij(struct oneint_data *data, int i, int j, double xA, double xB)
{
	double result;
	double alpha = data->alpha;
	double beta = data->beta;
	double p = data->p;
	double mu = data->mu;
	
	if (i + j == 0)
		return (alpha - 2.0*alpha*alpha*(beta*beta/(p*p)*(xA-xB)*(xA-xB) + 1.0/(2.0*p))) * Sij(data, i, j, xA, xB);
	else if (i >= j) {
		i -= 1;
		result = -beta/p*(xA-xB)*Tij(data, i, j, xA, xB) + 2.0*mu*Sij(data, i+1, j, xA, xB);
		if (i != 0)
			result += i*(Tij(data, i-1, j, xA, xB)/(2.0*p) - beta/p*Sij(data, i-1, j, xA, xB));
		if (j != 0)
			result += j*0.5*Tij(data, i, j-1, xA, xB)/p;
	}
	else {
		j -= 1;
		result = alpha/p*(xA-xB)*Tij(data, i, j, xA, xB) + 2.0*mu*Sij(data, i, j+1, xA, xB);
		if (j != 0)
			result += j*(Tij(data, i, j-1, xA, xB)/(2.0*p) - alpha/p*Sij(data, i, j-1, xA, xB));
		if (i != 0)
			result += i*0.5*Tij(data, i-1, j, xA, xB)/p;
	}
	return result;
}

double Tab(int *ijklmn, double *A, double *B, double a, double b)
{
	double Tabx, Taby, Tabz;
	struct oneint_data data;
	data.p = a + b;
	data.mu = a*b/data.p;
	data.alpha = a;
	data.beta = b;
	
	Tabx =
		Tij(&data, ijklmn[0], ijklmn[1], A[0], B[0]) *
		Sij(&data, ijklmn[2], ijklmn[3], A[1], B[1]) *
		Sij(&data, ijklmn[4], ijklmn[5], A[2], B[2]);
	Taby =
		Sij(&data, ijklmn[0], ijklmn[1], A[0], B[0]) *
		Tij(&data, ijklmn[2], ijklmn[3], A[1], B[1]) *
		Sij(&data, ijklmn[4], ijklmn[5], A[2], B[2]);
	Tabz =
		Sij(&data, ijklmn[0], ijklmn[1], A[0], B[0]) *
		Sij(&data, ijklmn[2], ijklmn[3], A[1], B[1]) *
		Tij(&data, ijklmn[4], ijklmn[5], A[2], B[2]);
	return Tabx + Taby + Tabz;
}

double aoint_kinetic(struct basis_function *fi, struct basis_function *fj)
{
	int k, h;
	double t = 0.0;
	int set[] = {0, 0, 0, 0, 0, 0};
	
	// X
	set[0] = fi->ijk[0];
	set[1] = fj->ijk[0];
	// Y
	set[2] = fi->ijk[1];
	set[3] = fj->ijk[1];
	// Z
	set[4] = fi->ijk[2];
	set[5] = fj->ijk[2];
	
	for (k = 0; k < fi->f->nprim; k++)
		for (h = 0; h < fj->f->nprim; h++) {
			double a = fi->f->exp[k];
			double b = fj->f->exp[h];
			double N = fi->norm[k] * fj->norm[h];
			double c = N * fi->f->c[k] * fj->f->c[h];
			t += c * Tab(set, fi->a->r, fj->a->r, a, b);
		}
	return t;
}

/***********************************************************************
 ****                 POTENTIAL-ENERGY INTEGRALS                    ****
 ***********************************************************************/
double Sigma(struct oneint_data *data, int N, int *ijklmn)
{
	int i, j, t;
	double result = 0.0, Xpa, Xpc, Xpb;
	double Rpc2 = data->Rpc2;
	double p = data->p;
	double Kab = data->Kab;
	double *Pat = data->Pat;
	double *Aat = data->Aat;
	double *Bat = data->Bat;
	double *Cat = data->Cat;
	
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
		result += Xpa*Sigma(data, N, ijklmn);
		result -= Xpc*Sigma(data, N+1, ijklmn);
		if (i != 0) {
			ijklmn[2*t] -= 1;
			result += 0.5/p*i*Sigma(data, N, ijklmn);
			result -= 0.5/p*i*Sigma(data, N+1, ijklmn);
			ijklmn[2*t] += 1;
		}
		if (j != 0) {
			ijklmn[2*t+1] -= 1;
			result += 0.5/p*j*Sigma(data, N, ijklmn);
			result -= 0.5/p*j*Sigma(data, N+1, ijklmn);
			ijklmn[2*t+1] += 1;
		}
		ijklmn[2*t] += 1;
	}
	else { // downward step by j
		j -= 1;
		ijklmn[2*t+1] -= 1;
		result += Xpb*Sigma(data, N, ijklmn);
		result -= Xpc*Sigma(data, N+1, ijklmn);
		if (i != 0) {
			ijklmn[2*t] -= 1;
			result += 0.5/p*i*Sigma(data, N, ijklmn);
			result -= 0.5/p*i*Sigma(data, N+1, ijklmn);
			ijklmn[2*t] += 1;
		}
		if (j != 0) {
			ijklmn[2*t+1] -= 1;
			result += 0.5/p*j*Sigma(data, N, ijklmn);
			result -= 0.5/p*j*Sigma(data, N+1, ijklmn);
			ijklmn[2*t+1] += 1;
		}
		ijklmn[2*t+1] += 1;
	}
	return result;
}

double Vab(int *ijklmn, double *A, double *B, double *C, double a, double b)
{
	int i;
	struct oneint_data data;
	
	data.p = a + b;
	data.mu = a*b/data.p;
	data.alpha = a;
	data.beta = b;
	data.Aat = A;
	data.Bat = B;
	data.Cat = C;
	
	for (i = 0; i < 3; i++)
		data.Pat[i] = (a*A[i]+b*B[i])/data.p;
	
	data.Rpc2 = dist2(data.Pat, C);
	data.Kab = exp(-data.mu*dist2(A, B));
	
	return Sigma(&data, 0, ijklmn);
}

double aoint_potential(struct basis_function *fi, struct basis_function *fj)
{
	int k, h, nuc;
	double v = 0.0;
	int set[] = {0, 0, 0, 0, 0, 0};
	
	// X
	set[0] = fi->ijk[0];
	set[1] = fj->ijk[0];
	// Y
	set[2] = fi->ijk[1];
	set[3] = fj->ijk[1];
	// Z
	set[4] = fi->ijk[2];
	set[5] = fj->ijk[2];
	
	for (k = 0; k < fi->f->nprim; k++)
		for (h = 0; h < fj->f->nprim; h++) {
			double a = fi->f->exp[k];
			double b = fj->f->exp[h];
			double N = fi->norm[k] * fj->norm[h];
			double c = N * fi->f->c[k] * fj->f->c[h];
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
double Seij(struct oneint_data *data, int i, int j, int e, double xA, double xB)
{
	double mu = data->mu;
	double p = data->p;
	double alpha = data->alpha;
	double beta = data->beta;
	
	if (i + j + e == 0)
		return sqrt(M_PI/p)*exp(-mu*(xA-xB)*(xA-xB));
	
	if (i + j == 0) { // e != 0
		return Seij(data, i, j+1, e-1, xA, xB) + xB * Seij(data, i, j, e-1, xA, xB);
	}
	else if (i >= j) {
		i -= 1;
		return -beta/p*(xA-xB)*Seij(data, i, j, e, xA, xB) + 0.5/p *
			((i == 0 ? 0.0 : i*Seij(data, i-1, j, e, xA, xB)) +
			 (j == 0 ? 0.0 : j*Seij(data, i, j-1, e, xA, xB)) +
			 (e == 0 ? 0.0 : e*Seij(data, i, j, e-1, xA, xB)));
	}
	else {  // i < j
		j -= 1;
		return alpha/p*(xA-xB)*Seij(data, i, j, e, xA, xB) + 0.5/p *
			((i == 0 ? 0.0 : i*Seij(data, i-1, j, e, xA, xB)) +
			 (j == 0 ? 0.0 : j*Seij(data, i, j-1, e, xA, xB)) +
			 (e == 0 ? 0.0 : e*Seij(data, i, j, e-1, xA, xB)));
	}
}

double Seab(int *ijklmn, int *e, double *A, double *B, double a, double b)
{
	struct oneint_data data;
	
	data.p = a + b;
	data.mu = a*b/data.p;
	data.alpha = a;
	data.beta = b;
	
	return Seij(&data, ijklmn[0], ijklmn[1], e[0], A[0], B[0]) *
		   Seij(&data, ijklmn[2], ijklmn[3], e[1], A[1], B[1]) *
		   Seij(&data, ijklmn[4], ijklmn[5], e[2], A[2], B[2]);
}

double aoint_multipole(struct basis_function *fi, struct basis_function *fj, int *e)
{
	int k, h;
	int L1 = fi->f->L;
	int L2 = fj->f->L;
	double s = 0.0;
	int set[] = {0, 0, 0, 0, 0, 0};
	
	// X
	set[0] = fi->ijk[0];
	set[1] = fj->ijk[0];
	// Y
	set[2] = fi->ijk[1];
	set[3] = fj->ijk[1];
	// Z
	set[4] = fi->ijk[2];
	set[5] = fj->ijk[2];
	
	for (k = 0; k < fi->f->nprim; k++)
		for (h = 0; h < fj->f->nprim; h++) {
			double a = fi->f->exp[k];
			double b = fj->f->exp[h];
			double N = fi->norm[k] * fj->norm[h];
			double c = N * fi->f->c[k] * fj->f->c[h];
			s += c * Seab(set, e, fi->a->r, fj->a->r, a, b);
		}
	return s;
}


