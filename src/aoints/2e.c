/***********************************************************************
 * 1e.c
 * ====
 * 
 * Molecular integral evaluation module. Two-electron integrals.
 * Algorithm: Obara-Saika
 * 
 * For details, see, for instance,
 * T. Helgaker, P. Jorgensen, J. Olsen, "Molecular Electronic-Structure Theory".
 *
 * 2016-2018 Alexander Oleynichenko
 **********************************************************************/
 
#include <math.h>

#include "aoints.h"

#ifndef M_PI
  #define M_PI 3.14159265358979323846
#endif

/***********************************************************************
 *                              GLOBALS
 * see 1e.c for definitions
 ***********************************************************************/

// for thread safety
struct eri_data {
	double alpha;
	double p;
	double q;
	double Kab_xyz;
	double Kcd_xyz;
	double *Aat;
	double *Bat;
	double *Cat;
	double *Dat;
	double Pat[3];
	double Qat[3];
};


double Theta(struct eri_data *data, int N, int *ijkl)
{
	int i, j, k, l, t;
	double result = 0.0;
	double Xpa, Xpq, Xpb, Xqc, Xqd;
	
	double alpha = data->alpha;
	double p = data->p, q = data->q;

	double *Pat = data->Pat, *Qat = data->Qat;
	double *Aat = data->Aat, *Bat = data->Bat, *Cat = data->Cat, *Dat = data->Dat;
	
	if (ijkl[0] + ijkl[1] + ijkl[2] + ijkl[3] + ijkl[4] + ijkl[5] +
		ijkl[6] + ijkl[7] + ijkl[8] + ijkl[9] + ijkl[10] + ijkl[11] == 0) {
		return 2.0*pow(M_PI,2.5)/(p*q*sqrt(p+q)) * data->Kab_xyz * data->Kcd_xyz * boys(N, alpha*dist2(Pat, Qat));
	}
	
	if (ijkl[0] + ijkl[1] + ijkl[2] + ijkl[3] != 0) {
		i = ijkl[0];
		j = ijkl[1];
		k = ijkl[2];
		l = ijkl[3];
		t = 0;
	}
	else if (ijkl[4] + ijkl[5] + ijkl[6] + ijkl[7] != 0) {
		i = ijkl[4];
		j = ijkl[5];
		k = ijkl[6];
		l = ijkl[7];
		t = 1;
	}
	else {
		i = ijkl[8];
		j = ijkl[9];
		k = ijkl[10];
		l = ijkl[11];
		t = 2;
	}
	
	Xpa = Pat[t] - Aat[t];
	Xpq = Pat[t] - Qat[t];
	Xpb = Pat[t] - Bat[t];
	Xqc = Qat[t] - Cat[t];
	Xqd = Qat[t] - Dat[t];
	
	if (i != 0) {
		i -= 1;
		ijkl[4*t] -= 1;
		result += Xpa*Theta(data, N, ijkl);
		result += -alpha/p*Xpq*Theta(data, N+1, ijkl);
		if (i != 0) {
			ijkl[4*t] -= 1;
			result += 0.5*i/p*(Theta(data, N, ijkl) - alpha/p*Theta(data, N+1, ijkl));
			ijkl[4*t] += 1;
		}
		if (j != 0) {
			ijkl[4*t+1] -= 1;
			result += 0.5*j/p*(Theta(data, N, ijkl) - alpha/p*Theta(data, N+1, ijkl));
			ijkl[4*t+1] += 1;
		}
		if (k != 0) {
			ijkl[4*t+2] -= 1;
			result += 0.5*k/(p+q) * Theta(data, N+1, ijkl);
			ijkl[4*t+2] += 1;
		}
		if (l != 0) {
			ijkl[4*t+3] -= 1;
			result += 0.5*l/(p+q) * Theta(data, N+1, ijkl);
			ijkl[4*t+3] += 1;
		}
		ijkl[4*t] += 1;
	}
	else if (j != 0) {
		j -= 1;
		ijkl[4*t+1] -= 1;
		result += Xpb*Theta(data, N, ijkl);
		result += -alpha/p*Xpq*Theta(data, N+1, ijkl);
		if (i != 0) {
			ijkl[4*t] -= 1;
			result += 0.5*i/p*(Theta(data, N, ijkl) - alpha/p*Theta(data, N+1, ijkl));
			ijkl[4*t] += 1;
		}
		if (j != 0) {
			ijkl[4*t+1] -= 1;
			result += 0.5*j/p*(Theta(data, N, ijkl) - alpha/p*Theta(data, N+1, ijkl));
			ijkl[4*t+1] += 1;
		}
		if (k != 0) {
			ijkl[4*t+2] -= 1;
			result += 0.5*k/(p+q) * Theta(data, N+1, ijkl);
			ijkl[4*t+2] += 1;
		}
		if (l != 0) {
			ijkl[4*t+3] -= 1;
			result += 0.5*l/(p+q) * Theta(data, N+1, ijkl);
			ijkl[4*t+3] += 1;
		}
		ijkl[4*t+1] += 1;
	}
	else if (k != 0) {
		k -= 1;
		ijkl[4*t+2] -= 1;
		result += Xqc*Theta(data, N, ijkl);
		result += alpha/q*Xpq*Theta(data, N+1, ijkl);
		if (k != 0) {
			ijkl[4*t+2] -= 1;
			result += 0.5*k/q*(Theta(data, N, ijkl) - alpha/q*Theta(data, N+1, ijkl));
			ijkl[4*t+2] += 1;
		}
		if (l != 0) {
			ijkl[4*t+3] -= 1;
			result += 0.5*l/q*(Theta(data, N, ijkl) - alpha/q*Theta(data, N+1, ijkl));
			ijkl[4*t+3] += 1;
		}
		if (i != 0) {
			ijkl[4*t] -= 1;
			result += 0.5*i/(p+q) * Theta(data, N+1, ijkl);
			ijkl[4*t] += 1;
		}
		if (j != 0) {
			ijkl[4*t+1] -= 1;
			result += 0.5*j/(p+q) * Theta(data, N+1, ijkl);
			ijkl[4*t+1] += 1;
		}
		ijkl[4*t+2] += 1;
	}
	else {  // l != 0
		l -= 1;
		ijkl[4*t+3] -= 1;
		result += Xqd*Theta(data, N, ijkl);
		result += alpha/q*Xpq*Theta(data, N+1, ijkl);
		if (k != 0) {
			ijkl[4*t+2] -= 1;
			result += 0.5*k/q*(Theta(data, N, ijkl) - alpha/q*Theta(data, N+1, ijkl));
			ijkl[4*t+2] += 1;
		}
		if (l != 0) {
			ijkl[4*t+3] -= 1;
			result += 0.5*l/q*(Theta(data, N, ijkl) - alpha/q*Theta(data, N+1, ijkl));
			ijkl[4*t+3] += 1;
		}
		if (i != 0) {
			ijkl[4*t] -= 1;
			result += 0.5*i/(p+q) * Theta(data, N+1, ijkl);
			ijkl[4*t] += 1;
		}
		if (j != 0) {
			ijkl[4*t+1] -= 1;
			result += 0.5*j/(p+q) * Theta(data, N+1, ijkl);
			ijkl[4*t+1] += 1;
		}
		ijkl[4*t+3] += 1;
	}
	
	return result;
}

double ERI(int *ijkl, double *A, double *B, double *C, double *D,
					  double a, double b, double c, double d)
{
	int i;
	struct eri_data data;
	
	data.p = a + b;
	data.q = c + d;
	data.alpha = data.p*data.q/(data.p+data.q);
	data.Kab_xyz = exp(-a*b*dist2(A, B)/data.p);
	data.Kcd_xyz = exp(-c*d*dist2(C, D)/data.q);
	data.Aat = A;
	data.Bat = B;
	data.Cat = C;
	data.Dat = D;
	
	for (i = 0; i < 3; i++) {
		data.Pat[i] = (a*A[i]+b*B[i])/data.p;
		data.Qat[i] = (c*C[i]+d*D[i])/data.q;
	}
	
	return Theta(&data, 0, ijkl);
}

double aoint_eri(struct basis_function *fi,
				 struct basis_function *fj,
				 struct basis_function *fk,
				 struct basis_function *fl)
{
	double eri = 0.0;
	int a, b, c, d, Li, Lj, Lk, Ll;
	int inprim = fi->f->nprim;
	int jnprim = fj->f->nprim;
	int knprim = fk->f->nprim;
	int lnprim = fl->f->nprim;
	                   /* i  j  k  l  */
	int set[] = { /* X */ 0, 0, 0, 0, /* Y */ 0, 0, 0, 0, /* Z */ 0, 0, 0, 0};
	
	// X
	set[0] = fi->ijk[0];
	set[1] = fj->ijk[0];
	set[2] = fk->ijk[0];
	set[3] = fl->ijk[0];
	// Y
	set[4] = fi->ijk[1];
	set[5] = fj->ijk[1];
	set[6] = fk->ijk[1];
	set[7] = fl->ijk[1];
	// Z
	set[8] = fi->ijk[2];
	set[9] = fj->ijk[2];
	set[10] = fk->ijk[2];
	set[11] = fl->ijk[2];
	
	for (a = 0; a < inprim; a++)
		for (b = 0; b < jnprim; b++)
			for (c = 0; c < knprim; c++)
				for (d = 0; d < lnprim; d++) {
					double expa = fi->f->exp[a];
					double expb = fj->f->exp[b];
					double expc = fk->f->exp[c];
					double expd = fl->f->exp[d];
					/*double N = cart_norm(expa, fi->ijk) *
					           cart_norm(expb, fj->ijk) *
					           cart_norm(expc, fk->ijk) *
					           cart_norm(expd, fl->ijk);*/
					double N = fi->norm[a] * fj->norm[b] * fk->norm[c] * fl->norm[d];
					double coeff = fi->f->c[a] * fj->f->c[b] * fk->f->c[c] * fl->f->c[d];
					eri += N * coeff * ERI(set, fi->a->r, fj->a->r, fk->a->r, fl->a->r,
						expa, expb, expc, expd);
				}
	return eri;
}




