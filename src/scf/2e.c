#include <math.h>

#include "ints.h"

/***********************************************************************
 *                              GLOBALS
 * see 1e.c for definitions
 ***********************************************************************/
/*extern double alpha;
extern double beta;
extern double p, q;
extern double mu;
extern double Rpc2;
extern double Kab;

extern double Pat[3], Qat[3];
extern double *Aat, *Bat, *Cat, *Dat;

double Kab_xyz;
double Kcd_xyz;*/


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
	int set[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	
	Li = fi->f->L;
	Lj = fj->f->L;
	Lk = fk->f->L;
	Ll = fl->f->L;
	double (*Ni)(double) = N0, (*Nj)(double) = N0, (*Nk)(double) = N0, (*Nl)(double) = N0;
	if (Li == BFN_P) {
		Ni = N1;
		set[4*fi->m] = 1;
	}
	if (Lj == BFN_P) {
		Nj = N1;
		set[4*fj->m+1] = 1;
	}
	if (Lk == BFN_P) {
		Nk = N1;
		set[4*fk->m+2] = 1;
	}
	if (Ll == BFN_P) {
		Nl = N1;
		set[4*fl->m+3] = 1;
	}
	
	for (a = 0; a < inprim; a++)
		for (b = 0; b < jnprim; b++)
			for (c = 0; c < knprim; c++)
				for (d = 0; d < lnprim; d++) {
					double expa = fi->f->exp[a];
					double expb = fj->f->exp[b];
					double expc = fk->f->exp[c];
					double expd = fl->f->exp[d];
					double N = Ni(expa)*Nj(expb)*Nk(expc)*Nl(expd);
					double coeff = fi->f->c[a] * fj->f->c[b] * fk->f->c[c] * fl->f->c[d];
					eri += N * coeff * ERI(set, fi->a->r, fj->a->r, fk->a->r, fl->a->r,
						expa, expb, expc, expd);
				}
	return eri;
}





/*
double Theta(int N, int *ijkl)
{
	int i, j, k, l, t;
	double result = 0.0;
	double Xpa, Xpq, Xpb, Xqc, Xqd;
	
	if (ijkl[0] + ijkl[1] + ijkl[2] + ijkl[3] + ijkl[4] + ijkl[5] +
		ijkl[6] + ijkl[7] + ijkl[8] + ijkl[9] + ijkl[10] + ijkl[11] == 0) {
		return 2.0*pow(M_PI,2.5)/(p*q*sqrt(p+q)) * Kab_xyz * Kcd_xyz * boys(N, alpha*dist2(Pat, Qat));
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
		result += Xpa*Theta(N, ijkl);
		result += -alpha/p*Xpq*Theta(N+1, ijkl);
		if (i != 0) {
			ijkl[4*t] -= 1;
			result += 0.5*i/p*(Theta(N, ijkl) - alpha/p*Theta(N+1, ijkl));
			ijkl[4*t] += 1;
		}
		if (j != 0) {
			ijkl[4*t+1] -= 1;
			result += 0.5*j/p*(Theta(N, ijkl) - alpha/p*Theta(N+1, ijkl));
			ijkl[4*t+1] += 1;
		}
		if (k != 0) {
			ijkl[4*t+2] -= 1;
			result += 0.5*k/(p+q) * Theta(N+1, ijkl);
			ijkl[4*t+2] += 1;
		}
		if (l != 0) {
			ijkl[4*t+3] -= 1;
			result += 0.5*l/(p+q) * Theta(N+1, ijkl);
			ijkl[4*t+3] += 1;
		}
		ijkl[4*t] += 1;
	}
	else if (j != 0) {
		j -= 1;
		ijkl[4*t+1] -= 1;
		result += Xpb*Theta(N, ijkl);
		result += -alpha/p*Xpq*Theta(N+1, ijkl);
		if (i != 0) {
			ijkl[4*t] -= 1;
			result += 0.5*i/p*(Theta(N, ijkl) - alpha/p*Theta(N+1, ijkl));
			ijkl[4*t] += 1;
		}
		if (j != 0) {
			ijkl[4*t+1] -= 1;
			result += 0.5*j/p*(Theta(N, ijkl) - alpha/p*Theta(N+1, ijkl));
			ijkl[4*t+1] += 1;
		}
		if (k != 0) {
			ijkl[4*t+2] -= 1;
			result += 0.5*k/(p+q) * Theta(N+1, ijkl);
			ijkl[4*t+2] += 1;
		}
		if (l != 0) {
			ijkl[4*t+3] -= 1;
			result += 0.5*l/(p+q) * Theta(N+1, ijkl);
			ijkl[4*t+3] += 1;
		}
		ijkl[4*t+1] += 1;
	}
	else if (k != 0) {
		k -= 1;
		ijkl[4*t+2] -= 1;
		result += Xqc*Theta(N, ijkl);
		result += alpha/q*Xpq*Theta(N+1, ijkl);
		if (k != 0) {
			ijkl[4*t+2] -= 1;
			result += 0.5*k/q*(Theta(N, ijkl) - alpha/q*Theta(N+1, ijkl));
			ijkl[4*t+2] += 1;
		}
		if (l != 0) {
			ijkl[4*t+3] -= 1;
			result += 0.5*l/q*(Theta(N, ijkl) - alpha/q*Theta(N+1, ijkl));
			ijkl[4*t+3] += 1;
		}
		if (i != 0) {
			ijkl[4*t] -= 1;
			result += 0.5*i/(p+q) * Theta(N+1, ijkl);
			ijkl[4*t] += 1;
		}
		if (j != 0) {
			ijkl[4*t+1] -= 1;
			result += 0.5*j/(p+q) * Theta(N+1, ijkl);
			ijkl[4*t+1] += 1;
		}
		ijkl[4*t+2] += 1;
	}
	else {  // l != 0
		l -= 1;
		ijkl[4*t+3] -= 1;
		result += Xqd*Theta(N, ijkl);
		result += alpha/q*Xpq*Theta(N+1, ijkl);
		if (k != 0) {
			ijkl[4*t+2] -= 1;
			result += 0.5*k/q*(Theta(N, ijkl) - alpha/q*Theta(N+1, ijkl));
			ijkl[4*t+2] += 1;
		}
		if (l != 0) {
			ijkl[4*t+3] -= 1;
			result += 0.5*l/q*(Theta(N, ijkl) - alpha/q*Theta(N+1, ijkl));
			ijkl[4*t+3] += 1;
		}
		if (i != 0) {
			ijkl[4*t] -= 1;
			result += 0.5*i/(p+q) * Theta(N+1, ijkl);
			ijkl[4*t] += 1;
		}
		if (j != 0) {
			ijkl[4*t+1] -= 1;
			result += 0.5*j/(p+q) * Theta(N+1, ijkl);
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
	
	p = a + b;
	q = c + d;
	alpha = p*q/(p+q);
	Kab_xyz = exp(-a*b*dist2(A, B)/p);
	Kcd_xyz = exp(-c*d*dist2(C, D)/q);
	Aat = A;
	Bat = B;
	Cat = C;
	Dat = D;
	
	for (i = 0; i < 3; i++) {
		Pat[i] = (a*A[i]+b*B[i])/p;
		Qat[i] = (c*C[i]+d*D[i])/q;
	}
	
	return Theta(0, ijkl);
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
	int set[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	
	Li = fi->f->L;
	Lj = fj->f->L;
	Lk = fk->f->L;
	Ll = fl->f->L;
	double (*Ni)(double) = N0, (*Nj)(double) = N0, (*Nk)(double) = N0, (*Nl)(double) = N0;
	if (Li == BFN_P) {
		Ni = N1;
		set[4*fi->m] = 1;
	}
	if (Lj == BFN_P) {
		Nj = N1;
		set[4*fj->m+1] = 1;
	}
	if (Lk == BFN_P) {
		Nk = N1;
		set[4*fk->m+2] = 1;
	}
	if (Ll == BFN_P) {
		Nl = N1;
		set[4*fl->m+3] = 1;
	}
	
	for (a = 0; a < inprim; a++)
		for (b = 0; b < jnprim; b++)
			for (c = 0; c < knprim; c++)
				for (d = 0; d < lnprim; d++) {
					double expa = fi->f->exp[a];
					double expb = fj->f->exp[b];
					double expc = fk->f->exp[c];
					double expd = fl->f->exp[d];
					double N = Ni(expa)*Nj(expb)*Nk(expc)*Nl(expd);
					double coeff = fi->f->c[a] * fj->f->c[b] * fk->f->c[c] * fl->f->c[d];
					eri += N * coeff * ERI(set, fi->a->r, fj->a->r, fk->a->r, fl->a->r,
						expa, expb, expc, expd);
				}
	return eri;
}

*/
