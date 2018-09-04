/***********************************************************************
 * eht.c
 * =====
 * 
 * Extended Huckel Method implementation -- Wolfsberg-Helmholtz variant.
 * Is used in minichem only for producing SCF initial guess.
 * 2016-2018 Alexander Oleynichenko.
 * 
 * Algorithm notes:
 * 
 * FC = ESC >>> C(0) >>> SCF start
 *   S -- real AO overlap matrix (no approximations)
 *   F -- "empirically" constructed matrix
 *        Fii = - IP[AO_i]
 *        Fij = 1/2*K*(Fii + Fjj)*Sij, K = 1.75
 * 
 * IP's are estimated using Koopmans' theorem at HF/cc-pVTZ level.
 **********************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "scf.h"
#include "chem.h"
#include "basis.h"
#include "util.h"

/* Ionization potentials
 * Two-dimensional array of format:
 * [Z][] = {1s 2s 2p 3s 3p 3d ...}
 */
double IP[][5] = {
/* H  */  { -0.50,  0.00,  0.00, 0.0, 0.0},
/* He */  { -0.91,  0.00,  0.00, 0.0, 0.0},
/* Li */  { -2.48, -0.20,  0.00, 0.0, 0.0},
/* Be */  { -4.73, -0.31,  0.00, 0.0, 0.0},
/* B  */  { -7.69, -0.50, -0.32, 0.0, 0.0},
/* C  */  {-11.32, -0.71, -0.44, 0.0, 0.0},
/* N  */  {-15.62, -0.94, -0.57, 0.0, 0.0},
/* O  */  {-20.67, -1.24, -0.67, 0.0, 0.0},
/* F  */  {-26.38, -1.57, -0.76, 0.0, 0.0},
/* Ne */  {-32.77, -1.93, -0.85, 0.0, 0.0}
};

/*    Fii = IP[i]
 *    Fij = 1/2*K*(Fii + Fjj)*Sij, K = 1.75
 */
void guess_F_eht(double *F, double *S, struct basis_function *bfn, int N)
{
	int i, j;
	
	for (i = 0; i < N; i++) {
		int n = bfn[i].f->n;
		int l = bfn[i].f->L;
		int z = bfn[i].a->Z;
		F[i*N+i] = IP[z-1][n-1+l];
	}
	for (i = 0; i < N; i++)
		for (j = i+1; j < N; j++) {
			F[i*N+j] = 0.875 * (F[i*N+i] + F[j*N+j]) * S[i*N+j];
			F[j*N+i] = F[i*N+j];
		}
}








