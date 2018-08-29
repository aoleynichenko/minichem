/***********************************************************************
 * basis.c
 * =======
 * 
 * Parses the 'basis' section of an input file.
 * Sample:
 * | basis "H_STO-3G" SPHERICAL
 * | H    S
 * |       3.42525091             0.15432897       
 * |       0.62391373             0.53532814       
 * |       0.16885540             0.44463454       
 * | end
 **********************************************************************/

#include "basis.h"
#include "chem.h"
#include "lexer.h"
#include "util.h"

#include <stdlib.h>
#include <string.h>

/***********************************************************************
 * 
 * global variables
 * 
 **********************************************************************/

/* task size -- number of basis functions */
int nbfns;

/* number of shells */
int nshells;

/* explicit set of basis functions */
/* all basis functions are centered on atom at (x,y,z) */
struct basis_function *bfns;
struct shell *shells;


// TODO: remove ths function
void determine_principal_numbers(struct basis_set *);


/***********************************************************************
 * directive_basis
 * 
 * Parses the 'basis' section of an input file.
 **********************************************************************/
void directive_basis()
{
	extern struct elem_info ptable[];   // see ./chem.c
	struct elem_info *p = &ptable[0];
	char basis_name[30];
	int basis_type = CARTESIAN;
	
	//printf("DIRECTIVE BASIS\n");
	nextToken();
	if (ttype == '"')
		strcpy(basis_name, sval);
	else {
		strcpy(basis_name, "user defined");
		lexerPushBack();
	}
	nextToken();
	if (ttype == TT_WORD && sval != NULL) {
		if (!strcmp(sval, "spherical"))
			basis_type = SPHERICAL;
		else if (!strcmp(sval, "cartesian"))
			basis_type = CARTESIAN;
		else
			lexerPushBack();
	}
	//printf("BASIS = %s LOOP\n", basis_type == SPHERICAL ? "SPHERICAL" : "CARTESIAN");
	
	for (;;) {
		nextToken(); /* sval == NULL ---> seg fault */
		//if (sval == NULL)
			//printf("NULL SVAL\n");
		if (ttype == TT_KW_END) {
			//printf("END OF BASIS REACHED\n");
			break;
		}
		else if (ttype != TT_WORD)
			errquit("in basis input: end keyword or element symbol is required");
		else {
			struct elem_info *elem = searchBySym(sval);
			int L, i, j, n, nprim, ncontr;
			double buf[10][10];
			
			if (!elem)
				errquit("in basis input: unknown element");
			//printf("ELEMENT Z = %d\n", elem->Z);
			match(TT_WORD);
			if (!strcmp(sval, "s"))
				L = BFN_S;
			else if (!strcmp(sval, "sp"))
				L = BFN_SP;
			else if (!strcmp(sval, "p"))
				L = BFN_P;
			else if (!strcmp(sval, "d"))
				L = BFN_D;
			else
				errquit("in basis input: wrong angular momentum (only S, SP, P, D are allowed)");
				
			if (!elem->bas) {   /* bind new basis set to element */
				elem->bas = (struct basis_set *) malloc (sizeof(struct basis_set));
				elem->bas->size = 0;
				elem->bas->type = basis_type;
				strcpy(elem->bas->name, basis_name);
			}
			
			nextToken();
			setEOLSignificant(1);
			i = 0;
			j = 0;
			ncontr = 0;
			while (ttype == TT_NUMBER) {
				buf[i][j++] = nval;
				nextToken();
				if (ttype == TT_EOL) {  // begin new row of doubles
					nextToken();
					i++;
					ncontr = j - 1;
					j = 0;
				}
			}
			nprim = i;
			setEOLSignificant(0);
			lexerPushBack();
			
			for (i = 0; i < ncontr; i++) {
				n = elem->bas->size++;
				if (i == 0 && L == BFN_SP) {
					elem->bas->cgtfs[n].L = BFN_S;
					L = BFN_P;
				}
				else
					elem->bas->cgtfs[n].L = L;
				elem->bas->cgtfs[n].nprim = nprim;
				for (j = 0; j < nprim; j++) {
					elem->bas->cgtfs[n].exp[j] = buf[j][0];
					elem->bas->cgtfs[n].c[j]   = buf[j][i+1];
				}
			}
			
			/*n = elem->bas->size;
			// config new contracted function
			elem->bas->cgtfs[n].nprim = 0;
			elem->bas->cgtfs[n].L = L;
			nprim = 0;
			// get exp & coeff
			nextToken();
			while (ttype == TT_NUMBER) {
				double exp = nval;
				if (exp <= 0)
					errquit("in basis input: exponent should be positive");
				elem->bas->cgtfs[n].exp[nprim] = exp;
				// get contraction coefficients
				match(TT_NUMBER);
				elem->bas->cgtfs[n].c[nprim] = nval;
				nprim++;
				elem->bas->cgtfs[n].nprim++; // to the next primitive
				nextToken();
			}
			elem->bas->size++;
			lexerPushBack();*/
		}
	}
	//printf("END OF DIRECTIVE BASIS\n");
	
	for (; *p->sym; p++) {
		if (!p->bas)
			continue;
		determine_principal_numbers(p->bas);
	}
}


/* Creates atom-centered basis functions from molecular data.
 * Returns: vector of atom-centered functions bfns with length M.
 * This function should be executed by all processes.
 * */
void form_atom_centered_bfns(struct cart_mol *molecule, struct basis_function **bfns, struct shell **shs, int *M, int *nsh)
{
	int i, j;
	int K = 0;
	int shn = 0;
	struct basis_function *p;
	struct shell *s;

	for (i = 0; i < molecule->size; i++) {
		struct elem_info *e = searchByZ(molecule->atoms[i].Z);
		if (e) {
			for (j = 0; j < e->bas->size; j++) {
				K += 2*e->bas->cgtfs[j].L + 1;  // пока и так сойдет
				shn++;
			}
		}
	}
	p = (struct basis_function *) malloc(K * sizeof(struct basis_function));
	s = (struct shell *) malloc(shn * sizeof(struct shell));
	*bfns = p;
	*shs = s;
	*M = K;
	*nsh = shn;
	for (i = 0; i < molecule->size; i++) {
		struct elem_info *e = searchByZ(molecule->atoms[i].Z);
		if (!e)
			continue;
		for (j = 0; j < e->bas->size; j++) {  // add atom-centered bfn
			int L = e->bas->cgtfs[j].L;
			int m;
			
			(*s).size = 2*L+1;
			(*s).start = p;
			for (m = 0; m < 2*L+1; m++) {
				(*p).m = m;
				(*p).f = &e->bas->cgtfs[j];
				(*p).a = &molecule->atoms[i];
				p++;
			}
			s++;
		}
	}
}


// TODO: delete this function (together with the implementation of the EHT)
// estimates principal quantum number n for every GTO in set 'bs'
// работает кое-как, но пока и ладно, сойдет
void determine_principal_numbers(struct basis_set *bs)
{
	int i, j;
	double maxs = 0.0, maxp = 0.0;
	
	for (i = 0; i < bs->size; i++) {
		struct cgtf *orb = &bs->cgtfs[i];
		int L = orb->L;
		double maxalpha = 0.0;
		for (j = 0; j < orb->nprim; j++)
			if (orb->exp[j] > maxalpha)
				maxalpha = orb->exp[j];
		if (L == BFN_S && maxalpha > maxs)
			maxs = maxalpha;
		else if (L == BFN_P && maxalpha > maxs)
			maxp = maxalpha;
	}
	
	// мы нашли максимальные альфы для каждого типа орбиталей
	// теперь попробуем оценить главные квантовые числа, считая, что для
	// каждого следующего n альфы отличаются на 1-2 порядка
	for (i = 0; i < bs->size; i++) {
		struct cgtf *orb = &bs->cgtfs[i];
		int L = orb->L;
		double maxalpha = 0.0;
		double ratio;
		
		for (j = 0; j < orb->nprim; j++)
			if (orb->exp[j] > maxalpha)
				maxalpha = orb->exp[j];
		if (L == BFN_S)
			ratio = maxs / maxalpha;
		else if (L == BFN_P)
			ratio = maxp / maxalpha;
		
		if (ratio < 2)
			orb->n = L + 1;
		else if (ratio < 70)
			orb->n = L + 1 + 1;
		else  // > 70
			orb->n = L + 1 + 2;
	}
}


/***********************************************************************
 * print_basis_summary
 * 
 * Print information about the current basis set
 * (exp-coef tables and summary).
 **********************************************************************/
void print_basis_summary()
{
	int i, j, type = -1;
	extern struct elem_info ptable[];
	struct elem_info *p = &ptable[0];
	static char *angmom[] = {"S", "P", "D", "F", "G", "H", "I"};
	
	printf("                    Basis information\n");
	for (; *p->sym; p++) {
		if (!p->bas)
			continue;
		type = p->bas->type;
		printf("   %s\n  ----\n", p->sym);
		for (i = 0; i < p->bas->size; i++) {
			struct cgtf *cf = &p->bas->cgtfs[i];
			printf("  %d   %1d%s%15.8f%15.8f\n", i+1, cf->n, angmom[cf->L], cf->exp[0], cf->c[0]);
			for (j = 1; j < cf->nprim; j++)
				printf("        %15.8f%15.8f\n", cf->exp[j], cf->c[j]);
			printf("\n");
		}
		printf("\n");
	}
	
	printf("%s basis set summary:\n", type == CARTESIAN ? "Cartesian" : "Spherical");
	printf("*---------*-----------------*--------*------------*\n");
	printf("| Element |   Description   | Shells | Primitives |\n");
	printf("*---------*-----------------*--------*------------*\n");
	for (p = &ptable[0]; *p->sym; p++) {
		int i, nprim = 0;
		if (!p->bas)
			continue;
			
		for (i = 0; i < p->bas->size; i++)
			nprim += p->bas->cgtfs[i].nprim;
		printf("|   %-2s    | %-15s |   %-2d   |     %-3d    |\n", p->sym, p->bas->name, p->bas->size, nprim);
	}
	printf("*---------*-----------------*--------*------------*\n\n");
}














