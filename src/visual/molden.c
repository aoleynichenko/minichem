/***********************************************************************
 * molden.c
 * 
 * Interface to the MOLDEN visualization program:
 * http://www.cmbi.ru.nl/molden/
 * 
 * 2016-2018 Alexander Oleynichenko
 **********************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "visual.h"
#include "input.h"
#include "chem.h"
#include "basis.h"
#include "util.h"

void vectors_molden(struct cart_mol *mol, double *mo, double *en, int *occup, double n)
{
	int i, j, k;
	char filename[256];
	FILE *f = NULL;
	int geom_units;
	char name[64];
	
	rtdb_get("geom:units", &geom_units);
	rtdb_get("top:short_name", name);
	
	sprintf(filename, "%s.mos", name);
	f = fopen(filename, "w");
	if (!f)
		errquit("while writing vectors to molden format file: unable to create file");
	
	fprintf(f, "[Molden Format]\n");
	fprintf(f, "[Atoms] %s\n", geom_units == UNITS_ANGSTROMS ? "Angs" : "AU");
	
	for (i = 0; i < mol->size; i++) {
		int Z = mol->atoms[i].Z;
		double *r = mol->atoms[i].r;
		fprintf(f, "%3s%4d%4d%14.8f%14.8f%14.8f\n", searchByZ(Z)->sym, i+1, Z, r[0], r[1], r[2]);
	}
	fprintf(f, "[5D]\n");
	fprintf(f, "[GTO]\n");
	for (i = 0; i < mol->size; i++) {
		struct basis_set *bs = searchByZ(mol->atoms[i].Z)->bas;
		fprintf(f, "%4d\n", i+1);
		for (j = 0; j < bs->size; j++) {
			struct cgtf *fun = &bs->cgtfs[j];
			fprintf(f, " %s%3d\n", fun->L == BFN_S ? "s" : "p", fun->nprim);
			for (k = 0; k < fun->nprim; k++)
				fprintf(f, " %.8e  %.8e\n", fun->exp[k], fun->c[k]);
		}
		fprintf(f, "\n");
	}
	fprintf(f, "[MO]\n");
	for (i = 0; i < n; i++) {
		fprintf(f, "Sym=C1\n");
		fprintf(f, "Ene=%14.7f\n", 27.211396132*en[i]);
		fprintf(f, "Spin=Alpha\n");
		fprintf(f, "Occup=%.8f\n", (double) occup[i]);
		for (j = 0; j < n; j++) {
			int ind = i*n+j;
			fprintf(f, "%4d%13.8f\n", j+1, mo[ind]);
		}
	}
	fprintf(f, "\n");
	
	fclose(f);
}
