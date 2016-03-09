#pragma once

#include "basis.h"

struct elem_info {
	int Z;
	char sym[3];
	double m; /* default atomic mass */
	struct basis_set *bas;
};

/*extern struct elem_info ptable[];*/

struct atom {
	int Z;
	double r[3];
};

struct cart_mol {
	int size;
	int capacity;
	int charge;
	int mult;
	struct atom *atoms;
};

void mol_summary(struct cart_mol *molecule);
struct elem_info *searchBySym(char *s);
struct elem_info *searchByZ(int z);
void append_atom(struct cart_mol *m, int Z, double x, double y, double z);
int nelec(struct cart_mol *mol);
