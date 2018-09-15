/***********************************************************************
 * chem.h
 * ======
 * 
 * 2016-2018 Alexander Oleynichenko
 **********************************************************************/

#ifndef CHEM_H_INCLUDED
#define CHEM_H_INCLUDED

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
typedef struct atom Atom_t;

struct cart_mol {
	int size;
	int capacity;
	int charge;
	int mult;
	struct atom *atoms;
};
typedef struct cart_mol       Molecule_t;

void mol_summary(struct cart_mol *molecule);
struct elem_info *searchBySym(char *s);
struct elem_info *searchByZ(int z);
void append_atom(struct cart_mol *m, int Z, double x, double y, double z);
void print_molecule(Molecule_t *geom, int units);
void distance_matrix(Molecule_t *geom, int units);
int atoms_are_equal(Atom_t *a, Atom_t *b);
int nelec(struct cart_mol *mol);
int nalphabeta(struct cart_mol *mol, int *Nalpha, int *Nbeta);

#endif /* CHEM_H_INCLUDED */

