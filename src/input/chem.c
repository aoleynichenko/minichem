#include "chem.h"
#include "input.h"
#include "util.h"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NELEMENTS 100

struct elem_info ptable[] = {
	{1,  "H",   1.00795,  NULL},
	{2,  "He",  4.002602, NULL},
	{3,  "Li",  6.9412,   NULL},
	{4,  "Be",  9.01218,  NULL},
	{5,  "B",   10.812,   NULL},
	{6,  "C",   12.0108,  NULL},
	{7,  "N",   14.0067,  NULL},
	{8,  "O",   15.9994,  NULL},
	{9,  "F",   18.9984,  NULL},
	{10, "Ne",  20.179,   NULL},
	{0,  NULL,  0.0,      NULL}
};

/* Case-insensitive strcmp */
int strcicmp(char const *a, char const *b)
{
    for (;; a++, b++) {
        int d = tolower(*a) - tolower(*b);
        if (d != 0 || !*a)
            return d;
    }
}

struct elem_info *searchBySym(char *s)
{
	struct elem_info *p = ptable;
	
	//printf("SEARCH BY SYM '%s'\n", s);
	for (; p->sym; p++)
		if (!strcicmp(p->sym, s))
			return p;
	return NULL;
}

struct elem_info *searchByZ(int z)
{
	struct elem_info *p = ptable;
	//printf("SEARCH BY Z '%d'\n", z);
	for (; p->sym; p++)
		if (p->Z == z)
			return p;
	return NULL;
}

void append_atom(struct cart_mol *m, int Z, double x, double y, double z)
{
	struct atom *newat;
	int i;
	
	if (m->capacity == 0) {
		m->capacity = 8;
		m->atoms = (struct atom *) malloc (8 * sizeof(struct atom));
	}
	if (m->size == m->capacity) { /* realloc */
		newat = (struct atom *) malloc (sizeof(struct atom) * m->size*2);
		m->capacity *= 2;
		for (i = 0; i < m->size; i++)
			newat[i] = m->atoms[i];
		free(m->atoms);
		m->atoms = newat;
	}
	m->atoms[m->size].Z = Z;
	m->atoms[m->size].r[0] = x;
	m->atoms[m->size].r[1] = y;
	m->atoms[m->size].r[2] = z;
	//printf(" +  %d %.8f %.8f %.8f\n", Z, m->atoms[m->size].r[0], m->atoms[m->size].r[1], m->atoms[m->size].r[2]);
	m->size++;
}

void mol_summary(struct cart_mol *molecule)
{
	int el[NELEMENTS] = {0};
	int i;
	int nelec = -molecule->charge;
	double mass = 0.0;
	
	for (i = 0; i < molecule->size; i++)
		el[molecule->atoms[i].Z-1]++;
	printf("Formula:   ");
	for (i = 0; i < NELEMENTS; i++)
		if (el[i] != 0) {
			struct elem_info *ei = &ptable[i];
			mass += ei->m*el[i];
			nelec += ei->Z*el[i];
		}  // ???
	printf("\n");
	printf("Atoms:     %d\n", molecule->size);
	printf("Electrons: %d\n", nelec);
	printf("Mol mass:  %.3f\n", mass);
	printf("Charge:    %d\n", molecule->charge);
	printf("Spin mult: %d\n", molecule->mult);
	printf("Units:     %s\n", calc_info.geom_units == UNITS_ANGSTROMS ? "angstroms" : "atomic");
	/* verify charge & mult */
	//if (!((nelec % 2 == 0) && (molecule->mult % 2 == 1)))
		//errquit("illegal charge/multiplicity");
}

int nelec(struct cart_mol *mol)
{
	int i, N = 0;
	
	for (i = 0; i < mol->size; i++)
		N += mol->atoms[i].Z;
	N -= mol->charge;
	return N;
}

int nalphabeta(struct cart_mol *mol, int *Nalpha, int *Nbeta)
{
	int N = nelec(mol);
	int unpaired = mol->mult - 1;
	*Nbeta  = (N-unpaired)/2;
	*Nalpha = *Nbeta + unpaired;
	if (*Nalpha + *Nbeta != N) {
		printf("Nelec  = %d\n", N);
		printf("Nalpha = %d\n", *Nalpha);
		printf("Nbeta  = %d\n", *Nbeta);
		printf("Charge = %d\n", mol->charge);
		printf("Mult   = %d\n", mol->mult);
		errquit("uhf: illegal charge/multiplicity");
	}
	return N;
}



