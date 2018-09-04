/***********************************************************************
 * chem.c
 * ======
 * 
 * Periodic table of elements and related routines.
 * 
 * 2016-2018 Alexander Oleynichenko
 * 
 **********************************************************************/

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "chem.h"
#include "input.h"
#include "util.h"

#define NELEMENTS 120

struct elem_info ptable[] = {
    {1,  "H",   1.00795,  NULL },
    {2,  "He",  4.002602, NULL },
    {3,  "Li",  6.9412,   NULL },
    {4,  "Be",  9.01218,  NULL },
    {5,  "B",   10.812,   NULL },
    {6,  "C",   12.0108,  NULL },
    {7,  "N",   14.0067,  NULL },
    {8,  "O",   15.9994,  NULL },
    {9,  "F",   18.9984,  NULL },
    {10, "Ne",  20.179,   NULL },
    {11,  "Na", 22.9897,  NULL },
    {12,  "Mg", 24.305 ,  NULL },
    {13,  "Al", 26.9815,  NULL },
    {14,  "Si", 28.0855,  NULL },
    {15,  "P" , 30.9738,  NULL },
    {16,  "S" , 32.065 ,  NULL },
    {17,  "Cl", 35.453 ,  NULL },
    {18,  "Ar", 39.0983,  NULL },
    {19,  "K" , 39.948 ,  NULL },
    {20,  "Ca", 40.078 ,  NULL },
    {21,  "Sc", 44.9559,  NULL },
    {22,  "Ti", 47.867 ,  NULL },
    {23,  "V" , 50.9415,  NULL },
    {24,  "Cr", 51.9961,  NULL },
    {25,  "Mn", 54.938 ,  NULL },
    {26,  "Fe", 55.845 ,  NULL },
    {27,  "Co", 58.6934,  NULL },
    {28,  "Ni", 58.9332,  NULL },
    {29,  "Cu", 63.546  , NULL },
    {30,  "Zn", 65.39   , NULL },
    {31,  "Ga", 69.723  , NULL },
    {32,  "Ge", 72.64   , NULL },
    {33,  "As", 74.9216 , NULL },
    {34,  "Se", 78.96   , NULL },
    {35,  "Br", 79.904  , NULL },
    {36,  "Kr", 83.8    , NULL },
    {37,  "Rb", 85.4678 , NULL },
    {38,  "Sr", 87.62   , NULL },
    {39,  "Y" , 88.9059 , NULL },
    {40,  "Zr", 91.224  , NULL },
    {41,  "Nb", 92.9064 , NULL },
    {42,  "Mo", 95.94   , NULL },
    {43,  "Tc", 98      , NULL },
    {44,  "Ru", 101.07  , NULL },
    {45,  "Rh", 102.9055, NULL },
    {46,  "Pd", 106.42  , NULL },
    {47,  "Ag", 107.8682, NULL },
    {48,  "Cd", 112.411 , NULL },
    {49,  "In", 114.818 , NULL },
    {50,  "Sn", 118.71  , NULL },
    {51,  "Sb", 121.76  , NULL },
    {52,  "Te", 126.9045, NULL },
    {53,  "I" , 127.6   , NULL },
    {54,  "Xe", 131.293 , NULL },
    {55,  "Cs", 132.9055, NULL },
    {56,  "Ba", 137.327 , NULL },
    {57,  "La", 138.9055, NULL },
    {58,  "Ce", 140.116 , NULL },
    {59,  "Pr", 140.9077, NULL },
    {60,  "Nd", 144.24  , NULL },
    {61,  "Pm", 145     , NULL },
    {62,  "Sm", 150.36  , NULL },
    {63,  "Eu", 151.964 , NULL },
    {64,  "Gd", 157.25  , NULL },
    {65,  "Tb", 158.9253, NULL },
    {66,  "Dy", 162.5   , NULL },
    {67,  "Ho", 164.9303, NULL },
    {68,  "Er", 167.259 , NULL },
    {69,  "Tm", 168.9342, NULL },
    {70,  "Yb", 173.04  , NULL },
    {71,  "Lu", 174.967 , NULL },
    {72,  "Hf", 178.49  , NULL },
    {73,  "Ta", 180.9479, NULL },
    {74,  "W" , 183.84  , NULL },
    {75,  "Re", 186.207 , NULL },
    {76,  "Os", 190.23  , NULL },
    {77,  "Ir", 192.217 , NULL },
    {78,  "Pt", 195.078 , NULL },
    {79,  "Au", 196.9665, NULL },
    {80,  "Hg", 200.59,   NULL },
    {81,  "Tl", 204.3833, NULL },
    {82,  "Pb", 207.2,    NULL },
    {83,  "Bi", 208.9804, NULL },
    {84,  "Po", 209,      NULL },
    {85,  "At", 210,      NULL },
    {86,  "Rn", 222,      NULL },
    {87,  "Fr", 223,      NULL },
    {88,  "Ra", 226,      NULL },
    {89,  "Ac", 227,      NULL },
    {90,  "Th", 232.0381, NULL },
    {91,  "Pa", 231.0359, NULL },
    {92,  "U" , 238.0289, NULL },
    {93,  "Np", 237,      NULL },
    {94,  "Pu", 244,      NULL },
    {95,  "Am", 243,      NULL },
    {96,  "Cm", 247,      NULL },
    {97,  "Bk", 247,      NULL },
    {98,  "Cf", 251,      NULL },
    {99,  "Es", 252,      NULL },
    {100, "Fm", 257,      NULL },
    {101, "Md", 258,      NULL },
    {102, "No", 259,      NULL },
    {103, "Lr", 262,      NULL },
    {104, "Rf", 261,      NULL },
    {105, "Db", 262,      NULL },
    {106, "Sg", 266,      NULL },
    {107, "Bh", 264,      NULL },
    {108, "Hs", 277,      NULL },
    {109, "Mt", 268,      NULL },
    {110, "Ds", 0,        NULL },
    {111, "Rg", 272,      NULL },
    {112, "Cn", 0,        NULL },
    {113, "Nh", 0,        NULL },
    {114, "Fl", 0,        NULL },
    {115, "Mc", 0,        NULL },
    {116, "Lv", 0,        NULL },
    {117, "Ts", 0,        NULL },
    {118, "Og", 0,        NULL },
    {0,  NULL,  0.0,      NULL }
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
	
	for (; p->sym; p++)
		if (!strcicmp(p->sym, s))
			return p;
	return NULL;
}


struct elem_info *searchByZ(int z)
{
	struct elem_info *p = ptable;

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


/***********************************************************************
 * atoms_are_equal
 * 
 * compares two atoms, returns 1 if equal, 0 otherwise
 **********************************************************************/
int atoms_are_equal(Atom_t *a, Atom_t *b)
{
	static const double THR = 1e-14;
	
	return (a->Z == b->Z &&
			fabs(a->r[0] - b->r[0]) < THR &&
			fabs(a->r[1] - b->r[1]) < THR &&
			fabs(a->r[2] - b->r[2]) < THR);
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


/***********************************************************************
 * nalphabeta
 * 
 * Calculate number of alpha and beta electrons (Nalpha >= Nbeta)
 **********************************************************************/
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



