/***********************************************************************
 * input.h
 * =======
 * 
 * 2016-2018 Alexander Oleynichenko
 **********************************************************************/

#ifndef INPUT_H_INCLUDED
#define INPUT_H_INCLUDED

#include "chem.h"

void calc_info_defaults();
void compute(char *filename);

#define UNITS_ANGSTROMS 0
#define UNITS_ATOMIC    1


typedef struct {
	int nproc;
	char name[128];  /* name of calculation, e.g. HeH+ */
	int memory;
	int echo;
	
	int geom_units;
	struct cart_mol molecule;
	
	// output options
	int out_molden_vectors;
} calc_information;


extern calc_information calc_info;

/* interface to the runtime database (rtdb) */
int rtdb_set(char *key, char *fmt, ...);
int rtdb_get(char *key, ...);
void rtdb_print_meta();

#endif /* INPUT_H_INCLUDED */
