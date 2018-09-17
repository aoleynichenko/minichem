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

// molecule under consideration
extern struct cart_mol molecule;

// print level (for the top-level 'print' directive)
enum {
	PRINT_NONE,
	PRINT_LOW,
	PRINT_MEDIUM,
	PRINT_HIGH,
	PRINT_DEBUG
};
extern int print_level;  // default: PRINT_MEDIUM

/* interface to the runtime database (rtdb) */
int rtdb_set(char *key, char *fmt, ...);
int rtdb_get(char *key, ...);
void rtdb_print_meta();

#endif /* INPUT_H_INCLUDED */
