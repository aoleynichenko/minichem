#pragma once

#include "chem.h"

void calc_info_defaults();
void compute(char *filename);
void interactive_mode();

#define UNITS_ANGSTROMS 0
#define UNITS_ATOMIC    1

typedef struct {
	int nalpha;
	int nbeta;
	int type;
	void *data;
} WaveFunction_t;

typedef struct {
	int nproc;
	char name[128];  /* name of calculation, e.g. HeH+ */
	int memory;
	int echo;
	
	int geom_units;
	struct cart_mol molecule;
	WaveFunction_t wf;
	
	// output options
	int out_molden_vectors;
} calc_information;

#define WF_NOTHING 0
#define WF_RHF     1
#define WF_UHF     2

extern calc_information calc_info;


