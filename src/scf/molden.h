#pragma once

#include "../input/chem.h"
#include "../input/basis.h"

void vectors_molden(struct cart_mol *mol, double *mo, double *en, int *occup, double n, char *filename);
