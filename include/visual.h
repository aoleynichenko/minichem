/***********************************************************************
 * visual.h
 * ========
 * 
 * 2016-2018 Alexander Oleynichenko
 **********************************************************************/

#ifndef VISUAL_H_INCLUDED
#define VISUAL_H_INCLUDED

#include "chem.h"
#include "basis.h"

void vectors_molden(struct cart_mol *mol, double *mo, double *en, int *occup, double n, char *filename);

#endif /* VISUAL_H_INCLUDED */
