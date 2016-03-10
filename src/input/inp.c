/***********************************************************************
 * Kernel of minichem
 * 
 * 
 ***********************************************************************
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#include "basis.h"
#include "chem.h"
#include "input.h"
#include "lexer.h"
#include "../scf/scf.h"
#include "../util/util.h"

static int rank;

calc_information calc_info;

/* from lexer.c */
extern char *source; /* input file */
extern int fsize;

void directive_start();
void directive_memory();
void directive_echo();
void directive_geometry();
void directive_basis();
void directive_task();
void directive_charge();
void directive_out();
void calc_info_defaults();
void directive_nproc();

void compute(char *filename)
{
	struct cart_mol *mol;
	int i;
	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &calc_info.nproc);
	
	/* master (0) reads file and broadcasts data to slaves */
	if (rank == 0) {
		FILE *f = fopen(filename, "r");
		if (!f) {
			printf("Error: cannot open file '%s'\n", filename);
			errquit("an error occured while reading input");
		}
		load(f, filename);
	}
	
	MPI_Bcast(&fsize, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if (rank != 0) {
		if (source)
			free(source);
		source = (char *) malloc(sizeof(char) * fsize + 10);
	}
	MPI_Bcast(source, fsize, MPI_BYTE, 0, MPI_COMM_WORLD);
	if (rank != 0) {
		initLexer();
		source[fsize-1] = '\0';
	}
	
	/* some defaults which must be set before every new task */
	calc_info_defaults();
	scf_init();
	strcpy(calc_info.name, filename);
	
	nextToken();
	while (ttype != TT_EOF) {
		switch (ttype) {
		case TT_KW_START:
			directive_start();
			break;
		case TT_KW_MEMORY:
			directive_memory();
			break;
		case TT_KW_ECHO:
			directive_echo();
			break;
		case TT_KW_OUT:
			directive_out();
			break;
		case TT_KW_CHARGE:
			directive_charge();
			break;
		case TT_KW_NPROC:
			directive_nproc();
			break;
		case TT_KW_GEOMETRY:
			directive_geometry();
			setbuf(stdout, NULL);
			mol = &calc_info.molecule;
			/*printf("size = %d\n", mol->size);
			for (i = 0; i < mol->size; i++) {
				printf("%d %g %g %g\n", mol->atoms[i].Z, mol->atoms[i].r[0], mol->atoms[i].r[1], mol->atoms[i].r[2]);
			}*/
			break;
		case TT_KW_BASIS:
			directive_basis();
			if (rank == 0) {
				print_basis_summary();
				line_separator();
			}
			break;
		case TT_KW_SCF:
			directive_scf();
			break;
		case TT_KW_TASK:
			directive_task();
			break;
		}
		nextToken();
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
}

void directive_start()
{
	nextToken();
	if (ttype != TT_WORD)
		errquit("in input file: identifier expected");
	strcpy(calc_info.name, sval);
}

void directive_out()
{
	nextToken();
	for (;;) {
		if (ttype == TT_KW_END)
			break;
		else if (ttype == TT_EOF)
			errquit("reached unexpected end of input file in section 'out'");
		else if (ttype == TT_WORD) {
			if (!strcmp(sval, "molden"))
				calc_info.out_molden_vectors = 1;
			else
				errquit("unknown keyword in 'out' section");
		}
		nextToken();
	}
}

void directive_nproc()
{
	int nproc;
	char ompv[6];
	
	match(TT_NUMBER);
	nproc = (int) nval;
	calc_info.nproc = nproc;
	omp_set_num_threads(nproc);
	
	switch (_OPENMP) {
		case 200505: strcpy(ompv, "2.5"); break;
		case 200805: strcpy(ompv, "3.0"); break;
		case 201107: strcpy(ompv, "3.1"); break;
		case 201307: strcpy(ompv, "4.0"); break;
		case 201511: strcpy(ompv, "4.5"); break;
		default:   strcpy(ompv, "undef"); break;
	}
	
	printf("\n                *****************************\n");
	printf("                *           OPENMP          *\n");
	printf("                *           ------          *\n");
	printf("                * Number of threads      %2d *\n", omp_get_max_threads());
	printf("                * Cores available        %2d *\n", omp_get_num_procs());
	printf("                * OpenMP version      %5s *\n", ompv);
	printf("                *****************************\n\n");
}

void directive_memory()
{
	int memsize;
	
	nextToken();
	if (ttype != TT_NUMBER)
		errquit("memory directive: number expected");
	memsize = (int) nval;
	nextToken();
	if (ttype != TT_WORD)
		errquit("memory directive: one of keywords mb, mw, gb expected");
	if (!strcmp(sval, "b"))
		calc_info.memory = memsize;
	else if (!strcmp(sval, "kb"))
		calc_info.memory = memsize*1024;
	else if (!strcmp(sval, "mb"))
		calc_info.memory = memsize*1024*1024;
	else if (!strcmp(sval, "mw"))
		calc_info.memory = memsize*1024*1024/2;
	else if (!strcmp(sval, "gb"))
		calc_info.memory = memsize*1024*1024*1024;
	else
		errquit("memory directive: one of keywords b, kb, mb, mw, gb expected");
	setmemavail(calc_info.memory);
}

void directive_echo()
{
	calc_info.echo = 1;
}

void directive_charge()
{
	match(TT_NUMBER);
	calc_info.molecule.charge = nval;
}

void directive_geometry()
{
	int i;
	
	for (;;) {
		nextToken(); /* sval == NULL ---> seg fault */
		if (ttype == TT_WORD) {
			if (!strcmp(sval, "units")) {
				nextToken();
				if (ttype != TT_WORD)
					errquit("in geometry input: expected 'angstroms' or 'atomic' after 'units' keyword");
				if (!strcmp(sval, "atomic"))
					calc_info.geom_units = UNITS_ATOMIC;
				else if (!strcmp(sval, "angstroms"))
					calc_info.geom_units = UNITS_ANGSTROMS;
				else
					errquit("in geometry input: expected 'angstroms' or 'atomic' after 'units' keyword");
			}
			/*else if (!strcmp(sval, "charge")) {
				match(TT_NUMBER);
				calc_info.molecule.charge = nval;
			}*/
			else if (!strcmp(sval, "mult")) {
				match(TT_NUMBER);
				calc_info.molecule.mult = nval;
			}
			else {
				double x, y, z;
				struct elem_info *elem = searchBySym(sval);
				if (!elem)
					errquit("unknown element");
				match(TT_NUMBER);
				x = nval;
				match(TT_NUMBER);
				y = nval;
				match(TT_NUMBER);
				z = nval;
				
				if (calc_info.geom_units == UNITS_ANGSTROMS) {
					x *= 1.889725989;
					y *= 1.889725989;
					z *= 1.889725989;
				}
				append_atom(&calc_info.molecule, elem->Z, x, y, z);
			}
		}
		else if (ttype == TT_KW_END)
			break;
		else {
			char buf[128];
			sprintf("unknown directive in input file: %s", sval == NULL ? "" : sval);
			errquit(buf);
		}
	}
	/*if (rank == 0) {
		mol_summary(&calc_info.molecule);
		line_separator();
	}*/
}

/* run calculation */
void directive_task()
{
	nextToken();
	if (ttype != TT_KW_SCF)
		errquit("only SCF calculations can be performed");
	nextToken();
	if (!(ttype == TT_WORD && sval && !strcmp(sval, "energy"))) /* task scf energ */
		lexerPushBack();
	scf_energy(&calc_info.molecule);
}

void calc_info_defaults()
{
	calc_info.echo = 0;
	strcpy(calc_info.name, "John-A-Pople");
	calc_info.nproc = 1;
	calc_info.memory = 200*1024*1014;  /* 200 mb per thread */
	/* molecule */
	calc_info.geom_units        = UNITS_ANGSTROMS;
	calc_info.molecule.size     = 0;
	calc_info.molecule.capacity = 0;
	calc_info.molecule.charge   = 0;
	calc_info.molecule.mult     = 1; /* singlet */
	/* output */
	calc_info.out_molden_vectors = 0;
}







