/***********************************************************************
 * inp.c
 * =====
 * 
 * Kernel of minichem.
 * minichem works as a simple interpreter, it collects data from the
 * input file and performs calculations when it finds the keyword 'task'
 * 
 * 2016-2018 Alexander Oleynichenko
 * 
 **********************************************************************/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#include "basis.h"
#include "chem.h"
#include "input.h"
#include "lexer.h"
#include "scf.h"
#include "util.h"


// MPI rank
static int rank;

/* molecule under consideration -- geometry, ... etc */
struct cart_mol molecule;

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


/***********************************************************************
 * compute
 * 
 * This subroutine perfoms parsing and 'execution' of the input file
 * 'filename'.
 * 
 * TODO: refactoring -- maybe rename it?
 **********************************************************************/
void compute(char *filename)
{
	struct cart_mol *mol;
	int i;
	int one = 1;
	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &one /*&calc_info.nproc*/);
	
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
			mol = &molecule;
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


// start -- name of the job/task/calculation...
void directive_start()
{
	nextToken();
	if (ttype != TT_WORD)
		errquit("in input file: identifier expected");
	rtdb_set("top:short_name", "%s", sval);
}


/***********************************************************************
 * directive_out
 * 
 * Flushes data for another programs (xyz, MOs, density, etc...)
 * Supported output formats:
 *  - Molden (vectors -- MO visualization)
 **********************************************************************/
void directive_out()
{
	nextToken();
	for (;;) {
		if (ttype == TT_KW_END)
			break;
		else if (ttype == TT_EOF)
			errquit("reached unexpected end of input file in section 'out'");
		else if (ttype == TT_WORD) {
			if (!strcmp(sval, "molden")) {
				rtdb_set("visual:molden", "%i", 1);
			}
			else
				errquit("unknown keyword in 'out' section");
		}
		nextToken();
	}
}


/***********************************************************************
 * directive_nproc
 * 
 * Set number of openmp threads.
 **********************************************************************/
void directive_nproc()
{
	int nproc;
	char ompv[6];
	
	match(TT_NUMBER);
	nproc = (int) nval;
	rtdb_set("top:nproc", "%i", nproc);
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


/***********************************************************************
 * directive_memory
 * 
 * Set max allowed memory (RAM) usage.
 **********************************************************************/
void directive_memory()
{
	int memsize;
	int calc_memory;
	
	nextToken();
	if (ttype != TT_NUMBER)
		errquit("memory directive: number expected");
	memsize = (int) nval;
	nextToken();
	if (ttype != TT_WORD)
		errquit("memory directive: one of keywords mb, mw, gb expected");
	if (!strcmp(sval, "b"))
		calc_memory = memsize;
	else if (!strcmp(sval, "kb"))
		calc_memory = memsize*1024;
	else if (!strcmp(sval, "mb"))
		calc_memory = memsize*1024*1024;
	else if (!strcmp(sval, "mw"))
		calc_memory = memsize*1024*1024/2;
	else if (!strcmp(sval, "gb"))
		calc_memory = memsize*1024*1024*1024;
	else
		errquit("memory directive: one of keywords b, kb, mb, mw, gb expected");
	
	rtdb_set("top:memory", "%i", calc_memory);
	setmemavail(calc_memory);
}


/***********************************************************************
 * directive_echo
 * 
 * Echo input file or not?
 **********************************************************************/
void directive_echo()
{
	rtdb_set("top:echo", "%i", 1);
}


/***********************************************************************
 * directive_charge
 * 
 * Set molecular (total) charge.
 * 
 * TODO: just for compatibility of input files with nwchem. To be
 * removed or it is a useful feature?
 **********************************************************************/
void directive_charge()
{
	match(TT_NUMBER);
	rtdb_set("geom:charge", "%i", nval);
}


/***********************************************************************
 * directive_geometry
 * 
 * Set molecular geometry (xyz).
 **********************************************************************/
void directive_geometry()
{
	int i;
	int geom_units = UNITS_ANGSTROMS;
	
	for (;;) {
		nextToken(); /* sval == NULL ---> seg fault */
		if (ttype == TT_WORD) {
			if (!strcmp(sval, "units")) {
				nextToken();
				if (ttype != TT_WORD)
					errquit("in geometry input: expected 'angstroms' or 'atomic' after 'units' keyword");
				if (!strcmp(sval, "atomic")) {
					geom_units = UNITS_ATOMIC;
					rtdb_set("geom:units", "%i", UNITS_ATOMIC);
				}
				else if (!strcmp(sval, "angstroms")) {
					geom_units = UNITS_ANGSTROMS;
					rtdb_set("geom:units", "%i", UNITS_ANGSTROMS);
				}
				else
					errquit("in geometry input: expected 'angstroms' or 'atomic' after 'units' keyword");
			}
			// TODO: enable this code?
			/*else if (!strcmp(sval, "charge")) {
				match(TT_NUMBER);
				calc_info.molecule.charge = nval;
			}*/
			else if (!strcmp(sval, "mult")) {
				match(TT_NUMBER);
				molecule.mult = nval;
				rtdb_set("geom:mult", "%i", nval);
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
				
				if (geom_units == UNITS_ANGSTROMS) {
					x *= 1.889725989;
					y *= 1.889725989;
					z *= 1.889725989;
				}
				append_atom(&molecule, elem->Z, x, y, z);
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
	// TODO: enable this code?
	/*if (rank == 0) {
		mol_summary(&calc_info.molecule);
		line_separator();
	}*/
}


/***********************************************************************
 * directive_task
 * 
 * Starts quantum chemistry calculation.
 **********************************************************************/
void directive_task()
{
	nextToken();
	if (ttype != TT_KW_SCF)
		errquit("only SCF calculations can be performed");
	nextToken();
	if (!(ttype == TT_WORD && sval && !strcmp(sval, "energy"))) /* task scf energ */
		lexerPushBack();
	scf_energy(&molecule);
}


/***********************************************************************
 * calc_info_defaults
 * 
 * Default settings for the minichem.
 **********************************************************************/
void calc_info_defaults()
{
	/* molecule */
	molecule.size     = 0;
	molecule.capacity = 0;
	molecule.charge   = 0;
	molecule.mult     = 1; /* singlet */
	
	// and the same to the RTDB
	/* top-level */
	rtdb_set("top:echo", "%i", 0);
	rtdb_set("top:short_name", "%s", "John-A-Pople");
	rtdb_set("top:nproc", "%i", 1);
	rtdb_set("top:memory", "%i", 200*1024*1024);
	/* molecule */
	rtdb_set("geom:units",  "%i", UNITS_ANGSTROMS);
	rtdb_set("geom:charge", "%i", 0);  /* neutral */
	rtdb_set("geom:mult",   "%i", 1);  /* singlet */
	/* output */
	rtdb_set("visual:molden", "%i", 0);
}







