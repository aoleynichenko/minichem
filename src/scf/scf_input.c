#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "lexer.h"
#include "util.h"
#include "input.h"
#include "scf.h"

void scf_print_opts(int doprint);
void scf_guess_inp();

void directive_scf()
{
	nextToken();
	for (;;) {
		if (ttype == TT_KW_END)
			break;
		else if (ttype == TT_EOF)
			errquit("reached unexpected end of input file in section 'scf'");
		else if (ttype == TT_WORD) {
			if (!strcmp(sval, "print"))
				scf_print_opts(1);
			else if (!strcmp(sval, "guess"))
				scf_guess_inp();
			else if (!strcmp(sval, "noprint"))
				scf_print_opts(0);
			else if (!strcmp(sval, "rhf")) {
				rtdb_set("scf:wf", "%i", SCF_RHF);
			}
			else if (!strcmp(sval, "uhf")) {
				rtdb_set("scf:wf", "%i", SCF_UHF);
			}
			// set multiplicity, NWChem compatible notation
			else if (!strcmp(sval, "singlet")) {
				rtdb_set("geom:mult", "%i", 1);
				molecule.mult = 1;
			}
			else if (!strcmp(sval, "doublet")) {
				rtdb_set("geom:mult", "%i", 2);
				molecule.mult = 2;
			}
			else if (!strcmp(sval, "triplet")) {
				rtdb_set("geom:mult", "%i", 3);
				molecule.mult = 3;
			}
			else if (!strcmp(sval, "quartet")) {
				rtdb_set("geom:mult", "%i", 4);
				molecule.mult = 4;
			}
			else if (!strcmp(sval, "quintet")) {
				rtdb_set("geom:mult", "%i", 5);
				molecule.mult = 5;
			}
			// convergence options
			else if (!strcmp(sval, "maxiter")) {
				match(TT_NUMBER);
				if (fabs(nval - (int)nval) < 1e-10) {
					rtdb_set("scf:maxiter", "%i", (int) nval);
				}
				else
					errquit("in scf:maxiter: integer value expected in input!");
			}
			else if (!strcmp(sval, "diis")) {
				nextToken();
				if (ttype == TT_NUMBER) {
					if (fabs(nval - (int)nval) < 1e-10) {
						rtdb_set("scf:diisbas", "%i", (int) nval);
					}
					else
						errquit("in scf:diis: integer value expected in input!");
				}
				else
					lexerPushBack();
				rtdb_set("scf:diis", "%i", 1);
			}
			else if (!strcmp(sval, "nodiis")) {
				rtdb_set("scf:diis", "%i", 0);
			}
			// conventional/direct
			else if (!strcmp(sval, "direct")) {
				rtdb_set("scf:direct", "%i", 1);
			}
            else if (!strcmp(sval, "nodirect")) {
				rtdb_set("scf:direct", "%i", 0);
			}
			else
				errquit("unknown keyword in scf section");
		}
		
		nextToken();
	}
}

void scf_print_opts(int doprint)
{
	nextToken();
	if (ttype != TT_WORD && ttype != '"')
		errquit("unexpected input in scf:[no]print");
	
	if (!strcmp(sval, "overlap")) {
		rtdb_set("aoints:print:overlap", "%i", doprint);
	}
	else if (!strcmp(sval, "kinetic")) {
		rtdb_set("aoints:print:kinetic", "%i", doprint);
	}
	else if (!strcmp(sval, "potential")) {
		rtdb_set("aoints:print:potential", "%i", doprint);
	}
	else if (!strcmp(sval, "final vectors analysis")) {
		rtdb_set("scf:print:vectors", "%i", doprint);
	}
	else if (!strcmp(sval, "eri") || !strcmp(sval, "ao2eints")) {
		rtdb_set("aoints:print:eri", "%i", doprint);
	}
}

void scf_guess_inp()
{
	nextToken();
	if (ttype != TT_WORD)
		errquit("unexpected input in scf:guess");
	
	if (!strcmp(sval, "core")) {
		rtdb_set("scf:guess", "%i", GUESS_BARE);
	}
	else if (!strcmp(sval, "eht")) {
		rtdb_set("scf:guess", "%i", GUESS_EHT);
	}
	else
		errquit("wrong input in scf:guess");
}




