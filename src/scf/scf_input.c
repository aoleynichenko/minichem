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
			else if (!strcmp(sval, "rhf"))
				scf_options.wavefuntype = SCF_RHF;
			else if (!strcmp(sval, "uhf"))
				scf_options.wavefuntype = SCF_UHF;
			else if (!strcmp(sval, "rohf"))
				scf_options.wavefuntype = SCF_ROHF;
			// set multiplicity, NWChem compatible notation
			else if (!strcmp(sval, "singlet"))
				calc_info.molecule.mult = 1;
			else if (!strcmp(sval, "doublet"))
				calc_info.molecule.mult = 2;
			else if (!strcmp(sval, "triplet"))
				calc_info.molecule.mult = 3;
			else if (!strcmp(sval, "quartet"))
				calc_info.molecule.mult = 4;
			else if (!strcmp(sval, "quintet"))
				calc_info.molecule.mult = 5;
			// convergence options
			else if (!strcmp(sval, "maxiter")) {
				match(TT_NUMBER);
				if (fabs(nval - (int)nval) < 1e-10)
					scf_options.maxiter = (int) nval;
				else
					errquit("in scf:maxiter: integer value expected in input!");
			}
			else if (!strcmp(sval, "diis")) {
				nextToken();
				if (ttype == TT_NUMBER) {
					if (fabs(nval - (int)nval) < 1e-10)
						scf_options.diisbas = (int) nval;
					else
						errquit("in scf:diis: integer value expected in input!");
				}
				else
					lexerPushBack();
				scf_options.diis = 1;
			}
			else if (!strcmp(sval, "nodiis"))
				scf_options.diis = 0;
			// conventional/direct
			else if (!strcmp(sval, "direct"))
				scf_options.direct = 1;
            else if (!strcmp(sval, "nodirect"))
                scf_options.direct = 0;
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
	
	if (!strcmp(sval, "overlap"))
		scf_options.print_1eov = doprint;
	else if (!strcmp(sval, "kinetic"))
		scf_options.print_1eke = doprint;
	else if (!strcmp(sval, "potential"))
		scf_options.print_1epe = doprint;
	else if (!strcmp(sval, "final vectors analysis"))
		scf_options.print_final_vectors = doprint;
	else if (!strcmp(sval, "eri") || !strcmp(sval, "ao2eints"))
		scf_options.print_2eri = doprint;
}

void scf_guess_inp()
{
	nextToken();
	if (ttype != TT_WORD)
		errquit("unexpected input in scf:guess");
	
	if (!strcmp(sval, "core"))
		scf_options.guess = GUESS_BARE;
	else if (!strcmp(sval, "eht"))
		scf_options.guess = GUESS_EHT;
	else
		errquit("wrong input in scf:guess");
}




