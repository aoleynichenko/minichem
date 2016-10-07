#include <stdio.h>

#include "input.h"
#include "lexer.h"
#include "qslex.h"
#include "../scf/scf.h"
#include "../util/util.h"

#define MAX_INPUT 512

void prompt()
{
	printf("> ");
}

int getcmd(char *to, int max)
{
	int c, n = 0;
	char* p = to;

	while ((c = getchar()) != '\n' && c != EOF) {
		if (n < MAX_INPUT) /* ignore redundant symbols */
			*p++ = c;
		n++;
	}
	*p++ = '\0';
	n;
}

void qs_instruction()
{
	lx_nextToken();
	if (lx_ttype == QS_STRING) {
		if (strcmp(lx_sval, "q") == 0 || strcmp(lx_sval, "exit") == 0) {
			exit(0);
		}
		else if (strcmp(lx_sval, "credits") == 0) {
			printf("Alexander Oleynichenko, ao2310@yandex.ru\n");
		}
		else if (strcmp(lx_sval, "help") == 0) {
			printf("\n");
		}
		else if (strcmp(lx_sval, "nw") == 0) {
			lx_match(QS_QUOTE);
			printf("File name: %s\n", lx_sval);
			FILE *nwf = fopen(lx_sval, "r");
			printf("File: %s\n", lx_sval);
			if (!nwf)
				errquit("file in nw format not found");
			close(nwf);
			compute(lx_sval);
		}
	}
	else
		return;
}

void interactive_mode()
{
	extern char *source; // see lexer.c
	char input[MAX_INPUT];

	printf("Minichem 0.0 - Interactive mode\n");
	printf("Type 'help' for help and 'q' to exit program\n\n");

	// init internal structures
	calc_info_defaults();
	scf_init();

	while (1) {
		prompt();
		if (getcmd(input, MAX_INPUT) == -1)
			errquit("while reading input line in the interactive mode");
		qs_lex(input);
		qs_instruction();
	}
}
