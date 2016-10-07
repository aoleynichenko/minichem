#include <stdio.h>

#include "input.h"
#include "lexer.h"
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
		printf("%s\n", input);
		source = input;
		printf("Source = >>>%s<<<\n", input);
		print_tokens();
	}
}