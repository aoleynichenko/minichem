#include <stdio.h>

#include "util.h"


void echo(char *filename)
{
	int c;
	FILE *f = fopen(filename, "r");
	
	if (!f) {
		char buf[300];
		sprintf(buf, "cannot open file '%s'", filename);
		errquit(buf);
	}
	printf(">>> %s\n", filename);
	while ((c = fgetc(f)) != EOF)
		putchar(c);
	fclose(f);
	line_separator();
}

void line_separator()
{
	printf("------------------------------------------------------------\n");
}
