#include <mpi.h>
#include <stdio.h>
#include <time.h>

#include "util.h"

void help()
{
	printf("Minichem 0.0\n");
	printf("Usage: [mpirun -np N] minichem [options] <input-files>\n");
	printf("where options include:\n");
	printf("    -h              print this help message\n");
	printf("    -v              print version and quit\n");
	printf("    -i              run in interactive mode\n");
	printf("    -o <file-name>  specify output file\n");
	printf("    -noecho         do not print input files\n");
	printf("for questions, ao2310@yandex.ru\n");
	printf("\n");
}

void print_header()
{
	time_t t = time(0);
	int size;
	
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	printf("Minichem 0.0\n");
	printf("By Alexander Oleynichenko, 2016\n");
	printf("Build:      %d\n", buildno());
	printf("Build date: %s %s\n", __DATE__, __TIME__);
	
	#if defined __ICC
	printf("Compiler:   Intel C Compiler %d\n", __ICC);
	#elif defined __GNUC__
	printf("Compiler:   gcc %d.%d.%d\n", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
	#else
	printf("Compiler:   undetected\n");
	#endif
	
	printf("Date:       %s", asctime(localtime(&t)));
	printf("Cores:      %d\n", size);
}

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
