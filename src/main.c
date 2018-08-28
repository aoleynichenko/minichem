/* Minichem - simple quantum chemistry program.
 * 
 * A. Oleynichenko
 * 2016 - 2018
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "util.h"
#include "input.h"

void print_usage();
void print_header();

int main(int argc, char **argv)
{
	int i;
	int nargc = argc;
	int ninp = 0;
	int noecho = 0;
	char **argvp = argv;
	char inputs[10][256];
	int rank, size;
	double t1, t2;
	
	/* MPI initialization */
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	t1 = MPI_Wtime();
	
	/* parse command line */
	if (rank == 0) {
		if (nargc == 1)
			errquit("no input files specified, -h for options");
		while (--nargc) {
			argvp++;
			if (!strcmp(*argvp, "-h") ||
				!strcmp(*argvp, "--h") ||
				!strcmp(*argvp, "-help") ||
				!strcmp(*argvp, "-?")) {
				print_usage();
				MPI_Finalize();
				return 0;
			}
			else if (!strcmp(*argvp, "-v")) {
				printf("Minichem 0.0\n\n");
				MPI_Finalize();
				return 0;
			}
			else if (!strcmp(*argvp, "-noecho")) {
				noecho = 1;
			}
			else if (!strcmp(*argvp, "-o")) {
				if (nargc > 1) {
					char output[256];
					argvp++;
					nargc--;
					strncpy(output, *argvp, 256);
					freopen(output, "w", stdout);
				}
				else {
					printf("-o option: no file name specified\n");
					MPI_Finalize();
					return 1;
				}
			}
			else {
				strncpy(inputs[ninp++], *argvp, 256);
			}
		}
		
		print_header();
		printf("Arguments:  ");
		/* print command line */
		for (i = 0; i < argc; i++)
			printf("%s ", argv[i]);
		printf("\n");
		line_separator();
	}
	
	/* 0 shares file names with others */
	MPI_Bcast(&ninp, 1, MPI_INT, 0, MPI_COMM_WORLD);
	for (i = 0; i < ninp; i++) {
		if (!noecho && rank == 0)
			echo(inputs[i]);
		MPI_Bcast(&inputs[i], 256, MPI_CHAR, 0, MPI_COMM_WORLD);
		compute(inputs[i]);  /* call input module */
		//MPI_Barrier(MPI_COMM_WORLD);
	}
	
	memstats();  // memory usage statistics
	if (rank == 0) {
		t2 = MPI_Wtime();
		printf("Total time: %.2f sec\n", t2 - t1);
	}
	
	MPI_Finalize();
}


void print_usage()
{
	printf("Minichem 1.0\n");
	printf("Usage: [mpirun -np N] minichem.x [options] <input-files>\n");
	printf("where options include:\n");
	printf("    -h              print this help message\n");
	printf("    -v              print version and quit\n");
	printf("    -o <file-name>  specify output file\n");
	printf("    -noecho         do not print input files\n");
	printf("for questions, alexvoleynichenko@gmail.com\n");
	printf("\n");
}


/***********************************************************************
 * print_header
 * 
 * Print version, compiler & platform info, settings for the parallel
 * execution.
 **********************************************************************/
void print_header()
{
	time_t t = time(0);
	int size;
	
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	printf("Minichem 1.0\n");
	printf("By Alexander Oleynichenko, 2018\n");
	printf("Build date:    %s %s\n", __DATE__, __TIME__);
	
	#if defined __ICC
	printf("Compiler:      Intel C Compiler %d\n", __ICC);
	#elif defined __GNUC__
	printf("Compiler:      gcc %d.%d.%d\n", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
	#else
	printf("Compiler:      undetected\n");
	#endif
	
	printf("Date:          %s", asctime(localtime(&t)));
	printf("# MPI threads: %d\n", size);
}





