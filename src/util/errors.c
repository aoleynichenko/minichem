#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "util.h"

void errquit(char *errmessage)
{
	int rank;
	time_t t = time(0);
	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0) {
		printf("------------------------------------------------------------\n");
		printf("Error: %s\n", errmessage);
		printf("Abort at %s", asctime(localtime(&t)));
	}
	MPI_Abort(MPI_COMM_WORLD, 1);
}
