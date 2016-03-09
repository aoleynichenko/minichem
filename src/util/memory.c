#include <mpi.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

#include "util.h"

int memavail = 200*1024*1014;  /* 200 mb per thread */;
int allocated = 0;
int total_allocated = 0;
int total_freed = 0;
int max_allocated = 0;

void setmemavail(int bytes)
{
	memavail = bytes;
}

/* Quantum chemistry allocation */
void *qalloc(size_t bytes)
{
	void *p = NULL;
	
	allocated += bytes;
	total_allocated += bytes;
	if (allocated > memavail) {
		char buf[512];
		sprintf(buf, "while dynamic memory allocation: not enough memory!\n\
Requested memory: %Lu bytes, available free memory: %Lu bytes.\n\
Available memory per thread: %Lu bytes.\n\
You may use the 'memory' directive to raise the memory usage threshold.\n\
For more information, see source code (util/memory.c/qalloc).\n",
			bytes, memavail-(allocated-bytes), memavail);
		errquit(buf);
		return NULL;
	}
	
	p = malloc(bytes);
	if (!p) {
		char buf[512];
		sprintf(buf, "while dynamic memory allocation: \
void *malloc(size_t) returned NULL pointer.!\n\
Requested memory: %Lu bytes.\n\
Available memory per thread: %Lu bytes.\n\
For more information, see source code (util/memory.c/qalloc).\n",
			bytes, memavail);
		errquit(buf);
		return NULL;
	}
	
	/* memory allocated successfully */
	if (allocated > max_allocated)
		max_allocated = allocated;
	return p;
}

/* The best way to free memory in minichem */
void qfree(void *p, size_t n)
{
	free(p);
	allocated -= n;
	total_freed += n;
}

/* print dynamic memory usage statistics */
void memstats()
{
	int rank;
	int totalloc = 0;
	int totfreed = 0;
	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Reduce(&total_allocated, &totalloc, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&total_freed, &totfreed, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	if (rank == 0) {
		printf("\nAllocated: %d bytes\nFreed:     %d bytes\n", totalloc, totfreed);
		printf("Max usage: %d bytes\n", max_allocated);
	}
}










