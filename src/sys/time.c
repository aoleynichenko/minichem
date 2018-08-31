/***********************************************************************
 * time.c
 * ========
 * 
 * Interface to the system-dependent functions for time measurements.
 * 
 * 2018 Alexander Oleynichenko
 **********************************************************************/

#include <sys.h>

#if __has_include(<mpi.h>)
#include <mpi.h>
double abs_time()
{
	return MPI_Wtime();
}

#else  /* no MPI */
#include <time.h>
double abs_time()
{
	return ((double) clock()) / CLOCKS_PER_SEC;
}

#endif
