/***********************************************************************
 * misc.c
 * ======
 * 
 * Miscellaneous utility functions.
 * 
 **********************************************************************/

#include <math.h>

// signum
int sgn(double x)
{
	if (fabs(x) < 1e-14) return 0;
	else if (x < 0)  return -1;
	else /* x > 0 */ return 1;
}


// distance between two atoms
double distance(double *A, double *B)
{
	return sqrt((A[0] - B[0])*(A[0] - B[0]) +
				(A[1] - B[1])*(A[1] - B[1]) +
				(A[2] - B[2])*(A[2] - B[2]));
}
