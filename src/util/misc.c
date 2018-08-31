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
