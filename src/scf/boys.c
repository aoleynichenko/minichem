#include <math.h>

#include "ints.h"

/* Very naive implementation of the Boys function Fn(x) */
double boys(int n, double x)
{
	if (x < 0.1) {    /* Taylor expansion */
		return 1.0/(2.0*n+1.0)
			   - x/(2.0*n+3.0)
			   + 0.5*x*x/(2.0*n+5.0)
			   - x*x*x/((2.0*n+7.0)*6.0)
			   + x*x*x*x/((2.0*n+9.0)*24.0)
			   - x*x*x*x*x/((2.0*n+11.0)*120.0)
			   + x*x*x*x*x*x/((2.0*n+13.0)*720.0);
	}
	if (n == 0) {
		return 0.5*sqrt(M_PI/x)*erf(sqrt(x));
	}
	else {   /* downward recursion - very bad! */
		return 0.5*((2.0*n-1.0)*boys(n-1, x)-exp(-x))/x;
	}
}
