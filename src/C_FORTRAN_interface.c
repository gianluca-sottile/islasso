/* $Id: C_FORTRAN_interface.c 313 2015-09-16 20:20:04Z mmaechler $
*
*  wrapper for calling R's random number generator from
*  the original FORTRAN code
*
*/

#include "islasso.h"

double F77_SUB(pnm)(double *x, double *mu, double *sigma){ return pnorm(*x, *mu, *sigma, 1, 0); }
double F77_SUB(dnm)(double *x, double *mu, double *sigma){ return dnorm(*x, *mu, *sigma, 0); }
double F77_SUB(qnm)(double *x, double *mu, double *sigma){ return qnorm(*x, *mu, *sigma, 1, 0); }
