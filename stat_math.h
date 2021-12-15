// vi: ts=2 sts=2 sw=2 et tw=100

#ifndef STAT_MATH
#define STAT_MATH

#include "numerical.h"

/* Returns the value n! as a floating-point number. */
double factrl(int);
/* Returns the binomial coecient (n choose k) as a floating point number */
double bico(int, int);
/* Returns the error function erf(x). */
double nr_erf(double);
/* Returns the complementary error function erfc(x). */
double nr_erfc(double);
/* value of the n-dimensional gaussian integral. */
double gauss_integral(double, int);
/* value of the n-dimensional complementary gaussian integral. */
double gauss_integral_c(double, int);
/* random number between 0 and 1 */
double uniform_deviate(long *);
/* random number with zero mean and unit varience */
double gauss_deviate(long *);
/* returns the coordinates of a random point uniformly distributed in a unit hypersphere */
/* r0 can be used to limit returned point to a shell of radius in the range r0-to-1 */
double hypsphere_deviate(double *, long *, double, int);
double gammln(double);

#endif
