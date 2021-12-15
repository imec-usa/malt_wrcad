
#include "stat_math.h"
#include <math.h>

/* all internal functions listed here. the external ones are in the header file */
double factln(int);
double beta(double, double);
double gammp(double, double);
double gammq(double, double);
void gser(double *, double, double, double *);
void gcf(double *, double, double, double *);

double hypsphere_deviate(double *p, long *idum, double r0, int dim)
{
  /* returns the coordinates of a random point uniformly distributed in a unit hypersphere */
  /* r0 can be used to limit returned point to a shell of radius in the range r0-to-1 */
  int j;
  double grad, rad, map;

  do {
    for (grad = 0, j = 0; j < dim; j++) {
      p[j] = gauss_deviate(idum);
      grad += p[j] * p[j];
    }
    grad = sqrt(grad);
    /* the mapping--only empirical until proven! */
    rad = erf(grad * 2.0 / M_PI);
    map = rad / grad;
    /* the following conditional looks terrible, but is actually quite efficient */
    /* therefore, don't really need to store values in bins */
  } while ((rad) < r0);
  for (j = 0; j < dim; j++)
    p[j] *= map;
  return rad;
}

double gammln(double xx)
/* Returns the value ln Gamma(xx) for xx > 0. */
{
  /* Internal arithmetic will be done in double precision, a nicety that you can omit if five figure
   */
  /* accuracy is good enough. */
  double x, y, tmp, ser;
  static double cof[6] = {76.18009172947146,  -86.50532032941677,    24.01409824083091,
                          -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5};
  int j;
  y = x = xx;
  tmp = x + 5.5;
  tmp -= (x + 0.5) * log(tmp);
  ser = 1.000000000190015;
  for (j = 0; j <= 5; j++)
    ser += cof[j] / ++y;
  return -tmp + log(2.5066282746310005 * ser / x);
}

double factrl(int n)
/* Returns the value n! as a floating-point number. */
{
  static int ntop = 4;
  static double a[170] = {1.0, 1.0, 2.0, 6.0, 24.0}; /* Fill in table only as required. */
  int j;

  if (n < 0)
    nrerror("Negative factorial in routine factrl");
  if (n > 169)
    return exp(gammln(n + 1.0));
  /* Larger value than size of table is required. Actually, this big a value is going to overflow */
  /* the double-precision many computers, but no harm in trying. */
  while (ntop < n) { /* Fill in table up to desired value. */
    j = ntop++;
    a[ntop] = a[j] * ntop;
  }
  return a[n];
}

double bico(int n, int k)
/* Returns the binomial coefficient (n choose k) as a floating point number */
{
  return floor(0.5 + exp(factln(n) - factln(k) - factln(n - k)));
  /* The floor function cleans up roundoff error for smaller values of n and k. */
}

double factln(int n)
/* Returns ln(n!). */
{
  static double a[171]; /* A static array is automatically initialized to zero. */
  if (n < 0)
    nrerror("Negative factorial in routine factln");
  if (n <= 1)
    return 0.0;
  if (n <= 170)
    return a[n] ? a[n] : (a[n] = gammln(n + 1.0)); /* In range of table. */
  else
    return gammln(n + 1.0); /* Out of range of table. */
}

double beta(double z, double w)
/* Returns the value of the beta function B(z, w). */
{
  return exp(gammln(z) + gammln(w) - gammln(z + w));
}

double gammp(double a, double x)
/* Returns the incomplete gamma function P(a, x). */
{
  double gamser = 0, gammcf, gln;

  if (x < 0.0 || a <= 0.0)
    nrerror("Invalid arguments in routine gammp");
  if (x < (a + 1.0)) { /* Use the series representation. */
    gser(&gamser, a, x, &gln);
    return gamser;
  } else { /* Use the continued fraction representation. */
    gcf(&gammcf, a, x, &gln);
    return 1.0 - gammcf; /* and take its complement. */
  }
}

double gammq(double a, double x)
/* Returns the incomplete gamma function Q(a, x)=1-P(a,x). */
{
  double gamser = 0, gammcf, gln;

  if (x < 0.0 || a <= 0.0)
    nrerror("Invalid arguments in routine gammq");
  if (x < (a + 1.0)) { /* Use the series representation */
    gser(&gamser, a, x, &gln);
    return 1.0 - gamser; /* and take its complement. */
  } else {               /* Use the continued fraction representation. */
    gcf(&gammcf, a, x, &gln);
    return gammcf;
  }
}

#define ITMAX 100
#define EPS 3.0e-28

void gser(double *gamser, double a, double x, double *gln)
/* Returns the incomplete gamma function P(a, x) evaluated by its series representation as gamser.
 */
/* Also returns ln Gamma(a) gln. */
{
  int n;
  double sum, del, ap;

  *gln = gammln(a);
  if (x <= 0.0) {
    if (x < 0.0)
      nrerror("x less than 0 in routine gser");
    *gamser = 0.0;
    return;
  } else {
    ap = a;
    del = sum = 1.0 / a;
    for (n = 1; n <= ITMAX; n++) {
      ++ap;
      del *= x / ap;
      sum += del;
      if (fabs(del) < fabs(sum) * EPS) {
        *gamser = sum * exp(-x + a * log(x) - (*gln));
        return;
      }
    }
    nrerror("a too large, ITMAX too small in routine gser");
    return;
  }
}

#undef ITMAX
#undef EPS

#define ITMAX 100      /* Maximum allowed number of iterations. */
#define EPS 3.0e-28    /* Relative accuracy. */
#define FPMIN 1.0e-300 /* Number near the smallest representable doubleing-point number. */

void gcf(double *gammcf, double a, double x, double *gln)
/* Returns the incomplete gamma function Q(a, x) evaluated by its continued fraction represen- */
/* tation as gammcf. Also returns ln Gamma(a) as gln */
{
  int i;
  double an, b, c, d, del, h;

  *gln = gammln(a);
  b = x + 1.0 - a; /* Set up for evaluating continued fraction */
  /* by modified Lentz's method (5.2) with b0 = 0. */
  c = 1.0 / FPMIN;
  d = 1.0 / b;
  h = d;
  for (i = 1; i <= ITMAX; i++) { /* Iterate to convergence. */
    an = -i * (i - a);
    b += 2.0;
    d = an * d + b;
    if (fabs(d) < FPMIN)
      d = FPMIN;
    c = b + an / c;
    if (fabs(c) < FPMIN)
      c = FPMIN;
    d = 1.0 / d;
    del = d * c;
    h *= del;
    if (fabs(del - 1.0) < EPS)
      break;
  }
  if (i > ITMAX)
    nrerror("a too large, ITMAX too small in gcf");
  *gammcf = exp(-x + a * log(x) - (*gln)) * h; /* Put factors in front. */
}

#undef ITMAX
#undef EPS
#undef FPMIN

double nr_erf(double x)
/* Returns the error function erf(x). */
{
  return x < 0.0 ? -gammp(0.5, x * x) : gammp(0.5, x * x);
}

double nr_erfc(double x)
/* Returns the complementary error function erfc(x). */
{
  return x < 0.0 ? 1.0 + gammp(0.5, x * x) : gammq(0.5, x * x);
}

double gauss_integral(double r, int n)
/* Returns the value of the n-dimensional gaussian integral. */
{
  return gammp(n / 2.0, r * r / 2.0);
}

double gauss_integral_c(double r, int n)
/* Returns the value of the complementary n-dimensional gaussian integral. */
/* i.e., 1-gauss_integral */ { return gammq(n / 2.0, r * r / 2.0); }

double gauss_deviate(long *idum)
{
  /* Returns a normally distributed deviate with zero mean and unit variance, */
  /* using uniform_deviate() as the source of uniform deviates. */

  static int iset = 0;
  static double gset;

  double fac, r, v1, v2;

  if (iset == 0) {
    do {
      /* pick two uniform numbers in the square extending */
      /* from -1.0 to +1.0 in each direction */
      v1 = 2.0 * uniform_deviate(idum) - 1.0;
      v2 = 2.0 * uniform_deviate(idum) - 1.0;
      r = v1 * v1 + v2 * v2;
    } while (r >= 1.0);

    /* Now make the Box-Muller transformation to get two normal deviates. */
    /* Return one and save the other for next time */
    fac = sqrt(-2.0 * log(r) / r);

    gset = v1 * fac;
    iset = 1;
    return v2 * fac;
  } else {
    iset = 0;
    return gset;
  }
}

#define IM 714025
#define IA 4096
#define IC 54773
#define SHUFFLE_LENGTH 98

double uniform_deviate(long *idum)
{
  /* Return a uniform random deviate between 0.0 and 1.0. */
  /* Set idum to any negative value to initialize or reinitialize the sequence. */

  static long iy, shuff_tab[SHUFFLE_LENGTH];
  static int iff = 0;

  int j;

  if ((*idum < 0) || (iff == 0)) {
    iff = 1;
    if ((*idum = (IC - (*idum)) % IM) < 0)
      *idum = -(*idum);

    /* initialize the shuffle table */
    for (j = 1; j <= SHUFFLE_LENGTH - 1; j++) {
      *idum = (IA * (unsigned long)(*idum) + IC) % IM;
      shuff_tab[j] = (*idum);
    }
    *idum = (IA * (unsigned long)(*idum) + IC) % IM;
    iy = (*idum);
  }
  j = 1 + (SHUFFLE_LENGTH - 1) * iy / IM;

  if ((j > SHUFFLE_LENGTH - 1) || (j < 0))
    nrerror("Error in random number generator uniform_deviate.");

  iy = shuff_tab[j];

  *idum = (IA * (unsigned long)(*idum) + IC) % IM;
  shuff_tab[j] = *idum;

  return (double)iy / IM;
}

#undef IM
#undef IA
#undef IC
#undef SHUFFLE_LENGTH
