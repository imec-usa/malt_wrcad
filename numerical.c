/* For the numerical recipes routines */
#include <stdio.h>
#include <stdlib.h>

void nrerror(const char *error_text)
{
  fprintf(stderr, "malt: Numerical Error\n%s\n", error_text);
  exit(EXIT_FAILURE);
}

void memerror()
{
  fprintf(stderr, "malt: Out of system memory\n");
  exit(EXIT_FAILURE);
}

int *ivector(int nl, int nh)
{
  int *v;
  v = malloc((nh - nl + 1) * sizeof *v);  // mem:bloodsucking
  if (!v)
    memerror();
  return v - nl;
}

double *vector(int nl, int nh)
{
  double *v;
  v = malloc((nh - nl + 1) * sizeof *v);  // mem:guaka
  if (!v)
    memerror();
  return v - nl;
}

double **matrix(int nrl, int nrh, int ncl, int nch)
{
  int i;
  double **m;

  m = malloc((nrh - nrl + 1) * sizeof *m);  // mem:unsurmountableness
  if (!m)
    memerror();
  m -= nrl;

  for (i = nrl; i <= nrh; i++) {
    m[i] = malloc((nch - ncl + 1) * sizeof *m[i]);  // mem:interantagonism
    if (!m[i])
      memerror();
    m[i] -= ncl;
  }
  return m;
}

void free_ivector(int *v, int nl, int nh)
{
  (void)nh;
  free(v + nl);  // mem:bloodsucking
}

void free_vector(double *v, int nl, int nh)
{
  (void)nh;
  free(v + nl);  // mem:guaka
}

void free_matrix(double **m, int nrl, int nrh, int ncl, int nch)
{
  int i;

  (void)nch;
  for (i = nrh; i >= nrl; i--) {
    free(m[i] + ncl);  // mem:interantagonism
  }
  free(m + nrl);  // mem:unsurmountableness
}
