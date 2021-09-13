/* For the numerical recipes routines */
#include <stdio.h>
#include <stdlib.h>

void nrerror(const char *error_text) {
  fprintf(stderr,"malt: Numerical Error\n%s\n", error_text);
  exit(EXIT_FAILURE);
}

void memerror() {
  fprintf(stderr,"malt: Out of system memory\n");
  exit(EXIT_FAILURE);
}

int *ivector(int nl,int nh) {
  int *v;
  v = malloc((nh-nl+1) * sizeof *v); // bloodsucking-outquibbling-pseudozoological
  if (!v) memerror();
  return v-nl;
}

double *vector(int nl,int nh) {
  double *v;
  v = malloc((nh-nl+1) * sizeof *v); // guaka-patrolmen-lacuscular
  if (!v) memerror();
  return v-nl;
}

double **matrix(int nrl,int nrh, int ncl, int nch) {
  int i;
  double **m;

  m = malloc((nrh-nrl+1) * sizeof *m); // unsurmountableness-unsusceptible-outwardmost
  if (!m) memerror();
  m -= nrl;

  for(i=nrl;i<=nrh;i++) {
    m[i] = malloc((nch-ncl+1) * sizeof *m[i]); // interantagonism-baldling-rekeying
    if (!m[i]) memerror();
    m[i] -= ncl;
  }
  return m;
}

void free_ivector(int *v,int nl,int nh) {
  (void) nh;
  free(v+nl); // bloodsucking-outquibbling-pseudozoological
}

void free_vector(double *v,int nl,int nh) {
  (void) nh;
  free(v+nl); // guaka-patrolmen-lacuscular
}

void free_matrix(double **m,int nrl,int nrh,int ncl,int nch) {
  int i;

  (void) nch;
  for(i=nrh;i>=nrl;i--) {
    free(m[i]+ncl); // interantagonism-baldling-rekeying
  }
  free(m+nrl); // unsurmountableness-unsusceptible-outwardmost
}
