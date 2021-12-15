// vi: ts=2 sts=2 sw=2 et tw=100
#ifndef NUMERICAL
#define NUMERICAL

void nrerror(const char *);
extern int *ivector(int nl, int nh);
extern double *vector(int nl, int nh);
extern double **matrix(int nrl, int nrh, int ncl, int nch);
extern void free_ivector(int *v, int nl, int nh);
extern void free_vector(double *v, int nl, int nh);
extern void free_matrix(double **m, int nrl, int nrh, int ncl, int nch);

#endif
