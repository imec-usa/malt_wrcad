
#ifndef MARG_OPT_YIELD
#define MARG_OPT_YIELD

#include "config.h"
#include "space.h"

typedef struct {
  short flag;
  unsigned short *points;
  double b;
  double *a;
} Plane;

typedef unsigned int corner_t;

void     makeiter(Configuration*, char);
int      checkiter(Configuration*);
int      tmargins(Configuration*, const Space *);
int      margins(Configuration *C, const Space *S, double *prhi, double *prlo)
  __attribute__((nonnull));
double addpoint_corners(const Configuration *C, const Space *S, corner_t *cornmin, double *pr,
                        const double *pc, const double *direction);
Plane  **plane_malloc(Plane **, int *, int, int);
void     plane_free(Plane **, int);
double **margpnts_malloc(double **, int *, int, int);
void     margpnts_free(double **, int);
void     intpickpnts(short *, Configuration *, const Space *, int, Plane **, double **, int *, int);
int      makeaplane(short *, Configuration *, const Space *, Plane **, double **, int *, int);
int      center(Configuration *, Space *, Plane **, int *, int, double *);
int      findface(int, double *, Plane **, int *, int, int);
double   det_dim(double **, int);
int      simplx(double *const *, int, int, int *, int *);
int      hull_dice(Configuration *, double *, Plane **, int);

#endif
