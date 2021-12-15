// vi: ts=2 sts=2 sw=2 et tw=100

#ifndef MARGINS
#define MARGINS

#include "config.h"

typedef struct marg {
  double angle;
  double x;
  double y;
} Marg;

int call_trace(Configuration *);
int call_marg(Configuration *);
int shmoo(Configuration *);
int margins2(Configuration *);
void shift(double *, double **, double **, int, int, int, int);

#endif
