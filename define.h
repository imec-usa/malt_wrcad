// vi: ts=2 sts=2 sw=2 et tw=100

#ifndef DEFINE
#define DEFINE

// typedef struct config Configuration;
#include "config.h"

typedef struct data {
  double *t;
  double tstep;
  double **x;
  double **upper;
  double **lower;
  int length;
  int dtlength;
} Data;

int call_def(Configuration *);

#endif
