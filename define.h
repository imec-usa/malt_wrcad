
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
