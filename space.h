#include <math.h>

#ifndef SPACE
#define SPACE

#include "config.h"

typedef struct space {
  double cornerhi;
  double cornerlo;
  double centerpnt;
} Space;

double maltspace(double, const Configuration *, int);
double physspace(double, const Configuration *, int);
int    initspace(Configuration *, Space *);

#endif
