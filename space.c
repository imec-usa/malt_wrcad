// vi: ts=2 sts=2 sw=2 et tw=100
#include "space.h"
#include "config.h"
#include <stdio.h>

/* malt is in malt space (i.e. logspace) */
/* addpoint converts to physical space for spice and back to malt space for malt */
/* malt *reports* physical space */

double maltspace(double a, const Configuration *C, int i)
{
  /* conditional so log space can be turned off */
  if (C->params[i].logs)
    a = C->params[i].nominal / C->params[i].sigabs * (log(a / C->params[i].nominal) + 1);
  else
    a /= C->params[i].sigabs;
  return a;
}

double physspace(double a, const Configuration *C, int i)
{
  /* conditional so log space can be turned off */
  if (C->params[i].logs)
    a = C->params[i].nominal * exp(a * C->params[i].sigabs / C->params[i].nominal - 1);
  else
    a *= C->params[i].sigabs;
  return a;
}

/* double physaxes(double a,Configuration *C,int i){ */
/*   /\* conditional so log space can be turned off *\/ */
/*   if(C->params[i].logs) */
/*     a=exp(a*C->params[i].sigabs/C->params[i].nominal); */
/*   else */
/*     a *=C->params[i].sigabs; */
/*   return a; */
/* } */

/* double physcrit(double a,Configuration *C,int i){ */
/*   /\* conditional so log space can be turned off *\/ */
/*   if(C->params[i].logs) */
/*     a=exp(a*C->params[i].sigabs/C->params[i].nominal); */
/*   return a; */
/* } */

int initspace(Configuration *C, Space *S)
{
  int i;

  /* translate values from config file to internal representation */
  for (i = 0; C->num_params_all > i; ++i) {
    /* if nom_min and nom_max are not defined, fake it */
    if (!C->params[i].isnommin)
      C->params[i].nom_min = C->params[i].min;
    if (!C->params[i].isnommax)
      C->params[i].nom_max = C->params[i].max;
    /* if static, fake it (overrides nom_min and nom_max) */
    if (C->params[i].staticc) {
      C->params[i].nom_min = C->params[i].nominal;
      C->params[i].nom_max = C->params[i].nominal;
    }
    /* *** try 1 *** */
    /* copy from sigma (percent) to sigabs (absolute value) */
    /* for logspace, sigma is percent of nominal. for linear space, sigma is of unity */
    /* C->params[i].sigabs=C->params[i].logs? */
    /*   C->params[i].sigma*C->params[i].nominal/100.0: */
    /*   C->params[i].sigma*1.0/100.0; */
    /* *** try 2 *** */
    /* *** need sigma to be absolute, meaning sigabs=sigma *** */
    /* C->params[i].sigabs=C->params[i].sigma; */
    /* *** try 9 *** */
    /* make sure there is only one of sigma and sig_abs */
    /* and apply the result to sigma, in units of sig_abs */
    if (C->params[i].sigma == 0.0 && C->params[i].sigabs != 0.0)
      ; /* noop */
    else if (C->params[i].sigma != 0.0 && C->params[i].sigabs == 0.0)
      C->params[i].sigabs = C->params[i].sigma * C->params[i].nominal / 100.0;
    else {
      fprintf(stderr, "Just one of sigma or sig_abs must be non-zero in the config: param = %s\n",
              C->params[i].name);
      return 0;
    }
    /* convert from physical space to malt space */
    S[i].centerpnt = maltspace(C->params[i].nominal, C, i);
    C->params[i].min = maltspace(C->params[i].min, C, i);
    C->params[i].max = maltspace(C->params[i].max, C, i);
    C->params[i].nom_min = maltspace(C->params[i].nom_min, C, i);
    C->params[i].nom_max = maltspace(C->params[i].nom_max, C, i);
    /* define the hi/lo corners that will be used for corners-distributed params */
    /* in corner analysis. The +/-1.0 is the value of 1-sigma in malt-space */
    S[i].cornerhi = S[i].centerpnt + 1.0;
    S[i].cornerlo = S[i].centerpnt - 1.0;
  }
  return 1;
}
