// vi: ts=2 sts=2 sw=2 et tw=100
/* calculate margins */
#include "margins.h"
#include "call_spice.h"
#include "config.h"
#include "gplot.h"
#include "malt.h"
#include "marg_opt_yield.h"
#include "space.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/* prints a file for each param pair, like margins2 */
int shmoo(Configuration *C)
{
  (void)C;
  printf("Not implemented at this time.\n");
  return 1;
}

static int newpoint(Configuration *C, const Space *S, double *pc, double *direction, double *angle,
                    double **pointsx, double **pointsy, int j, int i)
{
  double *pr = malloc(C->num_params * sizeof *pr);
  int ok = 0;
  /* param_x & param_y says which parameter number to do margins on */
  direction[C->_2D[i].param_x] = sin(angle[j]);
  direction[C->_2D[i].param_y] = cos(angle[j]);
  if (0.0 == addpoint_corners(C, S, NULL, pr, pc, direction)) {
    fprintf(stderr, "Nominal parameter values failed\n");
    goto fail;
  }
  /* save points */
  pointsx[i][j] = pr[C->_2D[i].param_x];
  pointsy[i][j] = pr[C->_2D[i].param_y];

  ok = 1;
fail:
  free(pr);
  return ok;
}

/* prints a file for each param pair */
int margins2(Configuration *C)
{
  int i, j, k;
  double *angle, **pointsx, **pointsy, dist;
  int all_good = 1, iterfile = 1;
  int num_all, jnew;
  double dnew;
  Space *S = malloc(C->num_params_all * sizeof *S);

  /* initialize */
  if (!initspace(C, S))
    return 0;
  /* create the iterate file */
  makeiter(C, '2');
  /* memory allocation */
  if (C->options._2D_iter < 3) {
    C->options._2D_iter = 3;
  }
  angle = malloc((C->options._2D_iter + 1) * sizeof *angle);  // mem:listener
  pointsx = malloc((C->num_2D) * sizeof *pointsx);            // mem:caracolite
  pointsy = malloc((C->num_2D) * sizeof *pointsy);            // mem:juvenilities
  for (i = 0; C->num_2D > i; ++i) {
    pointsx[i] = malloc((C->options._2D_iter + 1) * sizeof *pointsx[i]);  // mem:anacahuite
    pointsy[i] = malloc((C->options._2D_iter + 1) * sizeof *pointsy[i]);
  }  // mem:nonverticalness
  /* create pname file */
  pname(C);
  double *pc = malloc((C->num_params + C->num_params_corn) * sizeof *pc);  // mem:restrip
  double *direction = malloc(C->num_params * sizeof *direction);           // mem:overproficient
  for (i = 0; C->num_params > i; ++i)
    pc[i] = S[i].centerpnt;
  /* assign param_x and param_y */
  /* if the parameter name does not exist, assign -1 */
  for (i = 0; C->num_2D > i; ++i) {
    C->_2D[i].param_x = C->_2D[i].param_y = -1;
    for (j = 0; C->num_params > j; ++j) {
      if (!strcmp(C->_2D[i].name_x, C->params[j].name))
        C->_2D[i].param_x = j;
      else if (!strcmp(C->_2D[i].name_y, C->params[j].name))
        C->_2D[i].param_y = j;
    }
    /* check that there are no unknown parameters */
    if (C->_2D[i].param_x == -1) {
      fprintf(stderr, "malt: Undefined param_x name: %s\n", C->_2D[i].name_x);
      all_good = 0;
      goto cleanup;
    }
    if (C->_2D[i].param_y == -1) {
      fprintf(stderr, "malt: Undefined param_y name: %s\n", C->_2D[i].name_y);
      all_good = 0;
      goto cleanup;
    }
  }
  /* loop through each pair */
  for (i = 0; i < C->num_2D; i++) {
    /* 1) for x/2 iterations, do equal angles */
    num_all = C->options._2D_iter / 2;
    /* do at least eight points in this iteration */
    if (num_all < 8) {
      num_all = 8;
      if (num_all > C->options._2D_iter)
        num_all = C->options._2D_iter;
    }
    for (j = 0; j < num_all && iterfile; j++) {
      /* step through angle */
      angle[j] = 2 * M_PI * j / (double)num_all;
      /* initialize */
      for (k = 0; C->num_params > k; ++k)
        direction[k] = 0;
      if (!(all_good = newpoint(C, S, pc, direction, angle, pointsx, pointsy, j, i))) {
        goto cleanup;
      }
      /* check for the iterate.2 file */
      if (!checkiter(C))
        iterfile = 0;
    }
    /* duplicate the first point */
    angle[j] = 2 * M_PI;
    pointsx[i][j] = pointsx[i][0];
    pointsy[i][j] = pointsy[i][0];

    /* 2) for x/2 iterations, bisect angles corresponding to largest distances (pick one by one) */
    /* loop C->options._2D_iter-num_all times */
    for (; C->options._2D_iter > num_all; ++num_all) {
      /* find biggest distance */
      for (dist = 0, jnew = j = 0; num_all > j; ++j) {
        dnew = sqrt(pow(pointsx[i][j] - pointsx[i][j + 1], 2) +
                    pow(pointsy[i][j] - pointsy[i][j + 1], 2));
        if (dnew > dist) {
          dist = dnew;
          jnew = j;
        }
      }
      /* make room for new point */
      shift(angle, pointsx, pointsy, 1, num_all, jnew, i);
      /* add the new point */
      angle[jnew + 1] = (angle[jnew] + angle[jnew + 2]) / 2;
      if (!(all_good = newpoint(C, S, pc, direction, angle, pointsx, pointsy, jnew + 1, i))) {
        goto cleanup;
      }
      /* check for the iterate.2 file */
      if (!checkiter(C))
        iterfile = 0;
    }
  }
  /* print the file */
  /* names */
  lprintf(C, "# ");
  for (i = 0; i < C->num_2D; i++) {
    lprintf(C, "(%-8.8s %-8.8s)", C->params[C->_2D[i].param_x].name,
            C->params[C->_2D[i].param_y].name);
  }
  // clang-format off
    lprintf(C, "\n");  // FIXME: should be in the for loop by indentation? ~ntj
  /* points */
  for (j = 0; j <= C->options._2D_iter; j++) {
    for (i = 0; i < C->num_2D; i++) {
      lprintf(C, "  %8.3e %8.3e",
              physspace(pointsx[i][j], C, C->_2D[i].param_x),
              physspace(pointsy[i][j], C, C->_2D[i].param_y));
    }
    lprintf(C, "\n");
  }
  /* nominal points */
  lprintf(C, "\n");
  for (i = 0; i < C->num_2D; i++) {
    lprintf(C, "  %8.3e %8.3e",
            physspace(S[C->_2D[i].param_x].centerpnt, C, C->_2D[i].param_x),
            physspace(S[C->_2D[i].param_y].centerpnt, C, C->_2D[i].param_y));
  }
  // clang-format on
  lprintf(C, "\n");
  /* flush the data file */
  /* term_file_flush(C); */
  /* gnuplot the data */
  plot2(C, S);
  /* remove temp files */
cleanup:
  unlink(C->file_names.iter);
  unlink_pname(C);
  // free(S); // TODO: needed?
  free(angle);  // mem:listener
  for (i = 0; C->num_2D > i; ++i) {
    free(pointsx[i]);  // mem:anacahuite
    free(pointsy[i]);  // mem:nonverticalness
  }
  free(pointsx);    // mem:caracolite
  free(pointsy);    // mem:juvenilities
  free(pc);         // mem:restrip
  free(direction);  // mem:overproficient
  return all_good;
}

void shift(double *angle, double **pointsx, double **pointsy, int num_new, int num_all, int j,
           int i)
{
  int k;

  for (k = num_all; j < k; --k) {
    angle[k + num_new] = angle[k];
    pointsx[i][k + num_new] = pointsx[i][k];
    pointsy[i][k + num_new] = pointsy[i][k];
  }
}

int call_marg(Configuration *C)
{
  int all_good;
  double *prhi = malloc(C->num_params * sizeof *prhi);
  double *prlo = malloc(C->num_params * sizeof *prlo);

  /* initialize */
  Space *S = malloc(C->num_params_all * sizeof *S);
  if (!initspace(C, S))
    return 0;
  /* create pname file */
  pname(C);
  /* do it */
  if (!(all_good = margins(C, S, prhi, prlo)))
    goto cleanup;

  /* clean up temporary files */
cleanup:
  free(prhi);
  free(prlo);
  free(S);
  unlink_pname(C);
  return all_good;
}

int call_trace(Configuration *C)
{
  int all_good;

  /* initialize */
  Space *S = malloc(C->num_params_all * sizeof *S);
  if (!initspace(C, S))
    return 0;
  /* create pname file */
  pname(C);
  /* do it */
  all_good = tmargins(C, S);
  /* clean up temporary files */
  /*cleanup:*/
  unlink_pname(C);
  return all_good;
}
