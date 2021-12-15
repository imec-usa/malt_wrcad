// vi: ts=2 sts=2 sw=2 et tw=100
/* optimization subroutines */
#include "marg_opt_yield.h"
#include "call_spice.h"
#include "config.h"
#include "malt.h"
#include "numerical.h"
#include "space.h"
#include <assert.h>
#include <math.h>
#include <pthread.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/wait.h>
#include <unistd.h>

void makeiter(Configuration *C, char func)
{
  FILE *fp;
  /* create the iterate file */
  sprintf(C->file_names.iter, "%s.iterate.%c", C->command, func);
  if ((fp = fopen(C->file_names.iter, "w")) == NULL) {
    fprintf(stderr, "malt: Cannot write to the %s file\n", C->file_names.iter);
    exit(EXIT_FAILURE);
  }
  fprintf(fp, "Delete this file to terminate the program gracefully.\n");
  fclose(fp);
}

int checkiter(Configuration *C)
{
  FILE *fp;
  int exists = 1;

  /* check for the iterate.o file */
  if ((fp = fopen(C->file_names.iter, "r")) == NULL) {
    exists = 0;
    lprintf(C, "The %s file was deleted. Program interrupted\n", C->file_names.iter);
  } else {
    fclose(fp);
  }
  return exists;
}

#define N (C->num_params)
#define K (C->num_params_corn)

typedef struct addpoint_state {
  char *returnn;
  char *call;
  double *pc;
  pid_t pid;
  int ord;
} addpoint_t;

#define ADDPOINT_INIT                                             \
  {                                                               \
    .returnn = NULL, .call = NULL, .pc = NULL, .pid = 0, .ord = 0 \
  }

/* Finds one point on the boundary of an operating area by binary search.
 *
 * Initializes `*state` and kicks off the wrspice process.
 *
 * Returns the PID of the wrspice process so kicked off, which is suitable for passing to waitpid.
 *
 * `ord` is the ordinal corresponding to the corner being calculated by this job. It is used for
 * uniquely identifying temporary files.
 */
static pid_t start_addpoint(const Configuration *C, const Space *S, addpoint_t *state,
                            const double *pc, const double *direction, int ord)
{
#define PO_SHIFT C->options.binsearch_accuracy * 0.0001
  double *po = malloc(N * sizeof *po);  // mem:hyperplastic
  /* need to know the starting point (pc[..]) and the search direction (direction[..]) */

  /* assign values for this corner */
  state->pc = malloc((N + K) * sizeof *pc);
  memcpy(state->pc, pc, (N + K) * sizeof *pc);
  state->ord = ord;
  for (int i = N; i < N + K; ++i) {
    state->pc[i] = (ord % 2) ? S[i].cornerhi : S[i].cornerlo;
    ord /= 2;
  }

  /* find the closest boundary and calculate the search points */
  double cbig = 0.0;
  for (int i = 0; N > i; ++i) {
    double bound = ((direction[i] > 0.0) ? C->params[i].min : C->params[i].max);
    double c = direction[i] / (state->pc[i] - bound);
    if (c > cbig)
      cbig = c;
  }
  /* po on the boundary */
  /* po shifted out a bit */
  for (int i = 0; N > i; ++i) {
    po[i] = state->pc[i] - direction[i] / cbig - PO_SHIFT * direction[i];
  }
  /* calculate binsearch accuracy */
  /* in 1-D, dist=(state->pc[]-po[])/centerpnt[ ] */
  /* in N-D, calculate rms distance */
  double dist2 = 0.0;
  for (int i = 0; N > i; ++i) {
    /* binsearch_accuracy is in units of sigma */
    double dum = (state->pc[i] - po[i]);
    dist2 += dum * dum;
  }
  double dist = sqrt(dist2) / C->options.binsearch_accuracy;

  state->returnn =
      resprintf(NULL, "%s.%c.%d.return", C->command, C->function, state->ord);  // mem:kamleika
  state->call =
      resprintf(NULL, "%s.%c.%d.call", C->command, C->function, state->ord);  // mem:workmanships
  /* write .call file, call spice */
  state->pid = start_spice(C, dist, state->pc, po, state->call, state->returnn);
  free(po);  // mem:hyperplastic
  return state->pid;
#undef PO_SHIFT
}

/* Cleans up after wrspice and collects the data, storing the resulting point in pr_temp[0..N], or
 * returning 0 (without modifying pr_temp[..]) if concavity is detected.
 *
 * This function must not be called until after the wrspice process is finished, lest it read
 * incomplete data.
 *
 * `pr_temp` is the destination where the boundary point will be stored.
 *
 * On success, the value of `ord` originally passed to `start_addpoint` will be stored in `*ord`.
 */
static int addpoint_done(const Configuration *C, double *pr_temp, int *ord, addpoint_t *state)
{
  FILE *fp;

  /* read .return file */
  if ((fp = fopen(state->returnn, "r")) == NULL) {
    fprintf(stderr, "malt: Cannot open %s for reading\n", state->returnn);
    exit(EXIT_FAILURE);
  }

  int concave;
  /* read in the new point */
  int r = fscanf(fp, "%d", &concave);
  assert(1 == r);

  if (!concave) {
    /* throw away the zeroeth array element */
    int r = fscanf(fp, "%*f");
    assert(0 == r);
    for (int i = 0; N > i; ++i) {
      double pr_i;
      int r = fscanf(fp, "%lf", &pr_i);
      assert(1 == r);
      pr_temp[i] = maltspace(pr_i, C, i);
    }
  }

  // set *ord for caller
  *ord = state->ord;

  // clean up
  fclose(fp);
  unlink(state->call);
  unlink(state->returnn);
  free(state->call);  // mem:workmanships
  state->call = NULL;
  free(state->returnn);  // mem:kamleika
  state->returnn = NULL;
  free(state->pc);
  state->pc = NULL;
  state->pid = 0;

  return (!concave);
}

/* Checks if `addpoint_done` has been called on a job yet or not. */
static bool addpoint_is_done(const addpoint_t *state) { return state->pid == 0; }

/* Finds one point on the boundary of an operating area at the most limiting corner.
 *
 * This function takes the intersection of all the corners to find the smallest point in (N)d space
 * that works for all corners.
 *
 * `pc` is the center of the operating area given as a point in (N+K)d space, where N is
 * C->num_params and K is C->num_params_corn.
 *
 * `direction` is a unit vector in (N)d space representing the opposite direction along which to
 * search. (When `direction[i]` is positive the search along the (i)th axis will be in the negative
 * direction, and vice versa.)
 *
 * If `cornmin` is not NULL, `*cornmin` will be set to the ordinal value of the limiting corner.
 *
 * If `pr` is not NULL, the boundary point at the most limiting corner will be stored in `pr[0..N]`.
 *
 * Returns the distance from the origin in units of sigma.
 */
double addpoint_corners(const Configuration *C, const Space *S, corner_t *cornmin, double *pr,
                        const double *pc, const double *direction)
{
/* maximum number of subprocesses to run concurrently */
#define MAX_SUBS (C->options.max_subprocesses)
  int num_corn = 1 << K;

  /* set the corners=1 parameters to the corner values each in turn & calc margins */
  /* if there are none such then calc margins just once */
  /* *** later, the configuration can decide which corners to include if you specify *** */

  // spawn up to as many processes as there are corners, but not more than MAX_SUBS
  // (unless MAX_SUBS is zero, in which case spawn as many as you want)
  int processes = (MAX_SUBS > 0 && num_corn > MAX_SUBS) ? MAX_SUBS : num_corn;
  addpoint_t *jobs = malloc(processes * sizeof *jobs);
  int job;  // number of jobs started so far
  for (job = 0; job < processes; ++job) {
    addpoint_t init = ADDPOINT_INIT;
    jobs[job] = init;
    pid_t w = start_addpoint(C, S, &jobs[job], pc, direction, job);
    assert(w > 0);
  }

  // clean up the processes as they finish (and replace them with new ones)
  double f_min = INFINITY;                        // guaranteed to be greater than f at least once
  double *pr_temp = malloc(N * sizeof *pr_temp);  // mem:astern
  while (job < num_corn) {
    // wait for any wrspice process to finish
    int status;
    pid_t w;
    do {
      w = waitpid(-1, &status, 0);
      if (-1 == w) {
        perror("malt: waitpid");
        exit(EXIT_FAILURE);
      }
    } while (!WIFEXITED(status));

    int err = WEXITSTATUS(status);
    if (err != 0) {
      fprintf(stderr, "malt: %s returned an error (%d)\n", C->options.spice_call_name, err);
      perror("malt");
    }

    // find the job that exited by looking up its PID
    int j;
    for (j = 0; j < processes; ++j) {
      if (jobs[j].pid == w) {
        break;
      }
    }
    assert(j < processes);

    // finalize the job and check if the margin is 0
    int ord;  // ordinal of the just-finished job
    if (0 == addpoint_done(C, pr_temp, &ord, &jobs[j])) {
      /* this error message must not be here cause it jacks up the optimize routine */
      /* fprintf(stderr, "Circuit failed at nominal (%s:%d)\n", __FILE__, __LINE__); */
      // set result to 0.0
      f_min = 0.0;
      break;
    }

    // start a new job, reusing this job's slot
    w = start_addpoint(C, S, &jobs[j], pc, direction, job++);
    assert(w > 0);

    /* f is the size of the margin for this corner */
    double f2 = 0.0;
    for (int i = 0; i < N; ++i) {
      f2 += pow((pc[i] - pr_temp[i]), 2.0);
    }
    double f = sqrt(f2);

    /* pick the first corner (f_min starts at +infinity) or least corner */
    if (f_min > f) {
      f_min = f;
      if (pr != NULL) {
        for (int i = 0; i < N; ++i)
          pr[i] = pr_temp[i]; /* copy temporary result to final result */
      }
      if (cornmin != NULL) {
        *cornmin = ord;
      }
    }
  }

  // wait for the remaining jobs to finish
  for (int j = 0; j < processes; ++j) {
    // if the process has already finished: skip it
    if (addpoint_is_done(&jobs[j])) {
      continue;
    }

    // else: wait for THIS wrspice process to finish
    int status;
    do {
      if (-1 == waitpid(jobs[j].pid, &status, 0)) {
        perror("malt: waitpid");
        exit(EXIT_FAILURE);
      }
    } while (!WIFEXITED(status));

    int err = WEXITSTATUS(status);
    if (err != 0) {
      fprintf(stderr, "malt: %s returned an error (%d)\n", C->options.spice_call_name, err);
      perror("malt");
    }

    // finalize the job and check if the margin is 0
    int ord;  // ordinal of the just-finished job
    if (0 == addpoint_done(C, pr_temp, &ord, &jobs[j])) {
      // set result to 0.0
      // don't print this error for optimize, cause it is only just a convexity thing
      if (C->function != 'o')
        fprintf(stderr, "Circuit failed at nominal (%s:%d)\n", __FILE__, __LINE__);
      f_min = 0.0;
      // but keep waiting for the rest of the jobs
      continue;
    }

    /* f is the size of the margin for this corner */
    double f2 = 0.0;
    for (int i = 0; i < N; ++i) {
      f2 += pow((pc[i] - pr_temp[i]), 2.0);
    }
    double f = sqrt(f2);

    /* pick the first corner (f_min starts at +infinity) or least corner */
    if (f_min > f) {
      f_min = f;
      if (pr != NULL) {
        for (int i = 0; i < N; ++i)
          pr[i] = pr_temp[i]; /* copy temporary result to final result */
      }
      if (cornmin != NULL) {
        *cornmin = ord;
      }
    }
  }
  free(pr_temp);  // mem:astern

  /* return the distance from the origin in units of sigma, for whoever wants it */
  return f_min;
}

enum Direction {
  DOWN = 0,
  UP = 1,
};

/* Calculates the intersection of 1-D margins at all corners, storing the high margins in
 * `prhi[..]` and low margins in `prlo[..]` */
int margins(Configuration *C, const Space *S, double *prhi, double *prlo)
{
  int i, j;
  int ret = 0;
  double *pc = malloc((N + K) * sizeof *pc);          // mem:lumberer
  double *direction = malloc(N * sizeof *direction);  // mem:diallings
  double *pr = malloc(N * sizeof *pr);
  corner_t cornmin[2];

  /* Are there any included parameters? */
  if (N == 0) {
    lprintf(C, "There are no included parameters. Nothing to do here.\n");
    goto fail;
  }
  /* Header */
  lprintf(C, "\nParameter                 Nominal  Sigma    Logs? Corners?(L/H)\n");
  /* Names & Nominals */
  for (i = 0; i < N + K; i++) {
    lprintf(C, "%3d) %-19.19s %8.3f %8.4f %4d  ", i + 1, C->params[i].name,
            physspace(S[i].centerpnt, C, i), C->params[i].sigabs, C->params[i].logs ? 1 : 0);
    if (C->params[i].corners) {
      lprintf(C, "%8.3f %8.3f\n", physspace(S[i].cornerlo, C, i), physspace(S[i].cornerhi, C, i));
    } else {
      lprintf(C, "       0\n");
    }
  }
  /* 1D Margins */
  lprintf(C, "\n1D Margins           ");
  /* pad it out with whitespace */
  for (j = 0; K > j; j++) {
    lprintf(C, " ");
  }
  lprintf(C, "Min_Corner  ");
  /* pad it out with whitespace */
  for (j = 0; K > j; j++) {
    lprintf(C, " ");
  }
  lprintf(C, "Margin_in_Sigma     Margin_Parameter_Values\nParameter");
  /* pad it out with whitespace */
  for (j = 0; K > j; j++) {
    lprintf(C, " ");
  }
  lprintf(C, "            Low High    ");
  /* pad it out with whitespace */
  for (j = 0; K > j; j++) {
    lprintf(C, " ");
  }
  lprintf(C, "Low      High       Low      Nominal   High\n");
  /* initialize */
  for (i = 0; N > i; ++i) {
    for (j = 0; N > j; ++j) {
      pc[j] = S[j].centerpnt;
      direction[j] = 0.0;
    }
    /* lower margin & upper margin*/
    for (enum Direction d = DOWN; d <= UP; ++d) {
      direction[i] = (d == UP) ? -1.0 : 1.0;  // note that up is -1 and down is 1
      if (addpoint_corners(C, S, &cornmin[d], pr, pc, direction) == 0.0) {
        fprintf(stderr, "Circuit failed for nominal parameter values\n");
        goto fail;
      }
      /* used by opt & yield */
      if (d == UP) {
        prhi[i] = pr[i];
      } else {
        prlo[i] = pr[i];
      }
    }
    /* print a line */
    /* parameter name */
    lprintf(C, "%3d) %-19.19s", i + 1, C->params[i].name);
    /* the corners */
    for (j = 0; j < K; j++) {
      lprintf(C, "%s", cornmin[DOWN] & (1 << j) ? "H" : "L");
    }
    lprintf(C, " ");
    for (j = 0; j < K; j++) {
      lprintf(C, "%s", cornmin[UP] & (1 << j) ? "H" : "L");
    }
    /* the sigmas */
    lprintf(C, "     %7.2f%s %7.2f%s", (prlo[i] - S[i].centerpnt),
            (prlo[i] <= C->params[i].min) ? "*" : " ", (prhi[i] - S[i].centerpnt),
            (prhi[i] >= C->params[i].max) ? "*" : " ");
    /* the margins */
    lprintf(C, "   %8.3f%s %8.3f %8.3f%s\n", physspace(prlo[i], C, i),
            (prlo[i] <= C->params[i].min) ? "*" : " ", physspace(S[i].centerpnt, C, i),
            physspace(prhi[i], C, i), (prhi[i] >= C->params[i].max) ? "*" : " ");
  }
  ret = 1;
fail:
  free(pc);         // mem:lumberer
  free(direction);  // mem:diallings
  return ret;
}

int tmargins(Configuration *C, const Space *S)
{
  int i, j, k;
  FILE *fp;
  int ret = 0;

  /* Trace Margins */
  lprintf(C, "\nTrace Margins");
  lprintf(C, "\nParameter                Low     Nominal       High\n");
  double *pc = malloc((N + K) * sizeof *pc);          // mem:orchestrating
  double *direction = malloc(N * sizeof *direction);  // mem:rephase
  double *pr = malloc(N * sizeof *pr);
  double *prlo = malloc(N * sizeof *prlo);
  double *prhi = malloc(N * sizeof *prhi);

  /* calculate margins */
  for (i = 0; N > i; ++i) {
    for (j = 0; N > j; ++j) {
      pc[j] = S[j].centerpnt;
      direction[j] = 0.0;
    }
    /* upper margin */
    if (C->params[i].top_max) {
      strcpy(C->extensions.which_trace, C->params[i].name);
      strcat(C->extensions.which_trace, ".max");
      direction[i] = -1.0;
      if (addpoint_corners(C, S, NULL, pr, pc, direction) == 0.0) {
        fprintf(stderr, "Circuit failed for nominal parameter values\n");
        goto fail;
      }
      prhi[i] = pr[i];
    }
    /* lower margin */
    if (C->params[i].top_min) {
      strcpy(C->extensions.which_trace, C->params[i].name);
      strcat(C->extensions.which_trace, ".min");
      direction[i] = 1.0;
      if (addpoint_corners(C, S, NULL, pr, pc, direction) == 0.0) {
        fprintf(stderr, "Circuit failed for nominal parameter values\n");
        goto fail;
      }
      prlo[i] = pr[i];
    }
    /* print a line */
    if (C->params[i].top_max && C->params[i].top_min) {
      lprintf(C, "%-19.19s %8.3f%s...%8.3f ...%8.3f%s\n", C->params[i].name,
              physspace(prlo[i], C, i), (prlo[i] <= C->params[i].min) ? "*" : " ",
              physspace(S[i].centerpnt, C, i), physspace(prhi[i], C, i),
              (prhi[i] >= C->params[i].max) ? "*" : " ");
    } else if (C->params[i].top_max) {
      lprintf(C, "%-19.19s             %8.3f ...%8.3f%s\n", C->params[i].name,
              physspace(S[i].centerpnt, C, i), physspace(prhi[i], C, i),
              (prhi[i] >= C->params[i].max) ? "*" : " ");
    } else if (C->params[i].top_min) {
      lprintf(C, "%-19.19s %8.3f%s...%8.3f\n", C->params[i].name, physspace(prlo[i], C, i),
              (prlo[i] <= C->params[i].min) ? "*" : " ", physspace(S[i].centerpnt, C, i));
    }
  }

  /* make a spice control file to plot all this stuff */
  sprintf(C->file_names.plot, "%s%s.%c", C->command, C->extensions.plot, C->function);
  if (!(fp = fopen(C->file_names.plot, "w"))) {
    fprintf(stderr, "malt: Can not open the file '%s'\n", C->file_names.plot);
    exit(EXIT_FAILURE);
  }
  /* which plots to plot */
  fprintf(fp, "\n.control\n\n* which plots to plot\nnominal  = 1\nmax_min  = 1\nenvelope = 1\n\n");
  /* load the plots */
  /* the .nom file has the same base as the .envelope file, */
  fprintf(fp, "load %s.nom\nload %s\n", C->command, C->file_names.envelope);
  for (i = 0; N > i; ++i) {
    /* upper margin */
    if (C->params[i].top_max)
      fprintf(fp, "load %s.%s.max.pass\nload %s.%s.max.fail\n", C->command, C->params[i].name,
              C->command, C->params[i].name);
    /* lower margin */
    if (C->params[i].top_min)
      fprintf(fp, "load %s.%s.min.pass\nload %s.%s.min.fail\n", C->command, C->params[i].name,
              C->command, C->params[i].name);
  }
  /* nominal */
  fprintf(
      fp,
      "set group\n\nif nominal\nset title = \"param = all (nominal), node = all (define)\"\nplot");
  for (j = 0; C->num_nodes > j; ++j)
    fprintf(fp, " \\\ntran1.%s", C->nodes[j].name);
  fprintf(fp, "\nendif\n\nif max_min");
  /* max_min */
  for (k = 3, i = 0; N > i; ++i) {
    /* upper margin */
    if (C->params[i].top_max) {
      fprintf(fp, "\nset title = \"param = %s (max), node = all (pass)\"\nplot", C->params[i].name);
      for (j = 0; C->num_nodes > j; ++j)
        fprintf(fp, " \\\ntran%d.%s", k, C->nodes[j].name);
      k += 2;
    }
    /* lower margin */
    if (C->params[i].top_min) {
      fprintf(fp, "\nset title = \"param = %s (min), node = all (pass)\"\nplot", C->params[i].name);
      for (j = 0; C->num_nodes > j; ++j)
        fprintf(fp, " \\\ntran%d.%s", k, C->nodes[j].name);
      k += 2;
    }
  }
  fprintf(fp, "\nendif\n\n");
  /* envelope */
  fprintf(fp, "unset group\nset single\nif envelope\n");
  fprintf(fp, "set color2 = \"black\"\n");
  fprintf(fp, "set color3 = \"blue\"\n");
  fprintf(fp, "set color4 = \"black\"\n");
  fprintf(fp, "set color5 = \"green\"\n");
  fprintf(fp, "set color6 = \"red\"\n");
  fprintf(fp, "set color7 = \"green\"\n");
  fprintf(fp, "set color8 = \"red\"\n");
  for (k = 3, i = 0; N > i; ++i) {
    for (j = 0; C->num_nodes > j; ++j) {
      fprintf(fp, "\nset title = \"param = %s, node = %s\"\n", C->params[i].name, C->nodes[j].name);
      fprintf(fp, "plot \\\ntran2.hi%d \\\ntran1.%s \\\ntran2.lo%d", j, C->nodes[j].name, j);
      /* lower or upper margin */
      fprintf(fp, " \\\ntran%d.%s \\\ntran%d.%s", k, C->nodes[j].name, k + 1, C->nodes[j].name);
      /* both lower and upper margin */
      if (C->params[i].top_min && C->params[i].top_max) {
        fprintf(fp, " \\\ntran%d.%s \\\ntran%d.%s", k + 2, C->nodes[j].name, k + 3,
                C->nodes[j].name);
      }
    }
    k += 2;
    if (C->params[i].top_min && C->params[i].top_max)
      k += 2;
    fprintf(fp, "\n");
  }
  fprintf(fp, "endif\n\n.endc\n");
  fclose(fp);
  ret = 1;
fail:
  free(pc);         // mem:orchestrating
  free(direction);  // mem:rephase
  free(pr);
  free(prlo);
  free(prhi);
  return ret;
}

Plane **plane_malloc(Plane **plane, int *plnmemory, int inc, int dim)
{
  int i;
  Plane **new;

  /* from 0 to plnmemory-1 */
  if ((new = realloc(plane, (*plnmemory + inc) * sizeof *new)) == NULL) {  // mem:alacrity
    return NULL;
  }
  /* from 0 to plnmemory-1 */
  for (i = *plnmemory; i < (*plnmemory + inc); i++) {
    if ((new[i] = malloc(sizeof **new)) == NULL) {  // mem:pathogenicity
      return NULL;
    }
    /* from 0 to dim-1 */
    if ((new[i]->points = malloc(dim * sizeof *new[i]->points)) == NULL) {  // mem:impresting
      return NULL;
    }
    /* from 0 to dim-1 */
    if ((new[i]->a = malloc(dim * sizeof *new[i]->a)) == NULL) {  // mem:unlugubrious
      return NULL;
    }
  }
  *plnmemory += inc;
  return new;
}

void plane_free(Plane **plane, int plnmemory)
{
  int i;

  for (i = 0; i < plnmemory; i++) {
    free(plane[i]->points);  // mem:impresting
    free(plane[i]->a);       // mem:unlugubrious
    free(plane[i]);          // mem:pathogenicity
  }
  free(plane);  // mem:alacrity
}

double **margpnts_malloc(double **margpnts, int *pntmemory, int inc, int dim)
{
  int i;
  double **new;

  /* from 0 to pntmemory-1 */
  if ((new = realloc(margpnts, (*pntmemory + inc) * sizeof *new)) == NULL) {  // mem:contemporaneous
    return NULL;
  }
  /* from 0 to pntmemory-1 */
  for (i = *pntmemory; i < (*pntmemory + inc); i++) {
    /* from 0 to dim-1 */
    if ((new[i] = malloc(dim * sizeof *new[i])) == NULL) {  // mem:precisions
      return NULL;
    }
  }
  *pntmemory += inc;
  return new;
}

void margpnts_free(double **margpnts, int pntmemory)
{
  int i;

  for (i = 0; i < pntmemory; i++)
    free(margpnts[i]);  // mem:precisions
  free(margpnts);       // mem:contemporaneous
}

void intpickpnts(short *pntstack, Configuration *C, const Space *S, int depth, Plane **plane,
                 double **margpnts, int *plncount, int pntcount)
{
  for (pntstack[depth] = depth * 2; 2 + depth * 2 > pntstack[depth]; ++pntstack[depth])
    if (depth == N - 1)
      makeaplane(pntstack, C, S, plane, margpnts, plncount, pntcount);
    else
      intpickpnts(pntstack, C, S, depth + 1, plane, margpnts, plncount, pntcount);
}

int makeaplane(short *pntstack, Configuration *C, const Space *S, Plane **plane, double **margpnts,
               int *plncount, int pntcount)
{
  double distance, sumsqr, norm;
  int i, j, k, posd, negd, pntinpln;
  static int allocate = 1;
  static double b, *a, **aNmatrix;

  /* Dynamically allocate local vector and matrix, just once */
  if (allocate) {
    a = vector(0, N - 1);
    aNmatrix = matrix(1, N, 1, N);
    allocate = 0;
  }
  /* B */
  for (i = 0; N > i; ++i)
    for (j = 0; N > j; ++j)
      aNmatrix[i + 1][j + 1] = margpnts[pntstack[i]][j];
  b = det_dim(aNmatrix, N);
  /* A */
  for (k = 0; N > k; ++k) {
    for (i = 0; N > i; ++i) {
      for (j = 0; N > j; ++j)
        aNmatrix[i + 1][j + 1] = margpnts[pntstack[i]][j];
      aNmatrix[i + 1][k + 1] = -1.0;
    }
    a[k] = det_dim(aNmatrix, N);
  }
  /* check that all points (not in the plane) are on the same side of the plane */
  for (j = 0, posd = 0, negd = 0; !(posd && negd) && pntcount > j; j++) {
    for (i = 0, pntinpln = 0; !pntinpln && N > i; i++)
      pntinpln = (j == pntstack[i]);
    if (!pntinpln) {
      for (k = 0, distance = b; N > k; ++k)
        distance += a[k] * margpnts[j][k];
      if (distance > 0.0)
        posd = 1;
      else
        negd = 1;
    }
  }
  /* save the plane if it is good */
  if (!(posd && negd)) {
    /* normalize */
    for (k = 0, sumsqr = 0.0; k < N; k++)
      sumsqr += a[k] * a[k];
    norm = 1.0 / (sqrt(sumsqr));
    for (k = 0, b *= norm; k < N; k++)
      a[k] *= norm;
    for (k = 0, distance = b; k < N; k++)
      distance += a[k] * S[k].centerpnt;
    if (distance < 0.0)
      for (k = 0, b *= -1.0; k < N; k++)
        a[k] *= -1.0;
    /* store plane in global array */
    for (k = 0, plane[*plncount]->b = b; k < N; k++)
      plane[*plncount]->a[k] = a[k];
    for (k = 0, plane[*plncount]->flag = 0; k < N; k++)
      plane[*plncount]->points[k] = pntstack[k];
    ++*plncount;
  }
  return 1;
}

/* is the point inside the hull? */
/* if so, the distance to each hull plane should be positive */
int hull_dice(Configuration *C, double *pc, Plane **plane, int plncount)
{
  int pass;
  int i, k;
  double distance;

  for (i = 0, pass = 1; pass && plncount > i; i++) {
    for (k = 0, distance = plane[i]->b; N > k; ++k)
      distance += plane[i]->a[k] * pc[k];
    if (distance < 0.0)
      pass = 0;
  }
  return pass;
}

#define INSCRIBE 1e-10
int center(Configuration *C, Space *S, Plane **plane, int *tang, int plncount, double *radius)
{
  double **tab;
  int *right, *left;
  int x, y, tangent;

  /* *** is it wasteful to keep allocating over and over again? *** */
  /* allocate local vector and matrix */
  tab = matrix(1, N + 3, 1, plncount + 1 + 2 * N);  // mem:albronze
  left = ivector(1, N + 1);                         // mem:isolative
  right = ivector(1, plncount + 2 * N);             // mem:intrusional
  /* inscribe the sphere */
  /* initialize the input tabeau */
  for (y = 0; N >= y; ++y)
    tab[y + 1][1] = 0.0;
  tab[N + 2][1] = 1.0;
  for (x = 1; plncount >= x; ++x)
    tab[1][x + 1] = -plane[x - 1]->b;
  for (x = 1; N >= x; ++x) {
    tab[1][x * 2 + plncount] = -(C->params[x - 1].nom_max + INSCRIBE);
    tab[1][x * 2 + 1 + plncount] = C->params[x - 1].nom_min - INSCRIBE;
  }
  for (x = 1; plncount >= x; ++x)
    tab[N + 2][x + 1] = -1.0;
  for (x = plncount + 1; plncount + 2 * N >= x; ++x)
    tab[N + 2][x + 1] = 0.0;
  for (y = 1; N >= y; ++y)
    for (x = 1; plncount >= x; ++x)
      tab[y + 1][x + 1] = plane[x - 1]->a[y - 1];
  for (y = 1; N >= y; ++y)
    for (x = 1; N >= x; ++x) {
      tab[y + 1][x * 2 + plncount] = (y == x ? -1.0 : 0.0);
      tab[y + 1][x * 2 + 1 + plncount] = (y == x ? 1.0 : 0.0);
    }
  /* call the linear program */
  if (0 != simplx(tab, N + 1, plncount + 2 * N, right, left)) {
    tangent = 0;
    goto Bailed;
  }
  *radius = -tab[1][1];
  for (y = 1; N >= y; ++y)
    for (x = 1; plncount + 2 * N >= x; ++x)
      if (right[x] == y + plncount + 2 * N)
        S[y - 1].centerpnt = -tab[1][x + 1];
  for (x = 0, y = 1; N + 1 >= y; ++y)
    if (left[y] <= plncount)
      tang[x++] = left[y] - 1;
  tangent = x;
Bailed:
  /* free local memory */
  free_matrix(tab, 1, N + 3, 1, plncount + 1 + 2 * N);  // mem:albronze
  free_ivector(left, 1, N + 1);                         // mem:isolative
  free_ivector(right, 1, plncount + 2 * N);             // mem:intrusional
  return tangent;
}

int findface(int dim, double *facecenter, Plane **plane, int *tang, int plncount, int tangent)
{
  double **tab;
  double *sine;
  double cos, facevalue = 0.0;
  int *left, *intsect, *right;
  int x, y, z, i, large, share, noshare, count, face = -1;

  /* *** is it wasteful to keep allocating over and over again? *** */
  /* allocate local vector and matrix */
  tab = matrix(1, dim + 3, 1, dim + 3);  // mem:tactus
  left = ivector(1, dim + 1);            // mem:filets
  right = ivector(1, dim + 2);           // mem:retroreflector
  sine = vector(1, dim + 1);             // mem:females
  intsect = ivector(1, dim + 1);         // mem:massula
  /* find largest of the tangetial faces */
  for (i = 0; tangent > i; ++i) {
    /*   for (i=0; dim+1-pin > i; ++i) { */
    large = tang[i];
    if (!plane[large]->flag) {
      /* find the intersecting planes, which share N-1 points (noshare is 1) */
      for (x = 1, count = 0; plncount >= x; ++x) {
        for (y = 1, noshare = 0; noshare < 2 && dim >= y; ++y) {
          for (z = 1, share = 0; !share && dim >= z; ++z) {
            share = (plane[large]->points[z - 1] == plane[x - 1]->points[y - 1]);
          }
          if (!share)
            ++noshare;
        }
        if (noshare == 1)
          if (count++ <= dim)
            intsect[count] = x - 1;
      }
      if (count != dim) {
        /* planes are screwed up (e.g. if operating region has a planar surface), so delete it */
        /*  printf("How odd: count=%d dim=%d\n", count, dim); */
        plane[large]->flag = 1;
      } else {
        intsect[count + 1] = large;
        /* find the angles for the plane */
        for (y = 1, sine[dim + 1] = 0.0; dim >= y; ++y) {
          for (x = 1, cos = 0.0; dim >= x; ++x)
            cos += plane[intsect[y]]->a[x - 1] * plane[large]->a[x - 1];
          sine[y] = sqrt(1.0 - cos * cos);
        }
        for (y = 0; dim >= y; ++y)
          tab[y + 1][1] = 0.0;
        tab[dim + 2][1] = 1.0;
        for (x = 1; count + 1 >= x; ++x)
          tab[1][x + 1] = -plane[intsect[x]]->b;
        for (x = 1; count + 1 >= x; ++x)
          tab[dim + 2][x + 1] = -sine[x];
        for (y = 1; dim >= y; ++y)
          for (x = 1; count + 1 >= x; ++x)
            tab[y + 1][x + 1] = plane[intsect[x]]->a[y - 1];
        for (y = 0; dim + 1 >= y; ++y)
          tab[y + 1][dim + 3] = -tab[y + 1][dim + 2];
        tab[1][dim + 2] -= INSCRIBE;
        tab[1][dim + 3] -= INSCRIBE;
        /* call the linear program */
        if (0 != simplx(tab, dim + 1, count + 2, right, left)) {
          /* play with INSCRIBE, EPS1, EPS2, and plane->a (double vs. float) to fix this */
          /* Inscribe in face failed. Assign a tiny value to this facevalue */
          tab[1][1] = -INSCRIBE;
        }
        /* save face that is the biggest so far */
        if (-tab[1][1] > facevalue) {
          facevalue = (-tab[1][1]);
          face = large;
          for (y = 1; dim >= y; ++y)
            for (x = 1; dim + 2 >= x; ++x)
              if (right[x] == y + dim + 2)
                facecenter[y - 1] = -tab[1][x + 1];
        }
      }
    }
  }
  /* free local memory */
  free_matrix(tab, 1, dim + 3, 1, dim + 3);  // mem:tactus
  free_ivector(left, 1, dim + 1);            // mem:filets
  free_ivector(right, 1, dim + 2);           // mem:retroreflector
  free_vector(sine, 1, dim + 1);             // mem:females
  free_ivector(intsect, 1, dim + 1);         // mem:massula
  return face;
}

#undef K
#undef N

#undef INSCRIBE

#define TINY 1.0e-20
/* det_dim always uses the same second argument, thus the local vector can be allocated just once */
double det_dim(double **aa, int n)
{
  int i, imax = 0, j, k;
  double d, big, dum, sum, tem;
  static int allocate = 1;
  static double *vv;

  /* Dynamically allocate local vector, just once */
  if (allocate) {
    vv = vector(1, n);
    allocate = 0;
  }
  d = 1.0;
  for (i = 1; i <= n; i++) {
    big = 0.0;
    for (j = 1; j <= n; j++)
      if ((tem = fabs(aa[i][j])) > big)
        big = tem;
    if (big == 0.0) {
      nrerror("in routine det_dim");
    }
    vv[i] = 1.0 / big;
  }
  for (j = 1; j <= n; j++) {
    for (i = 1; i < j; i++) {
      sum = aa[i][j];
      for (k = 1; k < i; k++)
        sum -= aa[i][k] * aa[k][j];
      aa[i][j] = sum;
    }
    big = 0.0;
    for (i = j; i <= n; i++) {
      sum = aa[i][j];
      for (k = 1; k < j; k++)
        sum -= aa[i][k] * aa[k][j];
      aa[i][j] = sum;
      if ((dum = vv[i] * fabs(sum)) >= big) {
        big = dum;
        imax = i;
      }
    }
    if (j != imax) {
      for (k = 1; k <= n; k++) {
        dum = aa[imax][k];
        aa[imax][k] = aa[j][k];
        aa[j][k] = dum;
      }
      d = -d;
      vv[imax] = vv[j];
    }
    if (aa[j][j] == 0.0)
      aa[j][j] = TINY;
    if (j != n) {
      dum = 1.0 / (aa[j][j]);
      for (i = j + 1; i <= n; i++)
        aa[i][j] *= dum;
    }
  }
  /* tack on the determinant determination */
  for (j = 1; j <= n; j++)
    d *= aa[j][j];
  return d;
}

#undef TINY

#define EPS1 1e-9
#define EPS2 1e-9
#define FREEALL           \
  free_ivector(l2, 1, m); \
  free_ivector(l1, 1, n + 1);

static void simp1(double *const *a, int mm, int const *ll, int nll, int iabf, int *kp, double *bmax)
{
  int k;
  double test;

  *kp = ll[1];
  *bmax = a[mm + 1][*kp + 1];
  for (k = 2; k <= nll; k++) {
    if (iabf == 0)
      test = a[mm + 1][ll[k] + 1] - (*bmax);
    else
      test = fabs(a[mm + 1][ll[k] + 1]) - fabs(*bmax);
    if (test > 0.0) {
      *bmax = a[mm + 1][ll[k] + 1];
      *kp = ll[k];
    }
  }
}

static void simp2(double *const *a, int n, int const *l2, int nl2, int *ip, int kp, double *q1)
{
  int k, ii, i;
  double qp, q0, q;

  *ip = 0;
  for (i = 1; i <= nl2; i++) {
    if (a[l2[i] + 1][kp + 1] < -EPS2) {
      *q1 = -a[l2[i] + 1][1] / a[l2[i] + 1][kp + 1];
      *ip = l2[i];
      for (i = i + 1; i <= nl2; i++) {
        ii = l2[i];
        if (a[ii + 1][kp + 1] < -EPS2) {
          q = -a[ii + 1][1] / a[ii + 1][kp + 1];
          if (q < *q1) {
            *ip = ii;
            *q1 = q;
          } else if (q == *q1) {
            for (k = 1; k <= n; k++) {
              qp = -a[*ip + 1][k + 1] / a[*ip + 1][kp + 1];
              q0 = -a[ii + 1][k + 1] / a[ii + 1][kp + 1];
              if (q0 != qp)
                break;
            }
            if (q0 < qp)
              *ip = ii;
          }
        }
      }
    }
  }
}

static void simp3(double *const *a, int i1, int k1, int ip, int kp)
{
  int kk, ii;
  double piv;

  piv = 1.0 / a[ip + 1][kp + 1];
  for (ii = 1; ii <= i1 + 1; ii++)
    if (ii - 1 != ip) {
      a[ii][kp + 1] *= piv;
      for (kk = 1; kk <= k1 + 1; kk++)
        if (kk - 1 != kp)
          a[ii][kk] -= a[ip + 1][kk] * a[ii][kp + 1];
    }
  for (kk = 1; kk <= k1 + 1; kk++)
    if (kk - 1 != kp)
      a[ip + 1][kk] *= -piv;
  a[ip + 1][kp + 1] = piv;
}

int simplx(double *const *a, int m, int n, int *izrov, int *iposv)
{
  int i, ip, ir, is, k, kp, nl1, nl2;
  int *l1, *l2;
  double q1, bmax;
  /* allocate local memory */
  /* can this be done statically? */
  l1 = ivector(1, n + 1);
  l2 = ivector(1, m);
  nl1 = n;
  for (k = 1; k <= n; k++)
    l1[k] = izrov[k] = k;
  nl2 = m;
  for (i = 1; i <= m; i++) {
    if (a[i + 1][1] < 0.0)
      nrerror("in routine simplx");
    l2[i] = i;
    iposv[i] = n + i;
  }
  /* stage one */
  ir = 1;
  for (k = 1; k <= (n + 1); k++) {
    q1 = 0.0;
    for (i = 1; i <= m; i++)
      q1 += a[i + 1][k];
    a[m + 2][k] = -q1;
  }
  do {
    simp1(a, m + 1, l1, nl1, 0, &kp, &bmax);
    if (bmax <= EPS1 && a[m + 2][1] < -EPS1) {
      FREEALL
      return -1;
    } else if (bmax <= EPS1 && a[m + 2][1] <= EPS1) {
      for (ip = 1; ip <= m; ip++) {
        if (iposv[ip] == (ip + n)) {
          simp1(a, ip, l1, nl1, 1, &kp, &bmax);
          if (bmax > 0.0)
            goto one;
        }
      }
      ir = 0;
      break;
    }
    simp2(a, n, l2, nl2, &ip, kp, &q1);
    if (ip == 0) {
      FREEALL
      return -1;
    }
  one:
    simp3(a, m + 1, n, ip, kp);
    if (iposv[ip] >= (n + 1)) {
      for (k = 1; k <= nl1; k++)
        if (l1[k] == kp)
          break;
      --nl1;
      for (is = k; is <= nl1; is++)
        l1[is] = l1[is + 1];
      a[m + 2][kp + 1] += 1.0;
      for (i = 1; i <= m + 2; i++)
        a[i][kp + 1] = -a[i][kp + 1];
    }
    is = izrov[kp];
    izrov[kp] = iposv[ip];
    iposv[ip] = is;
  } while (ir);
  /* phase two */
  for (;;) {
    simp1(a, 0, l1, nl1, 0, &kp, &bmax);
    if (bmax <= 0.0) {
      FREEALL
      return 0;
    }
    simp2(a, n, l2, nl2, &ip, kp, &q1);
    if (ip == 0) {
      FREEALL
      return 1;
    }
    simp3(a, m, n, ip, kp);
    is = izrov[kp];
    izrov[kp] = iposv[ip];
    iposv[ip] = is;
  }
}

#undef EPS1
#undef EPS2
#undef FREEALL
