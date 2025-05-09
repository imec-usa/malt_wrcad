// vi: ts=2 sts=2 sw=2 et tw=100
#include "call_spice.h"
#include "config.h"
#include "malt.h"
#include "space.h"
#include <assert.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/wait.h>
#include <unistd.h>

/* Creates SPICE input files that are used by other routines.
 *
 * Returns 0 if opening any of the files fails.
 */
static int generic_spice_files(const Configuration *C)
{
  // these files go in the root of the working tree
  const char *wd = C->working_tree.ptr[0];
  char *filename = NULL;
  /* create generic files even if they already exist */
#define CREATE_FILE(name, contents, ...)                      \
  {                                                           \
    resprintf(&filename, "%s/" name, C->working_tree.ptr[0]); \
    FILE *fp = fopen(filename, "w");                          \
    if (fp == NULL) {                                         \
      error("Cannot create %s\n", filename);                  \
    }                                                         \
    fprintf(fp, contents, ##__VA_ARGS__);                     \
    fclose(fp);                                               \
  }

  CREATE_FILE(MALT_RUN_FILENAME, MALT_RUN)
  CREATE_FILE(MALT_BINSEARCH_FILENAME, MALT_BINSEARCH, wd)
  CREATE_FILE(MALT_PASSFAIL_FILENAME, MALT_PASSFAIL)

#undef CREATE_FILE
  free(filename);

  return 1;
}

/* remake pname file whenever included/excluded parameters change */
void pname(Configuration *C)
{
  FILE *fp = new_file_by_type(C, Ft_Pname);
  fprintf(fp, "* %s\n\n.control\n\n", malt_filename(C, Ft_Pname));
  for (int i = 0; C->num_params_all > i; ++i) {
    fprintf(fp, "%s = param[%i]\n", C->params[i].name, i + 1);
  }
  fprintf(fp, "\n.endc\n");
  fclose(fp);
}

/* Writes the input (.call) file and calls SPICE to perform a binary search.
 *
 * Returns the PID of the spawned SPICE process.
 *
 * `accuracy` is the tolerance of the binsearch algorithm
 * `pc` is the center point of the search (inner edge)
 * `po` is the outer edge of the search
 * `call` is the name of the input file to SPICE
 * `returnn` is the name of the output file.
 */
pid_t start_spice(const Configuration *C, double accuracy, double *pc, double *po, const char *call,
                  const char *returnn)
{
  FILE *fp;
  int i;
  static int dasht = 1;

  /* generic file generation */
  if (C->function == 'd') {
    generic_spice_files(C);
  }
  /* malt2spice file */
  if ((fp = fopen(call, "w")) == NULL) {
    fprintf(stderr, "malt: Cannot write to the '%s' file", call);
    exit(EXIT_FAILURE);
  }
  /* header stuff */
  fprintf(fp, "* %s\n\n.control\n\n", call);
  fprintf(fp, "set circuit = ( %s )\n", C->file_names.circuit);
  fprintf(fp, "set param   = ( %s )\n", (C->file_names.param) ? C->file_names.param : "no file");
  fprintf(fp, "set passf   = ( %s )\n", (C->file_names.passf) ? C->file_names.passf : "no file");
  fprintf(fp, "set pname   = ( %s )\n", malt_filename((Configuration*)C, Ft_Pname));
  fprintf(fp, "set return  = ( %s )\n", returnn);
  /* node math is legal (i.e. v(1)-v(2)) */
  /* ...so long as the component vectors also appear individually */
  fprintf(fp, "set node_name = (");
  for (i = 0; C->num_nodes > i; ++i)
    fprintf(fp, " %s ", C->nodes[i].name);
  fprintf(fp, ")\n");
  /* * * yield routine (or not) * * */
  /* skip the nominal sim for efficiency, after the initial margins */
  if (C->function == 'y' && C->func_init == 0)
    fprintf(fp, "dashc = 1\n");
  else
    fprintf(fp, "dashc = 0\n");
  /* * * trace routine (or not) * * */
  if (C->function == 't') {
    fprintf(fp, "dasht = %d\n", dasht);
    if (dasht == 1)
      dasht = 2;
    // TODO: find out if this is a data hazard for multiprocessing ~ntj
    const char *wd = lst_last((list_t*)&C->working_tree);
    fprintf(fp, "set n_return  = ( %s/nominals )\n", wd);
    fprintf(fp, "set p_return  = ( %s/%s.pass )\n", wd, C->extensions.which_trace);
    fprintf(fp, "set f_return  = ( %s/%s.fail )\n", wd, C->extensions.which_trace);
  } else if (C->function != 'd')
    fprintf(fp, "dasht = 0\n");
  /* * * define routine (or not) * * */
  if (C->function == 'd') {
    fprintf(fp, "set n_return  = ( %s/nominals )\n", lst_last((list_t*)&C->working_tree));
    /* nominal parameter values */
    for (i = 0; C->num_params_all > i; ++i)
      fprintf(fp, "param[%i]=%g\n", i + 1, C->params[i].nominal);
    fprintf(fp, "\nsource %s/%s\n\n.endc\n", C->working_tree.ptr[0], MALT_RUN_FILENAME);
  } else {
    fprintf(fp, "set envelope = ( %s )\n", C->file_names.env_call);
    /* binary search limits stuff */
    fprintf(fp, "pc[0]=0\npo[0]=%g\n", accuracy);
    /* pc on plane, po on the boundary */
    for (i = 0; C->num_params > i; ++i) {
      fprintf(fp, "pc[%i]=%g\n", i + 1, physspace(pc[i], C, i));
      fprintf(fp, "po[%i]=%g\n", i + 1, physspace(po[i], C, i));
      fprintf(fp, "pl[%i]=%d\n", i + 1, C->params[i].logs);
    }
    /* tack on the corner parameters */
    for (; (C->num_params + C->num_params_corn) > i; ++i) {
      fprintf(fp, "pc[%i]=%g\n", i + 1, physspace(pc[i], C, i));
      fprintf(fp, "po[%i]=%g\n", i + 1, physspace(pc[i], C, i));
      fprintf(fp, "pl[%i]=%d\n", i + 1, C->params[i].logs);
    }
    /* tack on the excluded parameters */
    for (; C->num_params_all > i; ++i) {
      fprintf(fp, "pc[%i]=%g\n", i + 1, C->params[i].nominal);
      fprintf(fp, "po[%i]=%g\n", i + 1, C->params[i].nominal);
      fprintf(fp, "pl[%i]=%d\n", i + 1, C->params[i].logs);
    }
    fprintf(fp, "\nsource %s/%s\n\n.endc\n", C->working_tree.ptr[0], MALT_BINSEARCH_FILENAME);
  }
  /* all routines */
  if (fclose(fp)) {
    fprintf(stderr, "malt: Error while writing to %s\n", call);
    perror("malt");
  }
  // empty buffers so they don't get flushed in the child as well
  fflush(stdout);
  pid_t wrspice = fork();
  if (0 == wrspice) {
    // child: call spice
    freopen("/dev/null", "w", stdin);  // hack to force batch mode is hacky
    if (C->options.spice_verbose) {
      freopen(".verbage", "w", stdout);
    } else {
      freopen("/dev/null", "w", stdout);
    }
    execlp(C->options.spice_call_name, "wrspice", "-b", call, (char *)NULL);
    // if execlp returns, it failed
    perror("malt: execlp");
    fprintf(stderr, "malt: (called %s)\n", C->options.spice_call_name);
    exit(EXIT_FAILURE);
  }
  // parent: return PID or 0 if unsuccessful
  if (wrspice == -1) {
    perror("malt: fork");
  }
  return wrspice;
}

/* Calls start_spice and waits for the wrspice process to finish before returning.
 *
 * define calls this function.
 */
void call_spice(const Configuration *C, double accuracy, double *pc, double *po, const char *call,
                const char *returnn)
{
  pid_t wrspice = start_spice(C, accuracy, pc, po, call, returnn);
  assert(wrspice > 0);

  int status;
  do {
    if (-1 == waitpid(wrspice, &status, 0)) {
      perror("malt: waitpid");
      exit(EXIT_FAILURE);
    }
  } while (!WIFEXITED(status));  // TODO: this is probably bogus ~ntj

  int err = WEXITSTATUS(status);
  if (err != 0) {
    fprintf(stderr, "malt: %s returned an error (%d)\n", C->options.spice_call_name, err);
    perror("malt");
  }
}
