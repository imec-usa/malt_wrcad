// vi: ts=2 sts=2 sw=2 et tw=100
/* parse the command line */
#include "malt.h"
#include "config.h"
#include "corners.h"
#include "define.h"
#include "margins.h"
#include "optimize.h"
#include <assert.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

void read_command_line(Args *args, int argc, char *argv[]);

int main(int argc, char *argv[])
{
  Args args = {0};

  /* parse command line */
  read_command_line(&args, argc, argv);

  /* create temporary log file */
  FILE *log = tmpfile();
  if (log == NULL) {
    char *err = resprintf(NULL, "malt: Can't open a temporary log file");
    perror(err);
    exit(EXIT_FAILURE);
  }

  /* find and parse the cascade of .toml configuration files */
  Configuration *C = Configure(&args, log);

  // TODO: resolve filenames properly, instead of changing directory here, which might be confusing
  chdir(lst_last(&C->working_tree));

  /* .cir */
  if (C->file_names.circuit != NULL) {
    info("SPICE (%s) file: '%s'\n", C->extensions.circuit, C->file_names.circuit);
  } else {
    error("The %s file is the SPICE deck and is required\n", C->extensions.circuit);
  }

  /* parameters */
  if (C->file_names.param != NULL) {
    info("Parameters file: '%s'\n", C->file_names.param);
  }

  /* passfail or envelope */
  if (C->file_names.passf != NULL) {
    info("Pass/fail file: '%s'\n", C->file_names.passf);
  } else if (C->file_names.envelope != NULL) {
    info("Envelope file: '%s'\n", C->file_names.envelope);
  } else if (args.function != 'd') {
    // if neither the .envelope or .passf file exists, only -d makes sense
    error("Cannot find an envelope or passfail file. Try running -d first?\n");
  }

  /* call the appropriate algorithm */
  switch (args.function) {
  case 'd':
    if (!call_def(C))
      fprintf(stderr, "Define routine exited on an error\n");
    break;
  case 'm':
    if (!call_marg(C))
      fprintf(stderr, "Margins routine exited on an error\n");
    break;
  case 't':
    if (!call_trace(C))
      fprintf(stderr, "Trace routine exited on an error\n");
    break;
  case '2':
    if (!margins2(C))
      fprintf(stderr, "2D Margins routine exited on an error\n");
    break;
  case 's':
    if (!shmoo(C))
      fprintf(stderr, "2D Shmoo routine exited on an error\n");
    break;
  case 'y':
    if (!marg_corners(C))
      fprintf(stderr, "Yield routine exited on an error\n");
    break;
  case 'o':
    if (!call_opt(C))
      fprintf(stderr, "Optimize routine exited on an error\n");
    break;
  }
  free(args.configuration);  // mem:mobster
  freeConfiguration(C);
}

__attribute__((noreturn)) static void usage(void)
{
  fprintf(stderr, "Usage: %s", MALTUSAGE);
  exit(EXIT_FAILURE);
}

void read_command_line(Args *args, int argc, char *argv[])
{
  int have_function = 0;
  int c;

  /* find options and arguments */
  while ((c = getopt(argc, argv, "hdmt2syovk")) != -1) {
    /* have an option */
    switch (c) {
    case 'h':
      fprintf(stderr, "%s", MALTHELP);
      exit(EXIT_SUCCESS);
    case 'd':
    case 'm':
    case 't':
    case '2':
    case 's':
    case 'y':
    case 'o':
      if (have_function) {
        /* only one function is allowed per command line */
        usage();
      }
      args->function = c;
      have_function = 1;
      break;
    case 'v':
      args->verbosity += 1;
      break;
    case 'k':
      args->keep_files = true;
      break;
    case '?':
      usage();
    default:
      /* should not be possible */
      abort();
    }
  }

  /* one of the function options must be provided */
  if (!have_function) {
    usage();
  }

  if (optind + 1 == argc) {
    /* we have a circuit name */
    args->configuration = strdup(argv[optind]);  // mem:mobster
  } else {
    usage();
  }
}

/* Print formatted output both to stdout (when C->print_terminal is true), but also to the
 * configured log file, C->log. */
void lprintf(const Configuration *C, const char *fmt, ...)
{
  va_list args;
  // print to terminal
  if (C->options.print_terminal) {
    va_start(args, fmt);
    vprintf(fmt, args);
    va_end(args);
  }

  // print to log file
  va_start(args, fmt);
  vfprintf(C->log, fmt, args);
  va_end(args);

  // flush logfile
  fflush(C->log);
}

/* Print formatted status message both to stderr and to log. */
void malt_status(FILE *log, const char *format, ...)
{
  va_list args;
  // print to terminal
  va_start(args, format);
  vfprintf(stderr, format, args);
  va_end(args);

  // print to log file
  va_start(args, format);
  vfprintf(log, format, args);
  va_end(args);

  // flush logfile
  fflush(log);
}

/* Works like sprintf, but (re)allocates beforehand
 * Implementation adapted from:
 * https://stackoverflow.com/questions/4899221/substitute-or-workaround-for-asprintf-on-aix/4899487#4899487
 */
char *resprintf(char **restrict strp, const char *fmt, ...)
{
  va_list args;

  // calculate required size of buffer
  va_start(args, fmt);
  int reserved = vsnprintf(NULL, 0, fmt, args);
  va_end(args);
  assert(reserved >= 0);

  char *ret = NULL;
  /* if the passed pointer is NULL, allocate anyway using ret as a temporary pointee */
  if (strp == NULL) {
    strp = &ret;
  } else {
    ret = *strp;
  }

  // allocate space
  *strp = realloc(*strp, reserved + 1);
  assert(*strp != NULL);

  // populate buffer
  va_start(args, fmt);
  int copied = vsnprintf(*strp, reserved + 1, fmt, args);
  va_end(args);
  assert(copied == reserved);

  return ret;
}
