// vi: ts=2 sts=2 sw=2 et tw=100

#ifndef MALT
#define MALT

#include "config.h"

typedef struct args {
  int function;
  char *circuit_name;
  int verbosity;
} Args;

void strchrcat(char *, int);
char *resprintf(char **restrict strp, const char *fmt, ...);
void lprintf(const Configuration *, const char *, ...);
void malt_status(FILE *, const char *, ...);

#define info(format, ...)                                       \
  {                                                             \
    malt_status(C->log, "malt: (info) " format, ##__VA_ARGS__); \
  }
#define warn(format, ...)                                          \
  {                                                                \
    malt_status(C->log, "malt: (warning) " format, ##__VA_ARGS__); \
  }
#define error(format, ...)                                       \
  {                                                              \
    malt_status(C->log, "malt: (error) " format, ##__VA_ARGS__); \
    exit(EXIT_FAILURE);                                          \
  }

#define MALTUSAGE "malt [-h] {-d|-m|-t|-2|-y|-o} <circuit_name>[.<ext_1>[.<ext_2>[...]]]\n"
#define MALTVERSION "3.0"
#define MALTHELP                                                            \
  "MALT " MALTVERSION "\n"                                                  \
  "  Parametric yield optimization utility for use with WRSpice\n"          \
  "SYNOPSIS\n"                                                              \
  "  " MALTUSAGE "OPTIONS\n"                                                \
  "  -h\tPrint help message and exit\n"                                     \
  "  -d\tDefine correct circuit operation\n"                                \
  "  -m\tCalculate individual parameter margins\n"                          \
  "  -t\tTrace nodes at marginal parameter values\n"                        \
  "  -2\tCalculate operating region in 2 dimensions\n"                      \
  "  -y\tCalculate circuit yield using corner analysis\n"                   \
  "  -o\tOptimize yield using inscribed hyperspheres\n"                     \
  "CIRCUIT NAME\n"                                                          \
  "  Used to find the .config, .cir, .param, .passf, and .envelope files\n" \
  "  Name extensions are delimited by periods\n"                            \
  "  The most specfic .cir, .param, .passf, and .envelope files are used\n" \
  "  The .config files are parsed from most general to most specific\n"

#endif
