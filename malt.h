// vi: ts=2 sts=2 sw=2 et tw=100

#ifndef MALT
#define MALT

#include "config.h"

typedef struct args {
  int function;
  bool keep_files;
  char *configuration;
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

#define MALTUSAGE "malt [-h] {-d|-m|-t|-2|-y|-o} [-k] CONFIG\n"
#define MALTVERSION "3.1.1"
#define MALTHELP                                                                      \
  "Malt " MALTVERSION "\n"                                                            \
  "  Parametric yield optimization utility for use with WRSpice\n"                    \
  "USAGE\n"                                                                           \
  "  " MALTUSAGE "\n"                                                                 \
  "COMMANDS\n"                                                                        \
  "  -h\tPrint help message and exit\n"                                               \
  "  -d\tDefine correct circuit operation\n"                                          \
  "  -m\tCalculate individual parameter margins\n"                                    \
  "  -t\tTrace nodes at marginal parameter values\n"                                  \
  "  -2\tCalculate operating region in 2 dimensions\n"                                \
  "  -y\tCalculate circuit yield using corner analysis\n"                             \
  "  -o\tOptimize yield using inscribed hyperspheres\n"                               \
  "\n"                                                                                \
  "OPTIONS\n"                                                                         \
  "  -k\tKeep (don't delete) additional temporary files\n"                            \
  "\n"                                                                                \
  "CONFIG\n"                                                                          \
  "  Path (relative to the current directory) of the most specific applicable\n"      \
  "  configuration (.toml) file for this analysis. When Malt runs, it searches in\n"  \
  "  the current directory and all ancestors for a file named Malt.toml, which is\n"  \
  "  taken as the base configuration. Then all .toml files of which CONFIG is a\n"    \
  "  path prefix are processed from least to most specific.\n"                        \
  "\n"                                                                                \
  "  Generated files such as parameter, pass/fail, and envelope files are searched\n" \
  "  for in similar fashion starting from the _malt directory, if one is present.\n"

#endif
