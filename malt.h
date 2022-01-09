
#ifndef MALT
#define MALT

//typedef struct config Configuration;
#include "config.h"

void strchrcat(char *, int);
char *resprintf(char **restrict strp, const char *fmt, ...);
void lprintf(const Configuration *, const char *, ...);

#define MALTUSAGE "malt [-h] {-d|-m|-t|-2|-y|-o} <circuit_name>[.<ext_1>[.<ext_2>[...]]]\n"
#define MALTHELP \
"MALT 3.0\n" \
"  Parametric yield optimization utility for use with WRSpice\n" \
"SYNOPSIS\n" \
"  " MALTUSAGE \
"OPTIONS\n" \
"  -h\tPrint help message and exit\n" \
"  -d\tDefine correct circuit operation\n" \
"  -m\tCalculate individual parameter margins\n" \
"  -t\tTrace nodes at marginal parameter values\n" \
"  -2\tCalculate operating region in 2 dimensions\n" \
"  -y\tCalculate circuit yield using corner analysis\n" \
"  -o\tOptimize yield using inscribed hyperspheres\n" \
"CIRCUIT NAME\n" \
"  Used to find the .config, .cir, .param, .passf, and .envelope files\n" \
"  Name extensions are delimited by periods\n" \
"  The most specfic .cir, .param, .passf, and .envelope files are used\n" \
"  The .config files are parsed from most general to most specific\n"

#endif
