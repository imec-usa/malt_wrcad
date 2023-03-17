// vi: ts=2 sts=2 sw=2 et tw=100

#ifndef CONFIG
#define CONFIG

#include <stdio.h>

#define LINE_LENGTH 1024

typedef struct node {
  const char *name;
  char units;
  double dt;
  double dx;
} Node;

typedef struct param {
  /* configurable */
  const char *name;
  double nominal;
  double sig_pct;
  double sigma;
  double min;
  double max;
  int top_min;
  int top_max;
  int staticc;
  double nom_min;
  double nom_max;
  int isnommin;
  int isnommax;
  int include;
  int corners;
  int logs;
} Param;

typedef struct _2D {
  const char *name_x;
  const char *name_y;
  int param_x;
  int param_y;
} _2D;

/* these are all internal, i.e. not defined in .config */
struct file_names {
  char *circuit;
  char *param;
  char *passf;
  char *envelope;
  char *config;
  char *env_call;
  char *plot;
  char *iter;
  char *pname;
};

struct extensions {
  const char *circuit;
  const char *param;
  const char *passf;
  const char *envelope;
  const char *config;
  const char *plot;
  char *which_trace;
};

struct options {
  int verbose;
  int threads;
  int max_subprocesses;
  int print_terminal;
  double binsearch_accuracy;
  int d_simulate;
  int d_envelope;
  int o_min_iter;
  int o_max_mem_k;
  int _2D_iter;
  /* TODO: add shmoo
  float s_granularity;
  */
  int y_search_depth;
  int y_search_width;
  int y_search_steps;
  int y_max_mem_k;
  double y_accuracy;
  int y_print_every;
  const char *spice_call_name;
};

typedef struct config {
  int function;
  int func_init;
  const char *command;

  struct file_names file_names;

  struct extensions extensions;

  struct options options;

  FILE *log;

  int num_nodes;
  const Node *nodes;

  int num_params;
  int num_params_corn;
  int num_params_all;
  Param *params;

  int num_2D;
  _2D *_2D;
} Configuration;

Configuration *Configure(char *command_name, int function, char **c_text, FILE *log);
void freeConfiguration(Configuration *C);

#endif
