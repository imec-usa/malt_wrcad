// vi: ts=2 sts=2 sw=2 et tw=100

#ifndef CONFIG
#define CONFIG

#include "list.h"
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
  char *env_call;
  char *data;
  char *plot;
  char *gnuplot;
  // TODO: delete these
  char *iter;
  char *pname;
};

struct extensions {
  const char *circuit;
  const char *param;
  const char *passf;
  const char *envelope;
  const char *env_call;
  const char *plot;
  const char *pname;
  char *which_trace;
};

enum filetype {
  Ft_Circuit,
  Ft_Parameters,
  Ft_Pname,
  Ft_PassFail,
  Ft_Envelope,
  Ft_EnvCall,
  Ft_Data,
  Ft_Plot,
  Ft_Gnuplot,
};

struct options {
  int spice_verbose;
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
  const char *command;
  FILE *log;

  int func_init;

  bool keep_files;

  list_t working_tree;

  struct file_names file_names;

  struct extensions extensions;

  struct options options;

  int num_nodes;
  const Node *nodes;

  int num_params;
  int num_params_corn;
  int num_params_all;
  Param *params;

  int num_2D;
  _2D *_2D;
} Configuration;

typedef struct args Args;

Configuration *Configure(const Args *args, FILE *log);
void freeConfiguration(Configuration *C);

const char *malt_filename(Configuration *C, enum filetype kind);
FILE *new_file_by_type(Configuration *C, enum filetype kind);

void unlink_pname(Configuration *C);

#endif
