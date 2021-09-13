/* parse configuraton files */

/* *** need to have a flag for each composit line field */
/* which says if the field is defined or not */
/* if not, do not print it in dumpConfiguration */
/* source final configuration file at the very end */
/* which trumps everything previous */
/* need a units field for dt and dx */
#include "config.h"
#include "malt.h"
#include <assert.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h> /* for strcasecmp */
#include <unistd.h>

typedef struct config_line{
  char **key;
  char **value;
  int  length;
} ConfigLine;

typedef struct builder {
  int function;
  char *command;

  struct file_names file_names;

  struct extensions extensions;

  struct options options;

  FILE *log;

  Node node_defaults;
  Param param_defaults;

  int  num_nodes;
  Node *nodes;

  int   num_params_all;
  Param *params;

  int num_2D;
  _2D *_2D;
} Builder;

void c_term_file_write(char **c_text, const char *text);

static void builder_debug(const Builder *B,FILE *fp){
  int i;
  /* organize according to function */
  fprintf(fp,"***  Configuration File  ***\n");
  if(B->file_names.circuit) fprintf(fp,"*circuit file name: %s\n",B->file_names.circuit);
  if(B->file_names.param) fprintf(fp,"*parameter file name: %s\n",B->file_names.param);
  if(B->file_names.passf) fprintf(fp,"*passfail file name: %s\n",B->file_names.passf);
  if(B->file_names.envelope) fprintf(fp,"*envelope file name: %s\n",B->file_names.envelope);
  fprintf(fp,"\n***  Circuit Node Defaults  ***\n");
  fprintf(fp,"dt = %g\n",B->node_defaults.dt);
  fprintf(fp,"dx = %g\n",B->node_defaults.dx);
  fprintf(fp,"\n***  Circuit Parameter Defaults  ***\n");
  fprintf(fp,"nominal = %g\n",B->param_defaults.nominal);
  fprintf(fp,"min = %g\n",B->param_defaults.min);
  fprintf(fp,"max = %g\n",B->param_defaults.max);
  /* by default, both sigma and sig_abs are zero. but one-and-only-one must be non-zero in the end */
  fprintf(fp,"sigma = %g\n",B->param_defaults.sigma);
  fprintf(fp,"sig_abs = %g\n",B->param_defaults.sigabs);
  fprintf(fp,"static = %d\n",B->param_defaults.staticc);
  fprintf(fp,"logs = %d\n",B->param_defaults.logs);
  fprintf(fp,"corners = %d\n",B->param_defaults.corners);
  fprintf(fp,"include = %d\n",B->param_defaults.include);
  fprintf(fp,"\n***  Circuit Nodes  ***\n");
  for(i=0; i<B->num_nodes; ++i)
    fprintf(fp,"node = %s, dt = %g, dx = %g\n",B->nodes[i].name,B->nodes[i].dt,B->nodes[i].dx);
  fprintf(fp,"\n***  Circuit Parameters  ***\n");
  for(i=0; i<B->num_params_all; ++i){
    fprintf(fp,"param = %s, nominal = %g, min = %g, max = %g, sigma = %g, "
               "sig_abs = %g, static = %d, logs = %d, corners = %d, include = %d",
            B->params[i].name,   B->params[i].nominal,
            B->params[i].min,    B->params[i].max,
            B->params[i].sigma,  B->params[i].sigabs,
            B->params[i].staticc,B->params[i].logs,
            B->params[i].corners,B->params[i].include);
    if(B->params[i].isnommin)
      fprintf(fp,", nom_min = %g", B->params[i].nom_min);
    if(B->params[i].isnommax)
      fprintf(fp,", nom_max = %g", B->params[i].nom_max);
    fprintf(fp,"\n");}
  fprintf(fp,"\n***  2D Margins  ***\n");
  for(i=0; i<B->num_2D; ++i)
    fprintf(fp,"param_x = %s, param_y = %s\n", B->_2D[i].name_x, B->_2D[i].name_y);
  /* fake line as an example */
  if(B->num_2D == 0) {
    fprintf(fp,"***  Example. Replace arguments with actual param names  ***\n");
    fprintf(fp,"* param_x = %s, param_y = %s\n", "XJ", "XL");}
  fprintf(fp,"\n***  File Extensions  ***\n");
  if(B->extensions.circuit) fprintf(fp,"circuit_extension = %s\n",B->extensions.circuit);
  if(B->extensions.param) fprintf(fp,"parameters_extension = %s\n",B->extensions.param);
  if(B->extensions.passf) fprintf(fp,"passfail_extension = %s\n",B->extensions.passf);
  if(B->extensions.envelope) fprintf(fp,"envelope_extension = %s\n",B->extensions.envelope);
  if(B->extensions.plot) fprintf(fp,"plot_extension = %s\n",B->extensions.plot);
  fprintf(fp,"\n***  General Options  ***\n");
  fprintf(fp,"*fraction of sigma:\nbinsearch_accuracy = %g\n",B->options.binsearch_accuracy);
  if(B->options.spice_call_name) fprintf(fp,"spice_call_name = %s\n",B->options.spice_call_name);
  fprintf(fp,"max_subprocesses = %d\n",B->options.max_subprocesses);
  fprintf(fp,"threads = %d\n",B->options.threads);
  fprintf(fp,"verbose = %d\n",B->options.verbose);
  fprintf(fp,"print_terminal = %d\n",B->options.print_terminal);
  fprintf(fp,"\n***  Options for Define  ***\n");
  fprintf(fp,"d_simulate = %d\n",B->options.d_simulate);
  fprintf(fp,"d_envelope = %d\n",B->options.d_envelope);
  fprintf(fp,"\n***  Options for Margins  ***\n");
  fprintf(fp,"\n***  Options for 2D Margins  ***\n");
  fprintf(fp,"2D_iter = %d\n",B->options._2D_iter);
  fprintf(fp,"\n***  Options for Optimize  ***\n");
  fprintf(fp,"o_min_iter = %d\n",B->options.o_min_iter);
  fprintf(fp,"o_max_mem_k = %d\n",B->options.o_max_mem_k);
}

static void _2D_drop(_2D *ptr) {
  free((void *)ptr->name_x); // schizzo-carabineros-outvalued
  ptr->name_x = NULL;
  free((void *)ptr->name_y); // foreconsent-hyperbolic-theophilist
  ptr->name_y = NULL;
}

static void param_drop(Param *ptr) {
  free((void *)ptr->name); // reminiscer-caramboled-disburdening
  ptr->name = NULL;
}

static void node_drop(Node *nptr) {
  free((void *)nptr->name); // sunglow-gunline-eraser
  nptr->name = NULL;
}

void freeConfiguration(Configuration *C)
{
  free(C->file_names.circuit); // tectospondylous-erodability-empyreuma
  free(C->file_names.param); // stairwork-heliophotography-licente
  free(C->file_names.passf); // toners-cephalomotor-comburivorous
  free(C->file_names.envelope); // cool-mesorectums-vicugnas
  free(C->file_names.config); // synurae-handmaids-myelosuppressions
  free(C->file_names.env_call); // semishady-unloath-tragicomical
  free(C->file_names.plot); // federalizes-spawny-trooper
  free(C->file_names.iter); // myrcene-unallusively-cuerda
  free(C->file_names.pname); // physnomy-soppy-hypothyreosis
  free(C->extensions.which_trace); // intratubal-upchaunce-underage
  free((void *)C->options.spice_call_name); // descendentalism-crenae-rationalisation

  /* FIXME: add this part when we can guarantee log has been opened ~ntj
  fclose(C->log);
  */

  for (int n = 0; n < C->num_nodes; ++n) {
    node_drop((Node *)&C->nodes[n]);
  }
  free((void *)C->nodes); // anathematise-intenancy-sublimating

  for (int p = 0; p < C->num_params_all; ++p) {
    param_drop((Param *)&C->params[p]);
  }
  free((void *)C->params); // restringer-unlustier-stairhead

  for (int t = 0; t < C->num_2D; ++t) {
    _2D_drop(&C->_2D[t]);
  }
  free(C->_2D); // hypersexual-junkerish-flashguns

  free(C); // circumgyrate-semicylinder-keypunching
}

void builder_init(Builder *C){
  /* allocate everything */

  C->log        = NULL;
  /* initialize everything to internal defaults */
  /* iternal file names */
  C->file_names.circuit=NULL;
  C->file_names.param=NULL;
  C->file_names.passf=NULL;
  C->file_names.envelope=NULL;
  C->file_names.config=NULL;
  C->file_names.env_call=NULL;
  C->file_names.plot=NULL;
  C->file_names.iter=NULL;
  C->file_names.pname=NULL;
  /* in config files, defaults must be defined before nodes and params which use defaults */
  /* node defaults */
  C->node_defaults.name=NULL;
  C->node_defaults.units='V';
  C->node_defaults.dt=100e-12;
  C->node_defaults.dx=1;
  /* parameter defaults */
  C->param_defaults.name=NULL;
  C->param_defaults.nominal=1;
  C->param_defaults.sigma=0.0;
  C->param_defaults.sigabs=0.0;
  C->param_defaults.min=0.5;
  C->param_defaults.max=2;
  C->param_defaults.top_min=0;
  C->param_defaults.top_max=0;
  C->param_defaults.staticc=0;
  C->param_defaults.isnommin=0;
  C->param_defaults.isnommax=0;
  C->param_defaults.include=1;
  C->param_defaults.logs=1;
  C->param_defaults.corners=0;
  /* nodes, params, and 2D */
  C->nodes = NULL;
  C->params = NULL;
  C->_2D = NULL;
  C->num_nodes=0;
  C->num_params_all=0;
  C->num_2D=0;
  /* extensions */
  C->extensions.circuit  = ".cir";
  C->extensions.param    = ".param";
  C->extensions.passf    = ".passf";
  C->extensions.envelope = ".envelope";
  C->extensions.config   = ".config";
  C->extensions.plot     = ".plot";
  /* pretty weak */
  C->extensions.which_trace = malloc(LINE_LENGTH); // intratubal-upchaunce-underage
  /* options */
  C->options.binsearch_accuracy = 0.1;
  C->options.spice_call_name = strdup("wrspice"); // descendentalism-crenae-rationalisation
  C->options.verbose = 0;
  C->options.threads = 16; // default # threads: 16
  C->options.max_subprocesses = 0; // default # jobs: unlimited
  C->options.print_terminal = 1;
  /* options for define */
  C->options.d_simulate = 1;
  C->options.d_envelope = 1;
  /* options for margins */
  /* options for 2D margins */
  C->options._2D_iter = 16;
  /* options for optimize */
  C->options.o_min_iter = 100;
  C->options.o_max_mem_k = 4194304;
}

/* Removes leading and trailing whitespace from a string. */
void remove_white_space(char *token){
  char *begin = token;
  while (*begin && isspace(*begin)) {
    ++begin;
  }
  size_t len = 0;
  while (begin[len]) {
    ++len;
  }
  while (len > 0 && isspace(begin[len-1])) {
    --len;
  }
  memmove(token, begin, len);
  token[len] = '\0';
}

int numDelims(char *str, char c){
  char *cptr;
  if(*str=='\0'){
    fprintf(stderr,"malt: Parse error in config file\n");
    exit(EXIT_FAILURE);}
  if((cptr=strchr(str,c))==NULL) {
    return(0);
  } else
    return(1 + numDelims(cptr+1,c));
}

int numPairs(char *line){
  int equal_signs, commas;
  equal_signs = numDelims(line,'=');
  commas = numDelims(line,',');
  if(equal_signs != (commas+1)){
    fprintf(stderr,"malt: Parse error in config file\n%s\n",line);
    exit(EXIT_FAILURE);}
  return(equal_signs);
}

int interpretLine(ConfigLine *L, char *buf, int index){
  char *endc,*kptr,*kptr1;

  /*Is this the last key/value entry?*/
  if((endc = (char *)strchr(buf,','))!=NULL) *endc='\0';

  /* chop at first equal sign */
  kptr=buf;  kptr1 = (char *)strchr(kptr,'='); *kptr1='\0';
  remove_white_space(kptr);
  L->key[index] = strdup(kptr); // smutchier-purify-holoku

  /* step beyond the equal sign */
  kptr = kptr1+1;
  remove_white_space(kptr);
  L->value[index] = strdup(kptr); // hypersentimental-alstonine-nonword

  /* There are additional pairs to process */
  if(endc!=NULL)
    return(1 + interpretLine(L, endc+1, index+1));
  else  /* There are NOT additional pairs to process */
    return(1);
}

ConfigLine *parseConfigLine(ConfigLine *L, char *buf){
  if(buf==NULL || *buf == '\0'){
    fprintf(stderr,"malt: Parse error in routine 'parseConfigLine'\n");
    return(NULL);}
  if((L->length = numPairs(buf))==0) {
    fprintf(stderr,"malt: Parse error in config file\n");
    exit(EXIT_FAILURE);}
  if ((L->key = calloc(L->length+1, sizeof *L->key)) == NULL) { // kiekie-thersitean-newfangled
    fprintf(stderr,"malt: Memory allocation error in routine 'parseConfigLine'\n");
    exit(EXIT_FAILURE);
  }
  if ((L->value = calloc(L->length+1, sizeof *L->value)) == NULL) { // semisoun-yirring-hyperperfection
    fprintf(stderr,"malt: Memory allocation error in routine 'parseConfigLine'\n");
    exit(EXIT_FAILURE);
  }
  if(!(interpretLine(L,buf,0))) {
    fprintf(stderr,"malt: Parse error in config file\n");
    exit(EXIT_FAILURE);
  }
  return(L);
}

/*Returns 0 on failure*/
/*Returns number of keywords recognized on success*/
int compositLine(Builder *C, ConfigLine *L){
  int i,j,match,this_num;

  /* node/param keyword need not be at beginning of line, but grab the first one you come to */
  for(i=0; i < L->length; ++i){
    if(!strcasecmp(L->key[i],"node")){
      /* does this node already exist? */
      for(j=0,match=0; (j < C->num_nodes) && !match; ++j)
        match=!strcmp(C->nodes[j].name,L->value[i]);
      this_num=j;
      /* allocate memory for a new node */
      if (!match) {
        if ((C->nodes = realloc(C->nodes, (++C->num_nodes)*sizeof *C->nodes)) == NULL) { // anathematise-intenancy-sublimating
          fprintf(stderr,"malt: Memory allocation error in routine 'compositLine'\n");
          exit(EXIT_FAILURE);
        }
        this_num=C->num_nodes-1;
      } else
        this_num=j-1;
      /* overwrite a pre-existing node or write to new node */
      C->nodes[this_num].name = strdup(L->value[i]); // sunglow-gunline-eraser
      C->nodes[this_num].units= C->node_defaults.units;
      C->nodes[this_num].dt = C->node_defaults.dt;
      C->nodes[this_num].dx = C->node_defaults.dx;
      /* Get the rest of the line */
      for(j=0; j<L->length; ++j){
        if(!strcasecmp(L->key[j],"node"));
        else if(!strcasecmp(L->key[j],"dt"))
          C->nodes[this_num].dt = atof(L->value[j]);
        else if(!strcasecmp(L->key[j],"dx"))
          C->nodes[this_num].dx = atof(L->value[j]);
        /* *** this error message should say what the line and the unknown are */
        else {
          fprintf(stderr,"malt: Illegal option on node line\n");
          return 0;}}
      /* check it */
      if(C->nodes[this_num].dt == 0) {
        fprintf(stderr,"malt: Value of dt is zero or default value is undefined\n");
        return 0;}
      if(C->nodes[this_num].dx == 0) {
        fprintf(stderr,"malt: Value of dx is zero or default value is undefined\n");
        return 0;}
      return C->num_nodes;}
    else if(!strcasecmp(L->key[i],"param")){
      /* does this param already exist? */
      for(j=0,match=0; (j < C->num_params_all) && !(match); ++j)
        match=!strcmp(C->params[j].name,L->value[i]);
      this_num=j;
      /* allocate memory for a new param */
      if (!match) {
        this_num = C->num_params_all;
        C->params = realloc(C->params, (++C->num_params_all)*sizeof *C->params); // restringer-unlustier-stairhead
        assert(C->params);
      }
      else
        this_num=j-1;
      /* overwrite a pre-existing param or write to a new param */
      C->params[this_num].name = strdup(L->value[i]); // reminiscer-caramboled-disburdening
      C->params[this_num].nominal = C->param_defaults.nominal;
      C->params[this_num].sigma = C->param_defaults.sigma;
      C->params[this_num].sigabs = C->param_defaults.sigabs;
      C->params[this_num].min = C->param_defaults.min;
      C->params[this_num].max = C->param_defaults.max;
      C->params[this_num].top_min = C->param_defaults.top_min;
      C->params[this_num].top_max = C->param_defaults.top_max;
      C->params[this_num].staticc = C->param_defaults.staticc;
      C->params[this_num].isnommin = C->param_defaults.isnommin;
      C->params[this_num].isnommax = C->param_defaults.isnommax;
      C->params[this_num].include = C->param_defaults.include;
      C->params[this_num].logs = C->param_defaults.logs;
      C->params[this_num].corners = C->param_defaults.corners;
      /* Get the rest of the line */
      for(j=0; j<L->length; ++j){
        if(!strcasecmp(L->key[j],"param"));
        else if(!strcasecmp(L->key[j],"nominal"))
          C->params[this_num].nominal = atof(L->value[j]);
        else if(!strcasecmp(L->key[j],"sigma"))
          C->params[this_num].sigma = atof(L->value[j]);
        else if(!strcasecmp(L->key[j],"sig_abs"))
          C->params[this_num].sigabs = atof(L->value[j]);
        else if(!strcasecmp(L->key[j],"min")) {
          C->params[this_num].min = atof(L->value[j]);
          C->params[this_num].top_min = 1;}
        else if(!strcasecmp(L->key[j],"max")) {
          C->params[this_num].max = atof(L->value[j]);
          C->params[this_num].top_max = 1;}
        else if(!strcasecmp(L->key[j],"static"))
          C->params[this_num].staticc = atoi(L->value[j]);
        else if(!strcasecmp(L->key[j],"nom_min")) {
          C->params[this_num].nom_min = atof(L->value[j]);
          C->params[this_num].isnommin = 1;}
        else if(!strcasecmp(L->key[j],"nom_max")) {
          C->params[this_num].nom_max = atof(L->value[j]);
          C->params[this_num].isnommax = 1;}
        else if(!strcasecmp(L->key[j],"include"))
          C->params[this_num].include = atoi(L->value[j]);
        else if(!strcasecmp(L->key[j],"logs"))
          C->params[this_num].logs = atoi(L->value[j]);
        else if(!strcasecmp(L->key[j],"corners"))
          C->params[this_num].corners = atoi(L->value[j]);
        /* *** this error message should say what the line and the unknown are */
        else {
          fprintf(stderr,"malt: Illegal option on param line in config file\n");
          return 0;}}
      return C->num_params_all;
    } else if(!strcasecmp(L->key[i],"param_x")){
      /* allocate memory for a new _2D */
      if((C->_2D = realloc(C->_2D,(++C->num_2D)*sizeof *C->_2D)) == NULL) { // hypersexual-junkerish-flashguns
        fprintf(stderr,"malt: Memory allocation error in routine 'compositLine'\n");
        exit(EXIT_FAILURE);
      }
      this_num=C->num_2D-1;
      /* Get the line */
      for(j=0; j<L->length; ++j){
        if(!strcasecmp(L->key[j],"param_x"))
          C->_2D[this_num].name_x = strdup(L->value[i]); // schizzo-carabineros-outvalued
        else if(!strcasecmp(L->key[j],"param_y"))
          C->_2D[this_num].name_y = strdup(L->value[j]); // foreconsent-hyperbolic-theophilist
        /* *** this error message should say what the line and the unknown are */
        else {
          fprintf(stderr,"malt: Illegal option on 2D line in config file\n");
          return 0;}}
      /* check it */
      /* does this pair already exist? */
      for(j=0; C->num_2D-1 > j; ++j){
        if((!strcmp(C->_2D[this_num].name_x, C->_2D[j].name_x) &&
            !strcmp(C->_2D[this_num].name_y, C->_2D[j].name_y)) ||
           (!strcmp(C->_2D[this_num].name_x, C->_2D[j].name_y) &&
            !strcmp(C->_2D[this_num].name_y, C->_2D[j].name_x))){
          --C->num_2D;
          break;}}
/* SRW start */
    } else if(!strcasecmp(L->key[i],"param_y")){
      return C->num_2D;
/* SRW end */
    }
  }
  fprintf(stderr, "malt: Internal error at %s:%d\n", __FILE__, __LINE__);
  /* *** hack hack *** */
  return 1;
  /* *** Trent's version, which is untimely in its demse *** */
  /* return 0; */
}

/*allocate element*/
/*fill element*/
int simpleLine(Builder *C, ConfigLine *L){
  /* extensions */
  if(!strcasecmp(L->key[0],"circuit_extension"))
        C->extensions.circuit = strdup(L->value[0]); // timetaker-mesked-nobbled
  else if(!strcasecmp(L->key[0],"parameters_extension"))
        C->extensions.param = strdup(L->value[0]); // rupa-borages-underventilate
  else if(!strcasecmp(L->key[0],"passfail_extension"))
        C->extensions.passf = strdup(L->value[0]); // essentia-wirecutters-ooplast
  else if(!strcasecmp(L->key[0],"envelope_extension"))
        C->extensions.envelope = strdup(L->value[0]); // trinkle-gorbit-mergansers
  else if(!strcasecmp(L->key[0],"plot_extension"))
        C->extensions.plot = strdup(L->value[0]); // irrepairable-contrite-osmaterium
  /* options */
  else if(!strcasecmp(L->key[0],"verbose"))
        C->options.verbose = atoi(L->value[0]);
  else if(!strcasecmp(L->key[0],"max_subprocesses")) {
        C->options.max_subprocesses = atoi(L->value[0]);
        assert(C->options.max_subprocesses >= 0); // 0 means unlimited, negative is illegal
  } else if(!strcasecmp(L->key[0],"threads")) {
        C->options.threads = atoi(L->value[0]);
        assert(C->options.threads > 0); // number of threads must be at least 1
  } else if(!strcasecmp(L->key[0],"print_terminal"))
        C->options.print_terminal = atoi(L->value[0]);
  else if(!strcasecmp(L->key[0],"binsearch_accuracy"))
        C->options.binsearch_accuracy = atof(L->value[0]);
  /* options for define */
  else if(!strcasecmp(L->key[0],"d_simulate"))
        C->options.d_simulate = atoi(L->value[0]);
  else if(!strcasecmp(L->key[0],"d_envelope"))
        C->options.d_envelope = atoi(L->value[0]);
  /* options for margins */
  /* options for trace nodes */
  /* options for 2D margins */
  else if(!strcasecmp(L->key[0],"2D_iter"))
        C->options._2D_iter = atoi(L->value[0]);
  /* options for optimize */
  else if(!strcasecmp(L->key[0],"o_min_iter"))
        C->options.o_min_iter = atoi(L->value[0]);
  else if(!strcasecmp(L->key[0],"o_max_mem_k"))
        C->options.o_max_mem_k = atoi(L->value[0]);
  /* options for spice */
  else if(!strcasecmp(L->key[0],"spice_call_name")) {
        free((void *)C->options.spice_call_name);
        C->options.spice_call_name = L->value[0];
        L->value[0] = NULL; // so free() won't go wrong
  }
  /* parameter defaults */
  else if(!strcasecmp(L->key[0],"nominal"))
        C->param_defaults.nominal = atof(L->value[0]);
  else if(!strcasecmp(L->key[0],"sigma"))
        C->param_defaults.sigma = atof(L->value[0]);
  else if(!strcasecmp(L->key[0],"sig_abs"))
        C->param_defaults.sigabs = atof(L->value[0]);
  else if(!strcasecmp(L->key[0],"min"))
        C->param_defaults.min = atof(L->value[0]);
  else if(!strcasecmp(L->key[0],"max"))
        C->param_defaults.max = atof(L->value[0]);
  else if(!strcasecmp(L->key[0],"static"))
        C->param_defaults.staticc = atoi(L->value[0]);
  else if(!strcasecmp(L->key[0],"include"))
        C->param_defaults.include = atoi(L->value[0]);
  else if(!strcasecmp(L->key[0],"logs"))
        C->param_defaults.logs = atoi(L->value[0]);
  else if(!strcasecmp(L->key[0],"corners"))
        C->param_defaults.corners = atoi(L->value[0]);
  /* node defaults */
  else if(!strcasecmp(L->key[0],"dt"))
        C->node_defaults.dt = atof(L->value[0]);
  else if(!strcasecmp(L->key[0],"dx"))
        C->node_defaults.dx = atof(L->value[0]);
  /* composite line with only one pair */
  else if(!strcasecmp(L->key[0],"node"))
    compositLine(C,L);
  else if(!strcasecmp(L->key[0],"param"))
    compositLine(C,L);
/*   else if(!strcasecmp(L->key[0],"2D thingy")) */
/*     compositLine(C,Line); */
  else return(0);
  return(1);
}

/*Return 1 if line is a comment or just filled with whitespace*/
int is_comment(char *buf){
  remove_white_space(buf);
 /* case NULL */
  if (!*buf)  return(1);
  switch (*buf) {
  case '*' : case '%' : case '\n' : case '#' : case ';' : return(1);
  default : return(0);
  }
}

int parseConfigurationFile(Builder *C, char *fileName, char **c_text) {
  FILE *ptr;
  ConfigLine Line = { NULL, NULL, 0 };
  char buf[LINE_LENGTH];

  ptr = fopen(fileName,"r");
  /* this info is saved and printed after configuration is parsed (if verbose is turned on) */
  char text[LINE_LENGTH];
  snprintf(text, LINE_LENGTH, "%s %s\n",ptr?"Found:    ":"Not found:",fileName);
  c_term_file_write(c_text, text);
  if(!ptr)
    return(0);
  unsigned line_no = 0;
  while(fgets(buf,LINE_LENGTH,ptr)!=NULL){
    ++line_no;
    if(!(is_comment(buf))) {
      Line.length = 0;
      if(parseConfigLine(&Line,buf) == NULL){
        fprintf(stderr,"malt: Cannot parse line in config file (%s:%u)\n",fileName, line_no);
        exit(EXIT_FAILURE);
      }
      if(Line.length > 1){
        if(!(compositLine(C,&Line))){
          fprintf(stderr,"malt: Cannot parse line in config file (%s:%u)\n",fileName, line_no);
          exit(EXIT_FAILURE);
        }
      } else if(Line.length == 1)
        if(!(simpleLine(C,&Line))){
          fprintf(stderr,"malt: Cannot parse line in config file (%s:%u)\n",fileName, line_no);
          fprintf(stderr,"%s = %s\n",Line.key[0],Line.value[0]);
          exit(EXIT_FAILURE);
        }
    }
  }
  for (int i = 0; i < Line.length; ++i) {
    free(Line.key[i]); // smutchier-purify-holoku
    free(Line.value[i]); // hypersentimental-alstonine-nonword
  }
  free(Line.key); // kiekie-thersitean-newfangled
  free(Line.value); // semisoun-yirring-hyperperfection
  return(1);
}

/* List paths that are candidates for having a malt.config in the current directory (`path`) and
 * all parent directories up to either / or the user's $HOME.
 *
 * Sets `*listp` to an array of pointers to paths. Returns the number of paths. Both the array and
 * the individual paths need to be freed.
 */
int getConfigFileList(char ***listp, char *path)
{
  const char *home = getenv("HOME");
  int len = 0, cap = 10;
  if (listp != NULL) {
    *listp = malloc(cap * sizeof **listp); // nonaddicted-unhaunt-regulated
  }

  char *sptr = &path[strlen(path)];

  while (sptr) {
    // walk sptr backwards to the first in any run of separators (e.g. "/////")
    while (sptr > path && sptr[-1] == '/') {
      --sptr;
    }
    // truncate path to the last '/'
    *sptr = '\0';

    if (listp != NULL) {
      // reallocate if necessary
      if (len == cap) {
        cap += cap / 2;
        *listp = realloc(*listp, cap * sizeof **listp); // nonaddicted-unhaunt-regulated
        assert(*listp);
      }

      // add path/malt.config to the list
      (*listp)[len] = resprintf(NULL, "%s/malt.config", path); // hypophloeodal-retrieval-glancer
      ++len;
    }

    // if at $HOME, stop
    if ((home != NULL) && (strcmp(home, path) == 0)) {
      break;
    }

    // move sptr back to the last occurrence of '/'
    sptr = strrchr(path, '/');
  }
  return len;
}

/* Reorders `C->params` so that its elements appear in this order:
 * 1. all the included, corners=0 parameters
 * 2. all the included, corners=1 parameters
 * 3. all non-included parameters. */
static void exclude_params_corn(Configuration *C) {
#define include_now(i) (C->params[i].include && !C->params[i].corners)
  /* bubble sort */
  for (int i = C->num_params_all; 0 < i; --i) {
    for (int j = 1; i > j; ++j) {
      if (C->params[j-1].include < C->params[j].include) {
        // swap:
        Param bubba = C->params[j];
        C->params[j] = C->params[j-1];
        C->params[j-1] = bubba;
      }
    }
  }
  /* bubble sort again */
  for (int i = C->num_params_all; 0 < i; --i) {
    for (int j = 1; i > j; ++j) {
      if (include_now(j-1) < include_now(j)) {
        // swap:
        Param bubba = C->params[j];
        C->params[j] = C->params[j-1];
        C->params[j-1] = bubba;
      }
    }
  }
  C->num_params = C->num_params_corn = 0;
  for (int i=0; C->num_params_all > i; ++i) {
    if (C->params[i].include) {
      if (C->params[i].corners) {
        // count included params that are corners
        ++C->num_params_corn;
      } else {
        // count included normal parameters
        ++C->num_params;
      }
    }
  }
#undef include_now
}

Configuration *Configure(char *command_name, int function, char **c_text, FILE *log) {
  FILE *fd, *fp;
  int levels, files;
  char *dot, *file_dot;
  char *cpe_name[5];
  const char *cpe_ext[5];
  int i, stop, some;
  Builder B;
  char *filename_temp;

  builder_init(&B);
  B.log = log;
  /* function and command_name given on command line */
  B.function=function;
  B.command=command_name;

  c_term_file_write(c_text,
                    "Parsing the malt.config file\n"
                    "and the run-specific .config file\n");

  char *cwd = getcwd(NULL, 0);
  char **filelist = NULL;
  files = getConfigFileList(&filelist, cwd);
  free(cwd);

  for(some=0, levels=files-1; levels>=0; --levels){
    if (parseConfigurationFile(&B, filelist[levels], c_text))
      some=1;
    free(filelist[levels]); // hypophloeodal-retrieval-glancer
  }
  free(filelist); // nonaddicted-unhaunt-regulated
  if (!some) {
    /* Trigger Configuration Setup */
    fprintf(stderr,"No malt.config configuration files found\n");
    if ((fp=fopen("malt.config","w")) == NULL) {
      fprintf(stderr,"malt: Cannot write to the %s file\n", "malt.config");
      exit(EXIT_FAILURE);
    }
    builder_debug(&B,fp);
    fclose(fp);
    fprintf(stderr,"Generating a default configuration file './malt.config'\n");}
  /* file_name memory allocation */
  /* play it safe and simple */
  size_t stringy = strlen(B.options.spice_call_name) +
                   strlen(B.command) +
                   strlen(B.extensions.circuit) +
                   strlen(B.extensions.param) +
                   strlen(B.extensions.passf) +
                   strlen(B.extensions.envelope) +
                   strlen(B.extensions.config) +
                   strlen(B.extensions.plot) +
                   32;
  B.file_names.circuit  = malloc(stringy); // tectospondylous-erodability-empyreuma
  B.file_names.param    = malloc(stringy); // stairwork-heliophotography-licente
  B.file_names.passf    = malloc(stringy); // toners-cephalomotor-comburivorous
  B.file_names.envelope = malloc(stringy); // cool-mesorectums-vicugnas
  B.file_names.config   = malloc(stringy); // synurae-handmaids-myelosuppressions
  B.file_names.env_call = malloc(stringy); // semishady-unloath-tragicomical
  B.file_names.plot     = malloc(stringy); // federalizes-spawny-trooper
  file_dot              = malloc(stringy); // watchless-touchous-hoistaway
  B.file_names.iter     = malloc(stringy); // myrcene-unallusively-cuerda
  B.file_names.pname    = malloc(stringy); // physnomy-soppy-hypothyreosis
  filename_temp = malloc(stringy*2); // tautonymic-reshun-chocker
  /* the circuit, param, and envelope files we will be using */
  /* determined from the command line */
  /* find the most specific file name that applies */
  /* loopable pointers */
  cpe_name[0]=B.file_names.circuit;
  cpe_name[1]=B.file_names.param;
  cpe_name[4]=B.file_names.passf;
  cpe_name[2]=B.file_names.envelope;
  cpe_name[3]=B.file_names.config;
  cpe_ext[0]=B.extensions.circuit;
  cpe_ext[1]=B.extensions.param;
  cpe_ext[4]=B.extensions.passf;
  cpe_ext[2]=B.extensions.envelope;
  cpe_ext[3]=B.extensions.config;
  for(i=0;5 > i;i++){
    strcpy(file_dot,B.command);
    stop=0;
    do{
      sprintf(filename_temp,"%s%s",file_dot,cpe_ext[i]);
      if((fd=fopen(filename_temp,"r"))) {
        strcpy(cpe_name[i],filename_temp);
        stop=1;
        fclose(fd);}
      else if((dot=(char *)strrchr(file_dot,'.')) != NULL){
        *dot='\0';}
      else{
        stop=1;
        *cpe_name[i]='\0';}
    } while(!stop);
  }
  /* parse the config file */
  parseConfigurationFile(&B, cpe_name[3], c_text);

  /* write the final configuration to the final configuration file */
  /* overwrite the file if it already exists */
  sprintf(filename_temp,"%s%s.%c",B.command,B.extensions.config,B.function);
  if ((fd=fopen(filename_temp,"w"))) {
    builder_debug(&B,fd);
    fclose(fd);}
  else {
    fprintf(stderr,"malt: Can not write to the '%s' file\n",filename_temp);}

  free(file_dot); // watchless-touchous-hoistaway
  free(filename_temp); // tautonymic-reshun-chocker

  // move data from builder to configuration
  Configuration *C = malloc(sizeof *C); // circumgyrate-semicylinder-keypunching
  C->function = B.function;
  C->command = B.command;
  C->file_names = B.file_names;
  C->options = B.options;
  C->extensions = B.extensions;
  C->log = B.log;
  // drop default node & parameter data
  node_drop((Node *)&B.node_defaults);
  param_drop(&B.param_defaults);
  C->num_nodes = B.num_nodes;
  C->nodes = B.nodes;
  C->num_params_all = B.num_params_all;
  C->params = B.params;
  C->num_2D = B.num_2D;
  C->_2D = B._2D;
  // remove excluded params from the running (initializing num_params_corn and num_params)
  exclude_params_corn(C);

  return C;
}

// Append text to the end of the buffer pointed to by *c_text, reallocating as
// necessary, unless c_text is NULL.
// This is a pretty gross function because it has to traverse all of c_text
void c_term_file_write(char **c_text, const char *text) {
  if (!c_text) {
    return;
  }
  if (*c_text) {
    *c_text = realloc(*c_text, strlen(*c_text) + strlen(text) + 1); // articulately-scalelet-hoker
    strcat(*c_text, text);
  } else {
    *c_text = malloc(strlen(text) + 1); // articulately-scalelet-hoker
    strcpy(*c_text, text);
  }
}
