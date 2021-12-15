#include "call_spice.h"
#include "config.h"
#include "define.h"
#include "malt.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

int  define(Configuration*);
int  readData(Configuration*, Data*, int *scramble);
void spiceBounds(Configuration *, Data *, int *scramble);
void spicePlot(Configuration *);
void dumpData(Data*);
void bound(Data*, Configuration*);
void dumpBounds(Configuration*, Data*, char*);
void viewSim(Configuration*);

int call_def(Configuration *C) {
  int all_good;

  /* create pname file, which needs to know currently excluded parameters */
  /* not applicable, include them all */
  pname(C);
  /* do it */
  all_good=define(C);
  /* clean up temporary files */
  unlink(C->file_names.pname);
  return all_good;
}

int define(Configuration *C)
{
  Data *theData;
  int *scramble = malloc(C->num_nodes * sizeof *scramble);
  int ret = 0;

  /* run spice which writes vectors to a file */
  if(C->options.d_simulate) {
    if(C->num_nodes == 0){
      fprintf(stderr,"No circuit nodes have been defined\n");
      exit(EXIT_FAILURE);
    }
    /* *** need to check for spice executable, etc */
    call_spice(C, 0, NULL, NULL, "define.call", "define.return"); // TODO: verify these NULLs are ok ~ntj
  }
  /* read the spice simulation vectors file */
  if(C->options.d_envelope){
    /* memory allocation */
    theData = malloc(sizeof *theData); // mem:beslur
    /* do not change node definitions between simulate and envelope, or you are hosed */
    /* load vectors from spice print file */
    if (!readData(C, theData, scramble))
      goto fail;
    /* make the envelope */
    bound(theData,C);
    /* write envelope to a spice-loadable file */
    sprintf(C->file_names.envelope,"%s%s",C->command,C->extensions.envelope);
    sprintf(C->file_names.env_call,"%s%s.",C->command,C->extensions.envelope);
    spiceBounds(C, theData, scramble);
  }
  /* now we want to call a spice script to view envelope */
  spicePlot(C);
  ret = 1;
fail:
  /* clean up generic temporary files */
  free(scramble);
  return ret;
}

void spiceBounds(Configuration *C, Data *D, int *scramble)
{
  int i,j,k=0,m;
  FILE *fp, *fp2;
  int step_div, step_mod, step;

  if (!(fp = fopen(C->file_names.envelope,"w"))) {
    fprintf(stderr,"malt: Can not open the file '%s'\n", C->file_names.envelope);
    exit(EXIT_FAILURE);
  }
  /* print the header material */
  fprintf(fp,"Title: Envelopes\n");
  fprintf(fp,"Plotname: Transient\n");
  fprintf(fp,"No. Variables: %d\n",1+2*C->num_nodes);
  fprintf(fp,"No. Points: %d\n",D->length);
  fprintf(fp,"Variables:\n 0 time S\n");
  for (i=0; C->num_nodes > i; i++){
    m = scramble[i];
    fprintf(fp," %d hi%d %c\n",++k,m,C->nodes[m].units);
    fprintf(fp," %d lo%d %c\n",++k,m,C->nodes[m].units);
  }
  /* print data and envelope */
  fprintf(fp,"Values:\n");
  for (i=0; D->length-D->dtlength > i; i++){
    fprintf(fp,"%d\t%e\n", i, D->t[i]);
    for (j=0; C->num_nodes > j; j++){
      m = scramble[j];
      fprintf(fp,"\t%e\n\t%e\n", D->upper[m][i], D->lower[m][i]);}
  }
  fclose(fp);
  /* * * * file that loads the above file */
  if (!(fp2 = fopen(C->file_names.env_call,"w"))) {
    fprintf(stderr,"malt: Can not open the file '%s'\n", C->file_names.env_call);
    exit(EXIT_FAILURE);
  }
  /* header stuff */
  fprintf(fp2,"* the envelope of acceptable values for each trace\n* called by binsearch generic\n.control\n");
  fprintf(fp2,"* node math is legal, i.e. v(1)-v(2)\n* ...so long as the nodes also appear separately\n");
  fprintf(fp2,"* lo and hi are the envelope definition\n\n");
  /* node_hi */
  fprintf(fp2,"set node_hi = (");
  for(i=0; i < C->num_nodes; ++i)
    fprintf(fp2," tran1.hi%d ",i);
  fprintf(fp2,")\n");
  /* node_lo */
  fprintf(fp2,"set node_lo = (");
  for(i=0; i < C->num_nodes; ++i)
    fprintf(fp2," tran1.lo%d ",i);
  fprintf(fp2,")\n");
  /* step values */
  fprintf(fp2,"compose step_value values ");
  /* step 10 times, equal number of time units per step */
  /* step D->length times if D->length < 10 */
  step_div=(D->length-D->dtlength)/10;
  step_mod=(D->length-D->dtlength)%10;
  for(i=0; i < 10; ++i) {
    step=(i < step_mod)?step_div+1:step_div;
    if(step)  fprintf(fp2," %d",step);}
  fprintf(fp2,"\nload %s\nsetplot constants\n",C->file_names.envelope);
  fprintf(fp2,".endc\n");
  fclose(fp2);
}

void spicePlot(Configuration *C){
  int i;
  FILE *fp;

  sprintf(C->file_names.plot,"%s%s.%c",C->command,C->extensions.plot,C->function);
  if (!(fp = fopen(C->file_names.plot,"w"))) {
    fprintf(stderr,"malt: Can not open the file '%s'\n",C->file_names.plot);
    exit(EXIT_FAILURE);
  }
  /* print the file */
  fprintf(fp,"\n.control\nload %s.nom\nload %s\nsetplot tran1\nset group\nplot",
          C->command,C->file_names.envelope);
  for (i=0; C->num_nodes > i; i++)
    fprintf(fp," %s",C->nodes[i].name);
  fprintf(fp,"\nunset group\nset single\n");
  fprintf(fp,"set color2 = \"black\"\nset color3 = \"blue\"\nset color4 = \"black\"\n");
  for (i=0; C->num_nodes > i; i++){
    fprintf(fp,"plot tran2.hi%d %s tran2.lo%d\n",i,C->nodes[i].name,i);}
  fprintf(fp,".endc\n");
  fclose(fp);
  /* run spice */
  char *sys = resprintf(NULL, "%s %s", C->options.spice_call_name, C->file_names.plot); // mem:absolutist
  system(sys);
  free(sys); // mem:absolutist
}

void bound(Data *D, Configuration *C){
  int i,j,k;
  double r;
  for(k=0; k<C->num_nodes; ++k){
    if(isfinite(r = floor((C->nodes[k].dt)/D->tstep))) {
      D->dtlength = (int)r;
    } else {
      fprintf(stderr,"malt: Value of tstep (%g) or dt (%g) is zero\n", D->tstep, C->nodes[k].dt);
      exit(EXIT_FAILURE);
    }
    for(i=0; i<D->length; ++i){
      if( (i+D->dtlength < D->length) ) {
        D->upper[k][i+D->dtlength] = D->x[k][i+D->dtlength] + C->nodes[k].dx;
        D->lower[k][i+D->dtlength] = D->x[k][i+D->dtlength] - C->nodes[k].dx;}
      for(j=-D->dtlength; j<D->dtlength; ++j)
        if( ((i+j)>=0) && ((i+j)<D->length) ){
          r = (C->nodes[k].dx)*sqrt(1.-(j*D->tstep)*(j*D->tstep)/((C->nodes[k].dt)*(C->nodes[k].dt)));
          D->upper[k][i+j] = (D->upper[k][i+j] < D->x[k][i] + r ) ? D->x[k][i] + r : D->upper[k][i+j];
          D->lower[k][i+j] = (D->lower[k][i+j] > D->x[k][i] - r ) ? D->x[k][i] - r : D->lower[k][i+j];
        }
    }
  }
}

int readData(Configuration *C, Data *D, int *scramble)
{
  /* *** fixed length string again */
  FILE *fp;
  char line[1024], name[1024], units;
  int index, i, j, node_num;
  float x;
  int ret = 1;

  char *returnn = resprintf(NULL, "%s.nom", C->command); // mem:quininic
  if (!(fp = fopen(returnn,"r"))) {
    fprintf(stderr, "malt: Can not open %s\n",returnn);
    ret = 0;
    goto cleanup;
  }
  /* size data */
  do {
    if (!fgets(line,1024,fp)) {
      ret = 0;
      goto cleanup;
    }
    /* read number of nodes */
    if (!strncmp("No. Variables:",line,14))
      sscanf(line,"No. Variables:%d",&node_num);
    /* read number of points */
    else if (!strncmp("No. Points:",line,11))
      sscanf(line,"No. Points:%d",&D->length);
  } while (strncmp("Variables:",line,10));
  /* make sure number of nodes agrees with expectations */
  if ((--node_num) != C->num_nodes) {
    fprintf(stderr, "malt: Number of nodes in %s file is not right.\nMaybe node math does not work\n",returnn);
    ret = 0;
    goto cleanup;
  }
  /* throw away time line */
  if (!fgets(line,1024,fp)) {
    ret = 0;
    goto cleanup;
  }
  /* pick off the name & units for each vector */
  /* name is needed to unscramble the order */
  for (i=0; C->num_nodes > i; i++) {
    if (!fgets(line,1024,fp)) {
      ret = 0;
      goto cleanup;
    }
    // a bit hacky way to match `v(phi.XICHECK0)` to `phi.XICHECK0`
    if (2 != sscanf(line, "%*d %s %c", name, &units)) {
      fprintf(stderr, "malt: Malformed return file %s\n"
                      "    %s\n"
                      "    ^ (on this line)\n", returnn, line);
      ret = 0;
      goto cleanup;
    }
    int match = 0;
    for (j=0; (C->num_nodes > j) && !match; j++) {
      if (0 == strcmp(C->nodes[j].name, name)) {
        scramble[i] = j;
        match = 1;
      }
    }
    if (!match) {
      fprintf(stderr, "malt: Unidentified node \"%s\" in return file %s.\n",
        name, returnn);
      fprintf(stderr, "malt: The following nodes are known:\n");
      for (j = 0; (C->num_nodes > j) && !match; j++) {
        fprintf(stderr, "    %2d: \"%s\"\n", j, C->nodes[j].name);
      }
      ret = 0;
      goto cleanup;
    }
    if (C->nodes[scramble[i]].units != units) {
      fprintf(stderr, "malt: Mismatched units of node %s: expected %c, found %c (in file %s).\n",
        name, C->nodes[scramble[i]].units, units, returnn);
      ret = 0;
      goto cleanup;
    }
  }
  /* skip down to vector values */
  do {
    if (!fgets(line,1024,fp)) {
      ret = 0;
      goto cleanup;
    }
  } while (strncmp("Values:",line,7));
  /* allocate D */
  D->t = calloc(D->length, sizeof *D->t); // mem:downhanging
  D->x = malloc(C->num_nodes * sizeof *D->x); // mem:overdistantly
  D->upper = malloc(C->num_nodes * sizeof *D->upper); // mem:vixenlike
  D->lower = malloc(C->num_nodes * sizeof *D->lower); // mem:carucate
  for(i=0; i < C->num_nodes; ++i) {
    D->x[i] = calloc(D->length, sizeof **D->x); // mem:pachyhaemia
    D->upper[i] = calloc(D->length, sizeof **D->upper); // mem:untrustness
    D->lower[i] = calloc(D->length, sizeof **D->lower); // mem:preabundantly
  }
  /* read data */
  for (j=0; D->length > j; j++) {
    /* index & time line */
    if (!fgets(line,1024,fp)) {
      ret = 0;
      goto cleanup;
    }
    sscanf(line,"%d %f",&index, &x);
    D->t[index] = (double)x;
    /* vector line */
    for (i=0; C->num_nodes > i; i++) {
      if (!fgets(line,1024,fp)) {
        ret = 0;
        goto cleanup;
      }
      D->x[scramble[i]][index] = atof(line);}
  }
  D->tstep = D->t[1] - D->t[0];
cleanup:
  free(returnn); // mem:quininic
  return ret;
}
