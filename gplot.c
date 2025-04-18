// vi: ts=2 sts=2 sw=2 et tw=100
/* gnuplot functions */
#include "gplot.h"
#include "config.h"
#include "malt.h"
#include "space.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* void viewSim(Configuration *C){ */
/*   int i; */
/*   FILE *strm; */

/*   if(!(strm = popen("gnuplot","w"))){ */
/*     fprintf(stderr,"Could not run the 'gnuplot' executable\n");} */
/*   else { */
/*     fprintf(strm,"unset key\n"); */
/*     fprintf(strm,"set xlabel 'time (s)'\n"); */
/*     for(i=0;i<C->num_nodes;i++){ */
/*       fprintf(strm,"set label 1 '%s' at graph 0.05, graph 0.95\n",C->nodes[i]->name); */
/*       fprintf(strm,"set label 2 'Level tolerance is %g' at graph 0.05, graph
 * 0.90\n",C->nodes[i]->dx); */
/*       fprintf(strm,"set label 3 'Jitter tolerance is %g' at graph 0.05, graph
 * 0.85\n",C->nodes[i]->dt); */
/*       /\* gnuplot is plotting data from the file 'gnudata' *\/ */
/*       fprintf(strm,"plot 'gnudata' index 0 using 1:%d with linespoints 2,'gnudata' index 0 using
 * 1:%d with linespoints 1,'gnudata' index 0 using 1:%d with linespoints 2,'gnudata' index 1 using
 * 1:%d with linespoints 1\n",3*i+2,3*i+3,3*i+4, i+2); */
/*       fflush(strm); */
/*       getchar(); */
/*     } */
/*     pclose(strm); */
/*   } */
/* } */

void plot2(Configuration *C, Space *S)
{
  int i;
  double a, b, c;
  FILE *fp;
  FILE *strm;

  /* write the gnuplot script to a file so it can be user-modified later */
  fp = new_file_by_type(C, Ft_Gnuplot);

  const char *scriptname = malt_filename(C, Ft_Gnuplot);
  info("Writing Gnuplot script to '%s'\n", scriptname);
  fprintf(fp, "# Gnuplot script\n# User-editable\n# Run 'gnuplot %s'\n\n", scriptname);
  fprintf(fp, "set size 0.75, 1.00\nunset key\n\n");
  /* for each plot... */
  for (i = 0; C->num_2D > i; ++i) {
    fprintf(fp, "%sset logscale x\n%sset logscale y\n",
            C->params[C->_2D[i].param_x].logs ? "" : "un",
            C->params[C->_2D[i].param_y].logs ? "" : "un");
    fprintf(fp, "set xlabel '%s'\nset ylabel '%s'\n", C->params[C->_2D[i].param_x].name,
            C->params[C->_2D[i].param_y].name);
    fprintf(fp, "set xrange [ %g: %g ]\nset yrange [ %g: %g ]\n",
            physspace(C->params[C->_2D[i].param_x].min, C, C->_2D[i].param_x),
            physspace(C->params[C->_2D[i].param_x].max, C, C->_2D[i].param_x),
            physspace(C->params[C->_2D[i].param_y].min, C, C->_2D[i].param_y),
            physspace(C->params[C->_2D[i].param_y].max, C, C->_2D[i].param_y));
    a = physspace(C->params[C->_2D[i].param_x].min, C, C->_2D[i].param_x);
    b = physspace(S[C->_2D[i].param_x].centerpnt, C, C->_2D[i].param_x);
    c = physspace(C->params[C->_2D[i].param_x].max, C, C->_2D[i].param_x);
    fprintf(fp, "%sset xtics ('%.3f' %.3f, '%.3f' %.3f, '%.3f' %.3f)\n",
            /* FIXME: The following doesn't make sense -- params[*] was a pointer but should never
               be null: C->params[C->_2D[i].param_x]?"":"# ",
             */
            "", a, a, b, b, c, c);
    a = physspace(C->params[C->_2D[i].param_y].min, C, C->_2D[i].param_y);
    b = physspace(S[C->_2D[i].param_y].centerpnt, C, C->_2D[i].param_y);
    c = physspace(C->params[C->_2D[i].param_y].max, C, C->_2D[i].param_y);
    fprintf(fp, "%sset ytics ('%.3f' %.3f, '%.3f' %.3f, '%.3f' %.3f)\n",
            /* FIXME: The following doesn't make sense -- params[*] was a pointer but should never
               be null: C->params[C->_2D[i].param_y]?"":"# ",
             */
            "", a, a, b, b, c, c);
    fprintf(fp, "plot '%s' using %d:%d with linespoints lc %d pt %d\n", malt_filename(C, Ft_Data),
            2 * i + 1, 2 * i + 2, i + 1, i + 1);
    /* ought to be able to say 'set term x11' instead of pause but our gnuplot is not allowing it
     * right now */
    fprintf(fp, "pause 3\n\n");
  }
  fclose(fp);
  /* execute the script */
  if (!(strm = popen("gnuplot", "w"))) {
    warn("Could not run the 'gnuplot' executable\n");
  } else {
    fprintf(strm, "load '%s'\n", scriptname);
    pclose(strm);
  }
}
