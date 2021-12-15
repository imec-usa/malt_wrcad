// vi: ts=2 sts=2 sw=2 et tw=100
/* optimization function */
#include "optimize.h"
#include "call_spice.h"
#include "config.h"
#include "malt.h"
#include "marg_opt_yield.h"
#include "space.h"
#include "stat_math.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#define N (C->num_params)

static int optimize(Configuration *C, Space *S, const double *prhi, const double *prlo)
{
  Plane **plane = NULL;
  double **margpnts = NULL;
  double *facecenter, *vect;
  int *tang;
  short *pntstack;
  int plncount = 0, plnmemory = 0, pntmemory = 0, pntcount = 0;
  int pinc, minc;
  int stop = 0, iterate = 0, max_iterate, max_plncount, pc_bytes;
  int tangent;
  double radius, dum, dum2;
  int flag, flag2 = 0;
  double radiushi, radiuslo;
  int i, j, k;
  int big, oldcount = 0, lowcount = 0;
  int match, ii;
  int all_good = 1;
  double *improve;
  double *pc = malloc((N + C->num_params_corn) * sizeof *pc);  // mem:unwomanish
  double *direction = malloc(N * sizeof *direction);           // mem:bergs
  double *pr = malloc(N * sizeof *pr);

  /* Exceptions */
  if (N == 1) {
    S[0].centerpnt = (prhi[0] + prlo[0]) / 2.0;
    if (S[0].centerpnt < C->params[0].nom_min)
      S[0].centerpnt = C->params[0].nom_min;
    if (S[0].centerpnt > C->params[0].nom_max)
      S[0].centerpnt = C->params[0].nom_max;
    lprintf(C, "center  %7.3f%s\n", physspace(S[0].centerpnt, C, 0),
            (C->params[0].nom_min < S[0].centerpnt && C->params[0].nom_max > S[0].centerpnt) ? " "
                                                                                             : "*");
    radiushi = prhi[0] - S[0].centerpnt;
    radiuslo = S[0].centerpnt - prlo[0];
    radius = (radiushi < radiuslo) ? radiushi : radiuslo;
    lprintf(C, "axes    %7.3f\n", radius);
    lprintf(C, "Yield ...%10.6e%%\n", gauss_integral(radius, N));
    lprintf(C, "1-Yield ...%10.6e%%\n", gauss_integral_c(radius, N));
    return 1;
  }
  lprintf(C,
          "\nOptimization in %d dimensions with %d corners"
          "\nPoints Simplexes M(Sigma)   Parameter_Values\n",
          N, 1 << C->num_params_corn);

  /* Exit criteria things */
  /* limit memory to a 32GiByte */
  if (C->options.y_max_mem_k > 33554432)
    C->options.y_max_mem_k = 33554432;
  /* translate max_mem into max_plncount */
  /* plncount memory budget */
  // Plane:         (short)flag=2, (short)points=2*N, double(b)=8, *a=8, double a[]=8*N
  // center matrix: (int)vector=8, (double)matrix=(N+3)*8
  pc_bytes = 2 + 8 + 8 + 8 + (2 + 8) * N + 8 * (N + 3);
  /* bytes per simplex: (double)mean=8, (double)powm=8, (short)t=N*2 */
  max_plncount = (int)(((long)C->options.y_max_mem_k * 1024) /
                       (long)pc_bytes);  // fits within a [signed] int, which is plenty
  max_iterate = 32767;  // because pntstack is of type [signed] short, which is plenty

  /* memory allocation */
  tang = malloc((N + 1) * sizeof *tang);                        // mem:quiets
  facecenter = malloc(N * sizeof *facecenter);                  // mem:intertrafficked
  vect = malloc(N * sizeof *vect);                              // mem:isostere
  pntstack = malloc(N * sizeof *pntstack);                      // mem:diaheliotropically
  improve = malloc((C->options.o_min_iter) * sizeof *improve);  // mem:knackwursts
  /* allocate first number of planes */
  pinc = (1 << N) + 1000;                                            // =2**N+1000
  if ((plane = plane_malloc(plane, &plnmemory, pinc, N)) == NULL) {  // mem:nucivorous
    perror("hurled whilest allocating planes\n");
    exit(EXIT_FAILURE);
  }
  minc = 2 * N + 100;
  margpnts = margpnts_malloc(margpnts, &pntmemory, minc, N);  // mem:crystallic
  /* initialize */

  /* store margins */
  pntcount = 2 * N;
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      margpnts[2 * i][j] = S[j].centerpnt;
      margpnts[2 * i + 1][j] = S[j].centerpnt;
    }
    margpnts[2 * i][i] = prlo[i];
    margpnts[2 * i + 1][i] = prhi[i];
  }
  /* pick initial combination of C->num_params*2 points and make the hull */
  intpickpnts(pntstack, C, S, 0, plane, margpnts, &plncount, pntcount);
  /*   lprintf(C, "simplex elements: %i total\n\n",plncount); */
  /* LOOP */
  do {
    ++iterate;
    if (!(tangent = center(C, S, plane, tang, plncount, &radius))) {
      fprintf(stderr, "malt: Numerical error in routine center");
      all_good = 0;
      goto cleanup;
    }
    if ((big = findface(N, facecenter, plane, tang, plncount, tangent)) == -1) {
      lprintf(C, "Simplex cannot inflate further: optimization complete\n");
      stop = 1;
    }
    if (!stop) {
      /* save search vector for quotation later */
      for (i = 0; N > i; ++i) {
        vect[i] = -plane[big]->a[i];
      }
      /* check that we are not against search boundary else flag plane */
      /* report the boundary that is too close */
      /* this is probably ineffective for high-dimensional space */
      for (i = 0, flag = -1; flag == -1 && N > i; ++i) {
        dum = facecenter[i];
        /* maybe this should get bigger for high-dimensional space */
        dum2 = C->options.binsearch_accuracy * S[i].centerpnt;
        if (dum > C->params[i].max - dum2)
          flag = i, flag2 = 1;
        if (dum < C->params[i].min + dum2)
          flag = i, flag2 = 0;
      }
      if (flag != -1) {
        plane[big]->flag = 1;
        /* i is indexed to 1 in the following statement */
        /* cause how it bust out the the for loop above */
        /* thusly we must subtract one from the index !!!!!!!!!!!!! */
        lprintf(C, "warning: Param '%s' limited by the binary search '%s'\n", C->params[i - 1].name,
                flag2 ? "max" : "min");
      } else {
        lowcount = oldcount = plncount;
        /* addpoint_corners needs to know the center (C->params[i].pc) */
        /*                  and the search direction (C->params[i].direction) */
        for (i = 0; N > i; ++i) {
          pc[i] = facecenter[i];
          direction[i] = plane[big]->a[i];
        }
        /* make sure there's room to store the new point */
        if (pntcount + 1 >= pntmemory)
          margpnts = margpnts_malloc(margpnts, &pntmemory, C->num_params * 2,
                                     C->num_params);  // mem:crystallic
        /* flag plane if not convex */
        if (addpoint_corners(C, S, NULL, margpnts[pntcount], pc, direction) == 0.0) {
          plane[big]->flag = 1;
          lprintf(C, "Convexity fail at point  ");
          for (i = 0; N > i; ++i) {
            lprintf(C, "%8.3f ", physspace(pc[i], C, i));
          }
          lprintf(C, "\n\n");
        } else {
          /* store the new point */
          ++pntcount;
          /* flag planes the new point is exterior to */
          for (j = 0; oldcount > j; ++j) {
            for (i = 0, dum = plane[j]->b; N > i; ++i) {
              dum += plane[j]->a[i] * margpnts[pntcount - 1][i];
            }
            if (dum < 0.0) {
              plane[j]->flag = 2;
              --lowcount;
            }
          }
          /* find new planes, using new point and old-plane-point combinations */
          for (k = 0; oldcount > k; ++k) {
            if (plane[k]->flag == 2) {
              for (i = 0; i < N; ++i) {
                for (j = 0; j < N; ++j) {
                  pntstack[j] = plane[k]->points[j];
                }
                pntstack[i] = pntcount - 1;
                /* allocate more planes in memory */
                if (plncount == plnmemory) {
                  pinc = plncount / 3 + 1;  // 30% more everytime
                  if ((plane = plane_malloc(plane, &plnmemory, pinc, N)) ==
                      NULL) {  // mem:nucivorous
                    perror("No more memory\n");
                    exit(EXIT_FAILURE);
                  }
                }
                makeaplane(pntstack, C, S, plane, margpnts, &plncount, pntcount);
              }
            }
          }
          /* flag duplicate planes (good e.g. if the operating region has a planar surface) */
          /* pick a plane */
          for (k = oldcount; plncount - 1 > k; ++k) {
            /* pick another plane */
            for (j = k + 1; plncount > j; ++j) {
              /* pick a point of the first plane */
              for (i = 0, match = 1; N > i && match; ++i) {
                /* does the point exist on the other plane? */
                for (ii = 0, match = 0; N > ii && !match; ++ii) {
                  match = (plane[k]->points[i] == plane[j]->points[ii]);
                }
              }
              if (match == 1) {
                plane[k]->flag = 2;
                /* printf("ODD: Duplicate planes were created!!\n"); */
              }
            }
          }
          /* delete old planes */
          /* now new planes may be flagged too */
          for (j = 0; plncount > j;) {
            if (plane[j]->flag == 2) {
              --plncount;
              plane[j]->b = plane[plncount]->b;
              plane[j]->flag = plane[plncount]->flag;
              for (i = 0; N > i; ++i) {
                plane[j]->a[i] = plane[plncount]->a[i];
                plane[j]->points[i] = plane[plncount]->points[i];
              }
            } else
              ++j;
          }
        }
      }
    }
    /* check if loop should terminate */
    /* maximum memory */
    if (!stop && plncount > max_plncount) {
      stop = 1;
      lprintf(C,
              "\nMemory usage has reached the limit of o_max_mem_k = %d (kB)\n\nOptimization "
              "Interrupted\n\n",
              C->options.y_max_mem_k);
    }
    /* maximum iterations */
    if (!stop && iterate == max_iterate) {
      stop = 1;
      lprintf(C,
              "\nIterations has reached the internal limit of %d\n\nOptimization Interrupted\n\n",
              max_iterate);
    }
    /* ITERATE file gone */
    if (!stop && !checkiter(C)) {
      stop = 1;
    }
    /* radius increasing slowly */
    if (!stop && (iterate > C->options.o_min_iter)) {
      if ((radius - improve[iterate % C->options.o_min_iter]) < C->options.binsearch_accuracy) {
        stop = 1;
        lprintf(C,
                "Margin has increased less than binsearch_accuracy=%.3f in the last o_min_iter=%d "
                "iterations:\nOptimization interrupted\n",
                C->options.binsearch_accuracy, C->options.o_min_iter);
      }
    }
    improve[iterate % C->options.o_min_iter] = radius;
    /* print the iteration */
    /* points, simplexes, M(sigma) */
    lprintf(C, "%4i %10i %7.2f  ", pntcount, plncount, radius);
    for (i = 0; N > i; ++i) {
      lprintf(C, "%8.3f%s", physspace(S[i].centerpnt, C, i),
              (C->params[i].nom_min < S[i].centerpnt && C->params[i].nom_max > S[i].centerpnt)
                  ? " "
                  : "*");
    }
    lprintf(C, "\n");

  } while (!stop);
  /* end loop */
  /* inscribe one more time */
  if (!(tangent = center(C, S, plane, tang, plncount, &radius))) {
    fprintf(stderr, "malt: Numerical error in routine center\n");
    all_good = 0;
    goto cleanup;
  }
  /* calculate criticalness */
  for (i = 0; N > i; ++i) {
    vect[i] = 0.0;
    for (j = 0; tangent > j; ++j) {
      vect[i] += (plane[tang[j]]->a[i]) * (plane[tang[j]]->a[i]);
    }
    vect[i] /= tangent;
  }

  /* some diagnostics */
  lprintf(C, "\nMemory Used (kB): %ld\n", ((long)plncount * (long)pc_bytes) / 1024 + 1);  // go long
  /* convexity check along critical vectors */
  lprintf(C, "\nConvexity: The distance ratio (simplex/ellipsoid) should be unity or greater"
             "\n                 Ratio:   Parameter_Values\n");
  for (j = 0; tangent > j; ++j) {
    for (i = 0; N > i; ++i) {
      /* find the closest boundary and calculate the search points */
      pc[i] = S[i].centerpnt;
      direction[i] = plane[tang[j]]->a[i];
    }
    /* find boundary */
    if (addpoint_corners(C, S, NULL, pr, pc, direction) == 0.0) {
      nrerror("Circuit failed for optimized parameter values\n");
    } else {
      for (i = 0, dum = 0.0; N > i; ++i) {
        dum += pow(pr[i] - S[i].centerpnt, 2);
      }
      dum = sqrt(dum);
      lprintf(C, "%22.3f:", (dum / radius));
      for (i = 0; N > i; ++i) {
        lprintf(C, "%8.3f ", physspace(pr[i], C, i));
      }
    }
    lprintf(C, "\n");
  }
  /* memory unallocation */
cleanup:
  free(tang);                          // mem:quiets
  free(facecenter);                    // mem:intertrafficked
  free(vect);                          // mem:isostere
  free(pntstack);                      // mem:diaheliotropically
  free(improve);                       // mem:knackwursts
  plane_free(plane, plnmemory);        // mem:nucivorous
  margpnts_free(margpnts, pntmemory);  // mem:crystallic
  free(pc);                            // mem:unwomanish
  free(direction);                     // mem:bergs
  free(pr);
  return all_good;
}

int call_opt(Configuration *C)
{
  int all_good = 1;

  /* initialize */
  Space *S = malloc(C->num_params_all * sizeof *S);
  if (!initspace(C, S))
    return 0;
  /* create the iterate file */
  makeiter(C, 'o');
  pname(C);
  double *prhi = malloc(C->num_params * sizeof *prhi);
  double *prlo = malloc(C->num_params * sizeof *prlo);
  /* Exceptions */
  if (N > 10) {
    fprintf(stderr, "%d non-corner parameters is greater than the 10 allowed.\n", N);
    all_good = 0;
    goto cleanup;
  }
  /* margins */
  if (!margins(C, S, prhi, prlo)) {
    all_good = 0;
    goto cleanup;
  }
  /* optimize */
  if (!optimize(C, S, prhi, prlo)) {
    all_good = 0;
    goto cleanup;
  }
  /* margins again at the end */
  if (!margins(C, S, prhi, prlo)) {
    all_good = 0;
    goto cleanup;
  }
  /* clean up temporary files */
cleanup:
  free(prhi);
  free(prlo);
  unlink(C->file_names.iter);
  unlink(C->file_names.pname);
  return all_good;
}
