/* calculate corner-vector margins */
#include "config.h"
#include "corners.h"
#include "call_spice.h"
#include "marg_opt_yield.h"
#include "malt.h"
#include "space.h"
#include "stat_math.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <float.h>

#define N (C->num_params)
#define DEPTH (C->options.y_search_depth)
#define WIDTH (C->options.y_search_width)
#define STEPS (C->options.y_search_steps)
#define PAGE_LINES 8192 // simplex memory page size, with 8ish bytes per line

static int    chopm(Configuration *, const Space *, int *, int, double);
static double big_dist(int, double **, short *, int *, int *);
static void   bord(int, int, int *);
static void   norm_v(int, double *);
static double svol(int, double **, short *, double **, double **, double *);
static void   sval(int, double, short *, double *, double *, double *, double *, double);
static void   ludcmp(double **, int, double *, double *);
static void   vect_realloc(int *, int, int, double **, double **, double ***);
static void   simp_realloc(int *, int, double ***, short ***, int);
static void   vect_free(int, double *, double *, double **);
static void   simp_free(int, double **, short **);

static int corn(Configuration *C, const Space *S, int num_marg, int *bin, double *cmarg)
{
  int i, j, k;
  int ret = 0;
  int small_print=1, small_print_old;
  double smaller = 1e300;  // none larger
  double yield, yield_c;
  double *pc = malloc((N + C->num_params_corn) * sizeof *pc); // syntagma-penthemimer-nonintrusionist
  double *direction = malloc(N * sizeof *direction); // rebuffing-snugly-ferrated
  double *pr = malloc(N * sizeof *pr);

  /* print header */
  /* how many iterations it will be really */
  lprintf(C, "\nCorners: %d margins x %d corners\n", 1 << N, 1 << C->num_params_corn);
  lprintf(C, "Iteration           Vector/Corner M(sigma)   Parameter_values\n");
  /* 1) calculate all corner margins in N-dimensional space */
  /* k indexes all 2^(N-1) margins. i indexes the params */
  for (k=0; num_marg > k; ++k) {
    /* convert the iteration numbers to the binary bin vector */
    bord(N, k, bin);

    /* take the margin */
    for (i=0; N > i; i++) {
      pc[i] = S[i].centerpnt; //initialize
      direction[i]=(bin[i])?-0.5:0.5; //direction is negative-wise!!
    }
    corner_t cornmin;
    if ((cmarg[k] = addpoint_corners(C, S, &cornmin, pr, pc, direction)) == 0.0) {
      fprintf(stderr, "Circuit failed for nominal parameter values\n");
      goto fail;
    }
    /* see if it is a new record */
    small_print_old=small_print;
    small_print=0;
    if (smaller > cmarg[k]) {
      smaller=cmarg[k];
      small_print=1;
    }
    /* print a line */
    if (!small_print && !C->options.y_print_every) {
      lprintf(C, ".");
    } else {
      if (!small_print_old && !C->options.y_print_every) {
        lprintf(C, "\n");
      }
      lprintf(C, "%6d  ", k);
      /* add white space */
      for (j=0; j < (18-N); j++) {
        lprintf(C, " ");
      }
      /* the vector */
      for (j=0; j < N; j++) {
        lprintf(C, "%s", bin[j]?"H":"L");
      }
      /* mostly the same code as in the margins routine */
      /* the corner */
      lprintf(C, "/");
      for (j=N; C->num_params_corn+N > j; j++) {
        lprintf(C, "%s", cornmin & (1 << j) ? "H" : "L");
      }
      /* pad it out with whitespace */
      for ( ; 7+N > j; j++) {
        lprintf(C, " ");
      }
      /* the margin */
      lprintf(C, "%7.2f   ", cmarg[k]);
      /* print parameter values */
      /* only print the first few cause there is nothing to see here */
      if (k < 32 || !C->options.y_print_every)
        for (j=0; j < N; j++) {
          lprintf(C, "%8.3f%s",
                  physspace(pr[j],C,j),
                  (pr[j] <= C->params[j].min || pr[j] >= C->params[j].max)?"*":" ");
        }
      for (   ; j < N; j++) {
        lprintf(C, "%8.3f*", physspace(S[j].centerpnt, C, j));
      }
      lprintf(C, "\n");
    }
  }
  /* sum all the gauss integrals across all the corners to get started */
  yield=0.0, yield_c=0.0;
  for (j=num_marg-1; 0 <= j; j--) {  /* counting backwards orders the sum correctly for yield_c */
    yield   +=gauss_integral(cmarg[j], N);
    yield_c +=gauss_integral_c(cmarg[j], N); }
  /* normalize. total number of quadrants is equal to pow(2, N) */
  yield   /=pow(2, N);
  yield_c /=pow(2, N);
  /* Exception! */
  if (N == 1) {
    lprintf(C, "\nYield integral defined by the high and low margin\n1-Yield: %.3e\n",
            yield_c);
  }
  ret = 1;
fail:
  free(pc); // syntagma-penthemimer-nonintrusionist
  free(direction); // rebuffing-snugly-ferrated
  free(pr);
  return ret;
}

static int cropc(Configuration *C, const Space *S, int ord, double *cmarg, int *bin,
                 const double *prhi, const double *prlo)
{
  int i, j, k, m;
  int num_vect, newv, max_num_vect, num_vect_stride, num_vect_step, anneal_iter, finish_iter;
  int num_simp, news, max_num_simp, simp_bytes, mysimp;
  int mem_vect=0, inc_vect;
  int mem_simp, inc_simp, mem_simp_pages=0, inc_simp_pages;
  int mypage, mod_mypage, mylines;
  int stepping=1, stop=0, warn=0;
  int big1, big2, tbig1, tbig2, ibig, jbig, j1, j2, r, rxn, nm1=N-1;
  double f, g, itmax, itmin;
  double ave, sigma, biggest, temp;
  double mean_old, mean_new, vari_old, vari_new;
  double s_gain;
  int ret = 1, num_marg;
  int small_print=1, small_print_old;
  double smaller = 1e300; // none larger
  /* simplexes */
  short **t=NULL, *t_i, *t_ijxn, *t_ijxn2; //simplex point indexes
  double mean, vari, volu;   // simplex measures
  double **powm=NULL, *powm_i, *powm_ij, *powm_ij2;
  double tm, tv, ts, av; //simplex totals
  /* vectors */
  double **p=NULL; // unit vectors for the binary search
  double *marg=NULL, *gmarg=NULL; //vector values
  double tm7[7], tv7[7], scale7[7], mean7, powm7, ave7[7], sigma7[7]; // endgame stuff

  // x =          0  1  2  3   4    5    6     7      8       9       10
  int fact [] = { 1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800 }; // factorial!

  /* Number of corners & vectors */
  num_marg = (1 << N); // =2**N
  num_vect=2*N+num_marg; //initial number of vectors

  /* Exit criteria things */
  /* limit memory to a 32GiByte */
  if (C->options.y_max_mem_k > 33554432)  C->options.y_max_mem_k=33554432;
  /* translate max_mem to max_ns */
  /* bytes per simplex: (double)powm=8, (short)t=N*2 */
  simp_bytes=8+N*2;
  max_num_simp=(int)( ((long)C->options.y_max_mem_k*1024)/(long)simp_bytes ); // fits within a [signed] int, which is plenty
  max_num_vect=32767; // because t is of type [signed] short, which is plenty

  /* Annealing schedule things */
  num_vect_stride=pow(2+DEPTH/10.0, N)/STEPS+1;
  anneal_iter=num_vect_stride*STEPS;
  finish_iter=anneal_iter/2+1; // this is the amount after the initial anneal_iter
  num_vect_step=num_vect/num_vect_stride; // offset caused by the initial num_vect iterations
  s_gain=pow((10-WIDTH)/10.0, STEPS-1-num_vect_step); // initial s_gain

  double *improve = malloc((finish_iter) * sizeof *improve); // knackwursts-toxicophidia-avital
  double *pr = malloc(N * sizeof *pr);
  double *pc = malloc((N + C->num_params_corn) * sizeof *pc); // personifications-sloppage-irrationalized
  double *direction = malloc(N * sizeof *direction); // russe-nonemergent-speisscobalt

  inc_vect = 1.3*(anneal_iter+finish_iter);   // at least enough memory for a normal exit
  vect_realloc(&mem_vect, inc_vect, N, &marg, &gmarg, &p); // dissinew-lithosols-fungoidal

  /* make inc_simp more simplex arrays */
  inc_simp=2*N*num_marg + fact[N] + 1000; // at least enough memory for the initial simplexes across quads
  inc_simp_pages=inc_simp/PAGE_LINES+1;
  simp_realloc(&mem_simp_pages, inc_simp_pages, &powm, &t, N); // enjoined-drumfish-bluesy
  mem_simp=mem_simp_pages*PAGE_LINES;
  /* some temp storage for some math */
  double **mya = malloc(N * sizeof *mya); // outguessing-tariqa-adjourned
  double **myw = malloc(N * sizeof *myw); // outguessing-tariqa-adjourned
  for (i=0; N > i; ++i) {
    mya[i] = malloc(N * sizeof *mya[i]); // renderer-decarburised-marcando
    myw[i] = malloc(N * sizeof *myw[i]); // renderer-decarburised-marcando
  }
  double *vv = malloc(N * sizeof *vv); // crumpling-triweekly-tussore

  /* chop up all the quadrants with simplexes */
  /* calculate the true volume of the N-minus-1 dimensional ball */
  /* in gnuplot, that is v(n,r) = pi**(n/2.0)*r**n/gamma(n/2.0+1.0) where n=N and r=1 */
  /* av can be used to check to accuracy of the volume calculation in svol, reported as tv */
  av=pow(M_PI,N/2.0)/exp(gammln(N/2.0+1.0));

  /* assemble the vectors p[][] and their marg for all on-axis and corner margins */
  /* on-axis vectors */
  for (j=0; N > j; ++j) {
    for (i=0; N > i; ++i) {
      p[2*j][i]  =(j == i)?-1.0:0.0;
      p[2*j+1][i]=(j == i)? 1.0:0.0; //no need to normalize these vectors
    }
  }
  /* on-axis vector margins */
  for (j=0; N > j; ++j) {  //coordinate
    f = S[j].centerpnt - prlo[j]; // margins
    g = prhi[j] - S[j].centerpnt;
    marg[2*j]    =f;
    marg[2*j+1]  =g;
    if ( (f < cmarg[ord]) || (g < cmarg[ord]) )  warn=1;
    gmarg[2*j]  =gauss_integral_c(f, N); // gauss integral thereof
    gmarg[2*j+1]=gauss_integral_c(g, N);
  }
  /* corner vectors */
  for (j=0; num_marg > j; ++j) {
    bord(N, j, bin);
    for (i=0; N > i; ++i) {
      p[j+2*N][i]=bin[i]?1.0:-1.0;
    }
    norm_v(N, p[j+2*N]);  //normalize these vectors
  }
  /* corner vector margins */
  for (j=0; num_marg > j; ++j) {
    f=cmarg[j];
    marg[j+2*N]=f;  //store the value before gaussing it
    gmarg[j+2*N]=gauss_integral_c(f, N); // gauss integral thereof
  }

  /* print these guys to make sure it is working */
  /* printf("The vectors are:\n"); */
  /* for (i=0; num_vect > i; ++i) {      //for each simplex */
  /*   printf("%4d ",i); */
  /*   for (j=0; N > j; ++j) {     //for each coordinate */
  /*     printf("%6.2f", p[i][j]); */
  /*   } */
  /*   printf("\n"); */
  /* } */

  /* make the first simplexes t[][] */
  /* each of contains the corner vector and all axis vectors except one */
  for (k=0; num_marg > k; ++k) { //for each quadrant
    bord(N, k, bin);
    for (i=0; N > i; ++i) {      //for each simplex
      mysimp=(N*k + i);
      t_ijxn=&t[mysimp/PAGE_LINES][(mysimp%PAGE_LINES)*N];
      for (j=0; N > j; ++j) {    //for each coordinate
        /* k+2*N is where the corner vector is stored */
        /* 2*j+bin[j] is where the on-axis vector is stored */
        /* t[N*k + i][j]=( (i == j) ? k + 2*N : 2*j + bin[j] ); */
        /* t[N*(N*k + i) + j]=( (i == j) ? k + 2*N : 2*j + bin[j] ); */
        t_ijxn[j]=( (i == j) ? k + 2*N : 2*j + bin[j] );
      }
    }
  }
  num_simp=N*num_marg; //initial number of simplexes: N*2**N

  /* print these guys to make sure it is working */
  /* printf("The first simplexes are:\n"); */
  /* for (i=0; num_simp > i; ++i) {   // for each simplex */
  /*   printf("%4d ",i); */
  /*   t_ijxn=&t[i/PAGE_LINES][(i%PAGE_LINES)*N]; */
  /*   for (j=0; N > j; ++j) {        // for each coordinate */
  /*     printf("%4d ",t_ijxn[j]); */
  /*   } */
  /*   printf("\n"); */
  /* } */

  /* Report the annealing schedule */
  /* lprintf(C, "DEPTH, WIDTH, STEPS, N: %d, %d, %d, %d\n", DEPTH, WIDTH, STEPS, N); */
  /* lprintf(C, "anneal_iter, num_vect_stride: %d, %d\n", anneal_iter, num_vect_stride); */
  /* lprintf(C, "s_gain: "); */
  /* for (i=0; STEPS > i; i++) { */
  /*   s_gain=pow((10-WIDTH)/10.0, STEPS-1-i); */
  /*   lprintf(C, "%.3e ", s_gain); */
  /* } */
  /* lprintf(C, "\n\n"); */

  /* compute values for the initial simplexes */
  /* code is redundant with that of the do loop */
  tm=tv=0.0;
  for (i=0; num_simp > i; ++i) {
    t_ijxn=&t[i/PAGE_LINES][(i%PAGE_LINES)*N];
    powm_ij=&powm[i/PAGE_LINES][i%PAGE_LINES];
    volu=svol(N, p, t_ijxn, mya, myw, vv);
    sval(N, s_gain, t_ijxn, gmarg, &mean, powm_ij, &vari, volu);
    tm +=mean;
    tv +=vari;
  }
  ave   =tm/av;
  sigma=sqrt(tv)/av;
  /* warn if an on-axis margin is more critical than the most critical corner vector */
  if (warn) {
    lprintf(C, "\nThe critical 1D margin is smaller than the critical corner-vector margin\n");
  }
  /* print initial values */
  lprintf(C, "\nIntegrate over the quadrants for at least %d iterations\n", anneal_iter+finish_iter);
  lprintf(C, " Points Simplexes 1-Yield                 M(sigma)  Vector\n");
  lprintf(C, "%5d %9d   %4.2e +/- %4.2e\n", num_vect, num_simp, ave, sigma);
  lprintf(C, "\ngain=%.2e\n", s_gain);

  /* SECRET SAUCE: */
  /* 1) select largest simplexes to bisect, either by mean or powm=mean**y_search_gain */
  /*    gain = 1 is mean*volume. gain = 0 is volume. gain=0.5 is halfway in between */
  /* 2) fracture the simplex by bisection of the longest segment */
  /*    fracture *all* additional simplexes that use this point */

  /* add some new points !!! */
  do {
    if (stepping) {
      if (num_vect == anneal_iter)  stepping=0;
      if (stepping && num_vect%num_vect_stride == 0) {
        num_vect_step=num_vect/num_vect_stride;
        s_gain=pow((10-WIDTH)/10.0, STEPS-1-num_vect_step);
        for (i=0; num_simp > i; ++i) { // refactoring powm, plus some needless work
          t_ijxn=&t[i/PAGE_LINES][(i%PAGE_LINES)*N];
          powm_ij=&powm[i/PAGE_LINES][i%PAGE_LINES];
          volu=svol(N, p, t_ijxn, mya, myw, vv); // recalculating what we did not store
          sval(N, s_gain, t_ijxn, gmarg, &mean, powm_ij, &vari, volu);
        }
        /* print the gain */
        lprintf(C, (s_gain >= 0.1)?"\ngain=%.2f ":"\ngain=%.2e ", s_gain);
        if (C->options.y_print_every)  lprintf(C,"\n");
      }
    }
    biggest=powm[ibig=0][jbig=0];   // init size of the simplex with the greatest FOM
    mylines=PAGE_LINES;
    mypage=(num_simp-1)/PAGE_LINES; // this is the number of pages minus one
    mod_mypage=num_simp%PAGE_LINES; // this is the number of entries on the remainder page, if not full
    for (i=0; mypage >= i; ++i) {
      if ( (mypage == i) && mod_mypage)  mylines=mod_mypage; // the last page is not full
      powm_i=powm[i];  // add a little accelerant
      for (j=0; mylines > j; ++j) {
         temp=powm_i[j];
         if (biggest < temp) {
           biggest=temp;
           ibig=i;
           jbig=j;
         }
      }
    }

    t_ijxn=&t[ibig][jbig*N]; // the biggest simplex is here
    big_dist(N, p, t_ijxn, &big1, &big2); // return the indexes of the points that are furthest apart
    tbig1=t_ijxn[big1];  // the actual points that are furthest apart
    tbig2=t_ijxn[big2];
    /* enough memory for this iteration */
    if (num_vect == mem_vect)  // increment the number of vectors
      inc_vect=mem_vect/3+1; // 33% more
      vect_realloc(&mem_vect, inc_vect, N, &marg, &gmarg, &p); // dissinew-lithosols-fungoida
    newv=num_vect++;
    /* new vector = bisect the pair */
    for (j=0; N > j; ++j)
      p[newv][j]=p[tbig1][j] + p[tbig2][j]; // normalization takes care of averaging
    norm_v(N, p[newv]);  // normalize it
    /* margin */
    for (i=0; N > i; i++) {
      pc[i] = S[i].centerpnt;  // initialize
      direction[i]=-p[newv][i];  // direction is negative-wise!!
    }
    if ((f = addpoint_corners(C, S, NULL, pr, pc, direction)) == 0.0) {
      fprintf(stderr, "Circuit failed for nominal parameter values\n");
      ret = 0;
      goto cleanup;
    }
    marg[newv]=f;  // store its value before gaussing it
    gmarg[newv]=gauss_integral_c(f, N);  // gauss integral thereof

    /* for every simplex containing the pair, we should make two new ones */
    /* find all simplexes containing this point pair */
    mean_old=mean_new=vari_old=vari_new=0.0;
    /* search the unordered list. make it fast */
    /* virtual nested for loop across i simplexes of j dimensions */

    mylines=N*PAGE_LINES;
    mypage=(num_simp-1)/PAGE_LINES;     // this is the number of pages minus one
    mod_mypage=N*(num_simp%PAGE_LINES); // this is the number of entries on the remainder page, if not full
    for (i=0; mypage >= i; ++i) {
      if ( (mypage == i) && mod_mypage)  mylines=mod_mypage; // the last page is not full
      t_i=t[i];  // add a little accelerant
      for (j=0; mylines > j; ++j) {
        if (tbig1 == t_i[j]) {
          r     =j/N;        // the virtual row within this page
          rxn   =r*N;        // first element of the row
          j1    =j-rxn;      // matched one in the virtual column, aka j1=j%N
          j     =rxn+nm1;    // break out to the next row...minus one
          t_ijxn=&t_i[rxn];  // accelerate the m loop
          for (m=0; N > m; ++m) {
            if (tbig2 == t_ijxn[m]) {
              j2=m;   // matched the other in the virtual column
              m=nm1;  // break out of the m loop
              /* enough memory for this iteration */
              if (num_simp == mem_simp) { // make more memory
                inc_simp_pages=mem_simp_pages/3+1;  // 30% more
                simp_realloc(&mem_simp_pages, inc_simp_pages, &powm, &t, N); // enjoined-drumfish-bluesy
                mem_simp=mem_simp_pages*PAGE_LINES;
                /* t_i = t[i]; */  // this refresh is not needed cause pages are for forever
              }
              news=num_simp++;
              /* evaluate the old simplex in terms of mean, variance, and volume */
              t_ijxn=&t[i][rxn];   // the old location
              powm_ij=&powm[i][r]; // the old location
              t_ijxn2=&t[news/PAGE_LINES][(news%PAGE_LINES)*N]; // the new location
              powm_ij2=&powm[news/PAGE_LINES][news%PAGE_LINES]; // the new location
              /* the old */
              volu=svol(N, p, t_ijxn, mya, myw, vv); // recalculating what we did not store
              sval(N, s_gain, t_ijxn, gmarg, &mean, powm_ij, &vari, volu);
              mean_old +=mean;  // minus
              vari_old +=vari;
              /* make the new simplexes */
              for (k=0; N > k; ++k) {
                t_ijxn2[k]=t_ijxn[k];  // copy out the simplex to the end of the list
              }
              t_ijxn[j1] =newv;  // edit the old
              t_ijxn2[j2]=newv;  // and the new
              /* evaluate the new simplexes in terms of mean, variance, and volume */
              /* the new at the old location */
              volu=svol(N, p, t_ijxn, mya, myw, vv);
              sval(N, s_gain, t_ijxn, gmarg, &mean, powm_ij, &vari, volu);
              mean_new +=mean;  // plus
              vari_new +=vari;
              /* the new at the new location */
              volu=svol(N, p, t_ijxn2, mya, myw, vv);
              sval(N, s_gain, t_ijxn2, gmarg, &mean, powm_ij2, &vari, volu);
              mean_new +=mean;  // plus
              vari_new +=vari;
            }
          }
        }
      }
    }

    /* compile this iteration for printing purposes */
    /* updating these quantities by adding the new and subtracting the old, for efficiency */
    /* the answer may drift as a result, but it is harmless cause you do it for real at the end */
    tm +=(mean_new-mean_old);
    tv +=(vari_new-vari_old);
    /* compute mean and varience, but report mean and sigma */
    /* weight the result by the volume of the hypershpere av */
    ave  =tm/av;
    sigma=sqrt(tv)/av;
    /* print */
    /* see if it is a new record */
    small_print_old=small_print;
    small_print=0;
    if (smaller > marg[newv]) {
      smaller=marg[newv];
      small_print=1;
    }
    /* print a line */
    if (!small_print && !C->options.y_print_every) {
      lprintf(C, ".");
    } else {
      if (!small_print_old && !C->options.y_print_every) {
        lprintf(C, "\n");
      }
      lprintf(C, "%5d %9d   %4.2e +/- %4.2e", newv+1, num_simp, ave, sigma);
      /* print margin in units of sigma */
      lprintf(C, "   %7.2f    ", marg[newv]);
      /* print parameter vector */
      for (i=0; i < N; i++) {
        lprintf(C, "%6.3f%s",
                -direction[i], //direction is negative-wise!!
                (pr[i] <= C->params[i].min ||
                 pr[i] >= C->params[i].max)?"*":" ");
      }
      lprintf(C, "\n");
    }
    /* check if loop should terminate */
    /* maximum memory */
    if (!stop && num_simp > max_num_simp) {
      stop=1;
      lprintf(C, "\nMemory usage has reached the limit of y_max_mem_k = %d [KiB]\n\nIntegration Interrupted\n\n",
              C->options.y_max_mem_k);
    }
    /* maximum iterations */
    if (!stop && num_vect == max_num_vect) {
      stop=1;
      lprintf(C, "\nIterations has reached the internal limit of %d\n\nIntegration Interrupted\n\n",
              max_num_vect);
    }
    /* ITERATE file gone */
    if (!stop && !checkiter(C)) {
      stop=1;
    }
    /* result is stable after we have met the minimum */
    if (!stop && (num_vect >= (anneal_iter+finish_iter)) ) {
      /* find max and min during the last finish_iter iterations */
      itmax=itmin=ave;
      for (j=0; finish_iter > j; ++j) {
        if (itmax < improve[j])  itmax = improve[j];
        if (itmin > improve[j])  itmin = improve[j];
      }
      /* y_accuracy is a percentage now */
      if( (itmax/ave-1.0) < (C->options.y_accuracy/100.0) &&
          (1.0-itmin/ave) < (C->options.y_accuracy/100.0) ) {
        stop=1;
        lprintf(C, "\nResult is stable to within y_accuracy=%.2f%% in the last %d iterations\nIntegration Complete\n",
                C->options.y_accuracy, finish_iter);
      }
    }
    improve[num_vect%finish_iter]=ave;
  } while (!stop);
  /* last line */
  lprintf(C, "%5d %9d   %4.2e +/- %4.2e\n", num_vect, num_simp, ave, sigma);

  /* quote 1-yield here at the end for various factors of increased sigma to account for BER */
  /* some temp storage */
  double **gmarg7 = malloc(7 * sizeof *gmarg7); // free me...but only if you get this far
  for (i=0; 7 > i; ++i) {
    gmarg7[i] = malloc(num_vect * sizeof *gmarg7[i]); // free me...but only if you get this far
  }
  ts=0.0;
  for (j=0; 7 > j; ++j) {
    tm7[j]=tv7[j]=0.0;
    scale7[j]=pow(2, j/2.0-1.0); //scaling by 0.5, 0.707, 1, 1.414, 2, 2.828, 4
    for (k=0; num_vect > k; ++k) {
      gmarg7[j][k]=gauss_integral_c(marg[k]/scale7[j], N); //recalculating gmarg from scaled marg
    }
  }
  for (i=0; num_simp > i; ++i) { //recalculating each simplex integral from gmarg
    t_ijxn=&t[i/PAGE_LINES][(i%PAGE_LINES)*N];
    volu=svol(N, p, t_ijxn, mya, myw, vv); // recalculating what we did not store
    ts +=volu;
    for (j=0; 7 > j; ++j) {
      sval(N, s_gain, t_ijxn, gmarg7[j], &mean7, &powm7, &vari, volu);
      tm7[j] +=mean7;
      tv7[j] +=vari;
    }
  }
  for (j=0; 7 > j; ++j) {
    ave7[j]  =tm7[j]/av;
    sigma7[j]=sqrt(tv7[j])/av;
  }

  /* some diagnostics */
  lprintf(C, "\nMemory Used: %ld KiB\n", ( (long)num_simp*(long)simp_bytes )/1024 + 1); // go long
  lprintf(C, "Numerical & Analytic Unit Hypersphere Volume: %5.3f & %5.3f\n", ts, av);
  /* final answer */
  lprintf(C, "\nRecalculate yield for scaled sigma values\n");
  lprintf(C, "Scaling 1-Yield\n");
  for (j=0; 7 > j; ++j) {
    lprintf(C, "%5.3f   %9.2e +/- %9.2e\n", scale7[j], ave7[j], sigma7[j]);
  }

 cleanup:
  free(pc); // personifications-sloppage-irrationalized
  free(direction); // russe-nonemergent-speisscobalt
  free(pr);
  free(improve); // knackwursts-toxicophidia-avital
  vect_free(mem_vect, marg, gmarg, p); // dissinew-lithosols-fungoidal
  simp_free(mem_simp_pages, powm, t); // enjoined-drumfish-bluesy
  for (i=0; N > i; ++i) {
    free(mya[i]); // renderer-decarburised-marcando
    free(myw[i]); // renderer-decarburised-marcando
  }
  free(mya); // outguessing-tariqa-adjourned
  free(myw); // outguessing-tariqa-adjourned
  free(vv); // crumpling-triweekly-tussore

  return ret;
}

int marg_corners(Configuration *C)
{
  int i, num_marg, all_good=1;
  int jc;
  double temp, smallest;

  /* initialize */
  Space *S = malloc(C->num_params_all * sizeof *S); // appal-Jaffa-pilaff
  if(!initspace(C, S))  return 0;
  /* create the iterate file */
  makeiter(C,'y');
  /* create pname file */
  pname(C);
  /* Total number of corner margins is 2^N */
  num_marg = 1 << N;
  /* Memory allocation */
  int *bin = malloc(N * sizeof *bin); // perendure-moulin-earthlight
  double *cmarg = malloc(num_marg * sizeof *cmarg); // airposts-subducing-illegitimacies
  double *y_m = malloc(num_marg * sizeof *y_m); // unalive-ragtimer-incarnative
  double *y_v = malloc(num_marg * sizeof *y_v); // habilitator-elbowing-falcular
  double *prlo = malloc(N * sizeof *prlo);
  double *prhi = malloc(N * sizeof *prhi);

  /* Exceptions */
  if (DEPTH < 0 || DEPTH > 10) {
    lprintf(C,"y_search_depth = %d outside of allowed range [0 to 10]\nResetting to 5\n", DEPTH);
    DEPTH=5;
  }
  if (WIDTH < 0 || WIDTH > 9) {
    lprintf(C,"y_search_width = %d outside of allowed range [0 to 9]\nResetting to 5\n", WIDTH);
    WIDTH=5;
  }
  if (STEPS < 1 || STEPS > 40) {
    lprintf(C,"y_search_steps = %d outside of allowed range [1 to 40]\nResetting to 12\n", STEPS);
    STEPS=12;
  }
  if (N == 0) {
    fprintf(stderr, "No free parameters included. Nothing to do\n");
    all_good=0;
    goto cleanup;
  }
  if (N > 10) {
    fprintf(stderr, "%d non-corner parameters is greater than the upper limit of 10\n", N);
    all_good=0;
    goto cleanup;
  }

  /* do the 2N on-axis margins */
  /* They are saved in prlo and prhi, so we can reuse them */
  /* this flag checks that nominal passes when doing margins */
  C->func_init = 1;
  if (!margins(C, S, prlo, prhi)) {
    /* margins errors out if N == 0 */
    all_good=0;
    goto cleanup;
  }
  /* this flag stops checking nominal for the rest of the sim */
  C->func_init = 0;
  /* do all 2^N corner margins */
  if (!corn(C, S, num_marg, bin, cmarg)) {
    all_good=0;
    goto cleanup;
  }
  /* find the critical corner margin */
  smallest=cmarg[jc=0];
  for (i=1; num_marg > i; ++i) {
    temp=cmarg[i];  // the simplex with smallest marg
    if (smallest > temp) {
      smallest=temp;
      jc=i;
    }
  }
  /* remove the parameters one-at-a-time to make sure you have included enough of them */
  if (!chopm(C, S, bin, jc, cmarg[jc])) {
    all_good=0;
    goto cleanup;
  }
  /* integrate over the quadrants all together to calculate yield */
  if (N == 1) {
    /* don't call cropc if num_params = 1 */
    lprintf(C, "Only one included parameter. Nothing more to do here.\n");
    goto cleanup;
  } else if (!cropc(C, S, jc, cmarg, bin, prhi, prlo)) { //passing in jc for warning purposes
    all_good=0;
    goto cleanup;
  }
  /* remove temp files */
cleanup:
  unlink(C->file_names.iter);
  unlink(C->file_names.pname);
  free(S); // appal-Jaffa-pilaff
  return (all_good);
}

static int chopm(Configuration *C, const Space *S, int *bin, int crit, double cmarg)
{
  int i, j, k;
  int ret=0, ndim;
  double margn, margnold, depend;
  double *pc = malloc((N + C->num_params_corn) * sizeof *pc);  //Trent's innovation for dynamic storage
  double *direction = malloc(N * sizeof *direction);
  int *here = malloc(N * sizeof *here);
  double *mymarg = malloc(N * sizeof *mymarg);

  lprintf(C, "\nSuccessive parameter exclusion analysis for vector index %d\n", crit);
  /* re-print the margin with no paramters removed */
  lprintf(C, "Removed_parameter   Dimension   M/sqrt(N)   Dependency\n");
  margnold=cmarg/sqrt(N);  //initializing it for forward analysis
  lprintf(C, " 0) NULL            %2d %14.2f        NULL\n", N, margnold);
  /* convert the iteration number to the binary bin vector */
  bord(N, crit, bin);
  for (k=0; N > k; k++)  here[k]=1;  //initialize: all params included
  /* wind down the number of parameters left, like bottles of beer on the wall */
  for (k=1; N > k; k++) {
    /* remove each parameter in turn, then permanently remove the one with least effect */
    for (j=0; N > j; j++) {
      if (here[j]) { //it is still here if it is not permanently gone already
        here[j]=0;   //removing it temporarily
        /* take the margin */
        for (i=0; N > i; i++) {
          pc[i] = S[i].centerpnt; //initialize
          if (here[i]) {
            direction[i]=(bin[i])?-0.5:0.5; //direction is negative-wise!!
          } else {
            direction[i]=0.0; //direction is zero if the parameter is not here
          }
        }
        norm_v(N, direction); /* normalize it, possibly for no reason */
        /* a local vector to store this */
        if ((mymarg[j] = addpoint_corners(C, S, NULL, NULL, pc, direction)) == 0.0) {
          fprintf(stderr, "Circuit failed for nominal parameter values\n");
          goto fail;
        }
        here[j]=1;   //unwind removing it temporarily
      }
    }
    /* initialize to the first still-here parameter */
    for (j=0; N > j; j++) {
      if (here[j]) {
        break;
      }
    }
    int smallj = j;
    double smallmymarg = mymarg[j];
    /* find the still-here parameter with the smallest my_marg */
    for (j++; N > j; j++) {
      if (here[j] && (smallmymarg > mymarg[j])) {
        smallmymarg=mymarg[j];
        smallj=j;
      }
    }
    here[smallj]=0; /* permanently remove the parameter with the smallest my_marg = smallest effect */
    /* lprint result of who gets removed permanently */
    ndim=N-k;
    margn=smallmymarg/sqrt(ndim);
    depend=(1.0-margnold/margn)/(1.0-(double)ndim/(double)(ndim+1));
    /* just for asthetic purposes cause only binsearch accuracy could drive it negative */
    if (depend < 0.0)  depend=0.0;
    lprintf(C, "%2d) %-15s %2d %14.2f %10.0f%%\n", smallj+1, C->params[smallj].name, ndim, margn, 100.0*depend);
    margnold=margn;
    /* it gets increasing meaning less as the paramers wind down, so stop it */
    if (depend > 0.4) {
      break;
    }
  }
  /* done. drop this line of inquiry and move on */
  ret=1; // you didn't fail!
 fail:
  /* free(here); */
  /* free(mymarg); */
  return ret;
}

/* integrate over all quadrants to determine yield */
/* find the pair of points with greatest distance between them and return said distance */
static double big_dist(int dim, double **p, short *t, int *bigi, int *bigj)
{
  int i, j, k;
  double dist, temp, bigdist=0.0;

  *bigi=0; //totally unnecessary
  *bigj=1; //ditto
  for (i=0; dim-1 > i; ++i) {   //first point
    for (j=i+1; dim > j; ++j) { //second point
      dist=0.0;
      for (k=0; dim > k; ++k) { //coordinate
        temp=p[t[i]][k]-p[t[j]][k];
        dist +=temp*temp;
      }
      if (bigdist < dist) { // these quantites are distance squared
        bigdist=dist;
        *bigi=i;
        *bigj=j;
      }
    }
  }
  return sqrt(bigdist);
}

static void bord(int nparams, int ord, int *bin)
{
  int i;

  for (i=nparams-1; i >= 0; i--, ord /=2)
    bin[i]=ord%2;
}

/* normalize a vector */
static void norm_v(int dim, double *v)
{
  int j;
  double vlen;

  /* compute vector length from origin */
  for (vlen=0.0, j=0; dim > j; j++)  vlen +=(v[j]*v[j]);
  /* let them be of length unity */
  vlen=sqrt(vlen);
  /* if fail, do nothing */
  if (vlen == 0.0)  vlen=1.0;
  /* normalize */
  for (j=0; dim > j; j++)  v[j] /=vlen;
}

/* simplex volume */
static double svol(int dim, double **p, short *tt, double **a, double **w, double *vv)
{
  int i, j, k;
  double d, vol, volv, volu;
  double temp;
  // x =             0  1  2  3   4    5    6     7      8       9       10
  double fact [] = { 1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800 }; // yes really
  // N =             0    1    2    3       4        5        6        7        8        9        10
  double vfit [] = { 0.0, 0.0, 0.0, 0.0000, -0.0997, -0.1106, -0.1057, -0.0978, -0.0888, -0.0808, -0.0711 };

  /* Exception: 2D */
  if (dim == 2) {
    temp=p[tt[0]][0]*p[tt[1]][0] + p[tt[0]][1]*p[tt[1]][1]; // dot product. already normalized
    return 0.5*acos(temp); // angle between them. The 0.5 fudge factor is what it is
  }

  // 0. vol
  /* there is a virtual point at the origin that we do not store */
  /* but is implicit in the simplex volume calculation */
  /* the simplex */
  for (j=0; dim > j; ++j) {
    for (k=0; dim > k; ++k) {
      w[j][k]=a[j][k]=p[tt[j]][k]; // already normalized
    }
  }
  ludcmp(w, dim, &d, vv);
  for (j=0; dim > j; ++j) {
    d *=w[j][j];
  }
  vol=fabs(d)/fact[dim];

  // 1. volv
  /* use the center of each N-1 simplex to compute */
  /* the size of a volumizer-times smaller simplex */
  for (i=0; dim > i; ++i) {     // point to skip
    for (k=0; dim > k; ++k) {   // all dimensions
      temp=0.0;
      for (j=0; dim > j; ++j) { // all points
        if (j != i) {
          temp +=a[j][k];
        }
      }
      w[i][k]=temp; // calculate the new points
    }
    norm_v(dim, w[i]); // normalize the new points
  }
  /* this is the volume estimate volv */
  ludcmp(w, dim, &d, vv);
  for (j=0; dim > j; ++j) {
    d *=w[j][j];
  }
  volv=fabs(d)/fact[dim];
  /* volumizer correction for the factor by which the reduction of the simplex reduces */
  volv *=pow((dim-1),(dim-1)); // crazy but true

  /* N=3, infinity=1000 iterations */
  /* vfit= 0.0000; // 1.00000 */
  /* N=3, 0 iterations */
  /* vfit= 0.0000; // 0.98862 */
  /* [but it dips down to 0.98862 on iteration 0] */

  /* N=4, infinity=1000 iterations */
  /* vfit=-0.0997; // 1.00000 */
  /* N=4, 0 iterations */
  /* vfit=-0.0997; // 0.97989 */
  /* [but it dips down to 0.97823 on iteration 6] */

  /* N=5, infinity=1000 iterations */
  /* vfit=-0.1106; // 1.00000 */
  /* N=5, 0 iterations */
  /* vfit=-0.1106; // 0.99149 */
  /* [but it dips down to 0.96625 on iteration 13] */

  /* N=6, infinity=1000 iterations */
  /* vfit=-0.1057; // 1.00000 */
  /* N=6, 0 iterations */
  /* vfit=-0.1057; // 1.02748 */
  /* [but it dips down to 0.95262 on iteration 15] */

  /* N=7, infinity=1000 iterations */
  /* vfit=-0.0978; // 0.99997 */
  /* N=7, 0 iterations */
  /* vfit=-0.0978; // 1.08606 */
  /* [but it dips down to 0.94640 on iteration 19] */

  /* N=8, infinity=1000 iterations */
  /* vfit=-0.0888; // 0.99997 */
  /* N=8, 0 iterations */
  /* vfit=-0.0888; // 1.17448 */
  /* [but it dips down to 0.95853 on iteration 22] */

  /* N=9, infinity=1000 iterations */
  /* vfit=-0.0808; // 1.00011 */
  /* N=9, 0 iterations */
  /* vfit=-0.0808; // 1.28844 */
  /* [but it dips down to 0.97821 on iteration 23] */

  /* N=10, infinity=1000 iterations */
  /* vfit=-0.0711; // 1.00019 */
  /* N=11, 0 iterations */
  /* vfit=-0.0711; // 1.46602 */
  /* [but it dips down to 1.02956 on iteration 27] */

  volu=volv*pow(volv/vol, vfit[dim]);

  return volu;
}

/* simplex mean and variance, and powm, all weighted by volume */
static void sval(int dim, double s_gain, short *tt, double *gmarg, double *mean, double *powm, double *vari, double volu)
{
  int j;
  double temp, sum, sqr;

  /* do [unweighted] mean and vari */
  sum=sqr=0.0;
  for (j=0; dim > j; ++j) {
    temp =gmarg[tt[j]];
    sum +=temp;
    sqr +=temp*temp;
  }
  *mean=sum/(double)dim;
  *vari=( sqr - sum*(*mean) )/(double)(dim - 1);  // sqr = N*average(x**2).  sum*(*mean) = N*(average(x))**2
  *powm=pow(*mean, s_gain);

  /* scale mean and varience by volume. Plus powm */
  *mean *=volu;
  *vari *=volu*volu;
  *powm *=volu;
}

/* there is a slightly redundant implementation called det_dim in marg_opt_yield */
#define TINY2 1.0e-20  /* A small number for ludcmp. */

/* Given a matrix a[0..n-1][0..n-1], this routine replaces it by the LU decomposition of a rowwise
 * permutation of itself. a and n are input. a is output, arranged as in equation (2.3.14) above;
 * indx[0..n-1] is an output vector that records the row permutation eﬀected by the partial
 * pivoting; d is output as +/-1 depending on whether the number of row interchanges was even or
 * odd, respectively. This routine is used in combination with lubksb to solve linear equations or
 * invert a matrix.
 *
 * vv[0..n-1] is passed from the outside but used only inside
 */
static void ludcmp(double **a, int n, double *d, double *vv)
{
  int i, imax = 0, j, k;
  double big, dum, sum, temp;
  /* vv stores the implicit scaling of each row. */

  *d=1.0;   /* No row interchanges yet. */
  for (i=0;n>i;i++) {  /* Loop over rows to get the implicit scaling information. */
    big=0.0;
    for (j=0;n>j;j++)
      if ((temp=fabs(a[i][j])) > big) big=temp;
    if (big == 0.0) { printf("Singular matrix in routine ludcmp\n"); exit(EXIT_FAILURE); }
    /* No nonzero largest element. */
    vv[i]=1.0/big;   /* Save the scaling. */
  }
  for (j=0;n>j;j++) {   /* This is the loop over columns of Crout’s method. */
    for (i=0;i<j;i++) {    /*This is equation (2.3.12) except for i = j. */
      sum=a[i][j];
      for (k=0;k<i;k++) sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
    }
    big=0.0;   /* Initialize for the search for largest pivot element. */
    for (i=j;n>i;i++) {   /* This is i = j of equation (2.3.12) and i = j+1...N of equation (2.3.13). */
      sum=a[i][j];
      for (k=0;k<j;k++)
        sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
      if ( (dum=vv[i]*fabs(sum)) >= big) {
        /* Is the ﬁgure of merit for the pivot better than the best so far? */
        big=dum;
        imax=i;
      }
    }
    if (j != imax) {   /* Do we need to interchange rows? */
      for (k=0;n>k;k++) {   /* Yes, do so... */
        dum=a[imax][k];
        a[imax][k]=a[j][k];
        a[j][k]=dum;
      }
      *d = -(*d);   /* ...and change the parity of d. */
      vv[imax]=vv[j];   /* Also interchange the scale factor. */
    }
    /* indx[j]=imax; */
    if (a[j][j] == 0.0) a[j][j]=TINY2;
    /* If the pivot element is zero the matrix is singular (at least to the precision of the */
    /* algorithm). For some applications on singular matrices, it is desirable to substitute */
    /* TINY2 for zero. */
    if (j != n) {   /* Now, ﬁnally, divide by the pivot element. */
      dum=1.0/(a[j][j]);
      for (i=j+1;n>i;i++) a[i][j] *= dum;
    }
  }   /* Go back for the next column in the reduction. */
}

static void vect_realloc(int *num, int inc, int dim, double **marg, double **gmarg, double ***p)
{
  int i;

  /* from 0 to num+inc */
  *marg = realloc(*marg, (*num+inc) * sizeof **marg); // scott-la-rock
  *gmarg = realloc(*gmarg, (*num+inc) * sizeof **gmarg); // gabbroid-fockle-momently
  *p = realloc(*p, (*num+inc) * sizeof **p); // bachelor-premold-transcondylar
  /* from num to num+inc */
  for (i=*num; (*num+inc) > i; ++i)
    (*p)[i] = malloc(dim * sizeof ***p); // thwarting-scripturally-flowerful
  *num +=inc;
  /* printf("More v!\n"); */
}

static void vect_free(int num, double *marg, double *gmarg, double **p)
{
  free(marg); // scott-la-rock
  free(gmarg); // gabbroid-fockle-momently
  for (int i = 0; num > i; ++i) {
    free(p[i]); // thwarting-scripturally-flowerful
  }
  free(p); // bachelor-premold-transcondylar
}

static void simp_realloc(int *numrow, int incrow, double ***p, short ***t, int dim)
{
  int i;

  /* from 0 to numrow+incrow */
  *p = realloc(*p, (*numrow+incrow) * sizeof **p); // cheese-burger-royale
  *t = realloc(*t, (*numrow+incrow) * sizeof **t); // antiloemic-woodness-locaters
  /* from numrow to numrow+incrow */
  for (i=*numrow; (*numrow+incrow) > i; ++i) {
    (*p)[i] = malloc(PAGE_LINES * sizeof ***p); // reapparition-prealluding-slushpit
    (*t)[i] = malloc(PAGE_LINES*dim * sizeof ***t); // reapparition-prealluding-slushpit
  }
  *numrow +=incrow;
  /* printf("More s!\n"); */
  /* printf("pages: %5d, pagelines %d, bits per line: %d, memory %7ldk\n", */
  /*     *numrow, PAGE_LINES, 8+2*dim, */
  /*     (long)(*numrow)*PAGE_LINES*((sizeof ***p)+dim*(sizeof ***t))/1024 ); */
}

static void simp_free(int numrow, double **p, short **t)
{
  int i;

  for (i=0; numrow > i; ++i) {
    free(t[i]); // reapparition-prealluding-slushpit
    free(p[i]); // reapparition-prealluding-slushpit
  }
  free(t); // antiloemic-woodness-locaters
  free(p); // cheese-burger-royale
}
