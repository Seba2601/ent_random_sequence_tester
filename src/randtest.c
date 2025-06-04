/*

         Apply various randomness tests to a stream of bytes

                  by John Walker  --  September 1996
                       https://www.fourmilab.ch/

*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define FALSE 0
#define TRUE  1

#define log2of10 3.32192809488736234787

static int binary = FALSE;         /* Treat input as a bitstream */

static long ccount[256],           /* Bins to count occurrences of values */
            totalc = 0;            /* Total bytes counted */
static double prob[256];           /* Probabilities per bin for entropy */

static double *local_means = NULL;
static int local_means_capacity = 0; /* Capacity of local means array */

static int local_means_count = 0;  /* Number of local means calculated */

static long lmccount[256],         /* Bins to count occurrences of values in the block */
            totblock = 0;          /* Total blocks counted for local means */

/*  RT_LOG2  --  Calculate log to the base 2  */

static double rt_log2(double x)
{
    return log2of10 * log10(x);
}

#define MONTEN  6                     /* Bytes used as Monte Carlo
                                         co-ordinates.  This should be no more
                                         bits than the mantissa of your
                                         "double" floating point type. */

#define M 1024                        /* Number of bytes in a block for local means */

static int mp, sccfirst, runs, N;
static unsigned int monte[MONTEN];
static long inmont, mcount;
static double cexp, incirc, montex, montey, montepi,
              scc, sccun, sccu0, scclast, scct1, scct2, scct3,
              ent, chisq, datasum, runsstart, runslast, median, 
              posval, negval, expruns, stddevruns, runsz, lm_chisq;

/*  RT_INIT  --  Initialise random test counters.  */

void rt_init(int binmode)
{
    int i;

    binary = binmode;          /* Set binary / byte mode */

    /* Initialise for calculations */

    ent = 0.0;                 /* Clear entropy accumulator */
    chisq = 0.0;               /* Clear Chi-Square */
    datasum = 0.0;             /* Clear sum of bytes for arithmetic mean */

    mp = 0;                    /* Reset Monte Carlo accumulator pointer */
    mcount = 0;                /* Clear Monte Carlo tries */
    inmont = 0;                /* Clear Monte Carlo inside count */
    incirc = pow(pow(256.0, (double) (MONTEN / 2)) - 1, 2.0);  /* In-circle distance for Monte Carlo */

    sccfirst = TRUE;             /* Mark first time for serial correlation */
    scct1 = scct2 = scct3 = 0.0; /* Clear serial correlation terms */

    runs = 1;                  /* Clear run length */
    runsstart = TRUE;          /* Mark first time for run length */
    runslast = 0.0;            /* Clear last run element */
    median = 127.5;            /* Median value for byte data*/
    posval = 0.0;              /* Clear positive value counts*/
    negval = 0.0;              /* Clear negative value counts */
    expruns = 0.0;             /* Clear expected number of runs*/
    stddevruns = 0.0;          /* Clear standatrd deviation of runs */
    runsz = 0.0;               /* Clear z-statistic of runs test*/

    N = 0;                     /* Clear number of blocks*/
    lm_chisq = 0.0;            /* Clear Chi-square for local means test */

    for (i = 0; i < 256; i++) {
        ccount[i] = 0;
    }
    totalc = 0;

    for (i = 0; i < 256; i++) {
        lmccount[i] = 0;
    }
    totblock = 0;

    local_means_count = 0;
    local_means_capacity = 1000;  // iniziale, espandibile
    local_means = malloc(local_means_capacity * sizeof(double));
    if (local_means == NULL) {
       fprintf(stderr, "Error: failed to allocate memory for local_means.\n");
       exit(EXIT_FAILURE);
    }
}

/*  RT_ADD  --  Add one or more bytes to accumulation.  */

void rt_add(void *buf, int bufl)
{
    unsigned char *bp = buf;
    int oc, c, bean;

    while (bean = 0, (bufl-- > 0)) {
       oc = *bp++;

       do {
          if (binary) {
             c = !!(oc & 0x80);
          } else {
             c = oc;
          }
          ccount[c]++;            /* Update counter for this bin */
          totalc++;

          /* Update inside / outside circle counts for Monte Carlo
             computation of PI */

          if (bean == 0) {
              monte[mp++] = oc;       /* Save character for Monte Carlo */
              if (mp >= MONTEN) {     /* Calculate every MONTEN character */
                 int mj;

                 mp = 0;
                 mcount++;
                 montex = montey = 0;
                 for (mj = 0; mj < MONTEN / 2; mj++) {
                    montex = (montex * 256.0) + monte[mj];
                    montey = (montey * 256.0) + monte[(MONTEN / 2) + mj];
                 }
                 if ((montex * montex + montey *  montey) <= incirc) {
                    inmont++;
                 }
              }
          }

          /* Update calculation of serial correlation coefficient */

          sccun = c;
          if (sccfirst) {
             sccfirst = FALSE;
             scclast = 0;
             sccu0 = sccun;
          } else {
             scct1 = scct1 + scclast * sccun;
          }
          scct2 = scct2 + sccun;
          scct3 = scct3 + (sccun * sccun);
          scclast = sccun;

          /* Update calculation for runs test */

          if(runsstart) {
             runsstart = FALSE;
          } else {
             if (c > median && runslast < median) {
                  runs++;
             }
             else if (c < median && runslast > median) {
                  runs++;
             }
          }
          runslast = c;

          if(c > median) {
             posval++;
          } else {
             negval++;
          }

          /* Update calculation of local means */

            if (totblock < M) {
               lmccount[c]++;
               totblock++;
            } else {
               double mean = 0.0;
               int i;

               for(i=0; i < 256; i++) {
                  mean += ((double) i) * lmccount[i];
               }
               mean /= M;

               if (local_means_count >= local_means_capacity) {
                  local_means_capacity *= 2;
                  double *new_ptr = realloc(local_means, local_means_capacity * sizeof(double));
               if (new_ptr == NULL) {
                  fprintf(stderr, "Error: failed to reallocate memory for local_means.\n");
                  free(local_means);
                  exit(EXIT_FAILURE);
               }
               local_means = new_ptr;
               }
               local_means[local_means_count++] = mean;

               for (i = 0; i < 256; i++) {
                  lmccount[i] = 0;
               }
               lmccount[c]++;
               totblock = 1;
               N++;
            }

          oc <<= 1;
       } while (binary && (++bean < 8));
    }
}

/*  RT_END  --  Complete calculation and return results.  */

void rt_end(double *r_ent, double *r_chisq, double *r_mean,
            double *r_montepicalc, double *r_scc, int *r_runs,
            double *r_runsz, int *r_N, double *r_lm_chisq)
{
    int i;

    /* Complete calculation of serial correlation coefficient */

    scct1 = scct1 + scclast * sccu0;
    scct2 = scct2 * scct2;
    scc = totalc * scct3 - scct2;
    if (scc == 0.0) {
       scc = -100000;
    } else {
       scc = (totalc * scct1 - scct2) / scc;
    }

    /* Scan bins and calculate probability for each bin and
       Chi-Square distribution.  The probability will be reused
       in the entropy calculation below.  While we're at it,
       we sum of all the data which will be used to compute the
       mean. */

    cexp = totalc / (binary ? 2.0 : 256.0);  /* Expected count per bin */

    /* rule-of-thumb check: χ² requires Eᵢ ≥ 5 for validity */
    if (cexp < 5.0) {
      fprintf(stderr,"Warning: expected count per bin (%.2f) < 5; "
         "χ² test may be unreliable.\n", cexp);
    }

    for (i = 0; i < (binary ? 2 : 256); i++) {
       double a = ccount[i] - cexp;;

       prob[i] = ((double) ccount[i]) / totalc;
       chisq += (a * a) / cexp;
       datasum += ((double) i) * ccount[i];
    }

    /* Calculate entropy */

    for (i = 0; i < (binary ? 2 : 256); i++) {
       if (prob[i] > 0.0) {
          ent += prob[i] * rt_log2(1 / prob[i]);
       }
    }

    /* Calculate Monte Carlo value for PI from percentage of hits
       within the circle */

    montepi = 4.0 * (((double) inmont) / mcount);

    /* Calculate expected number of runs and standard deviation of the number of runs for runs test*/

    expruns = (2 * posval * negval) / (posval + negval) + 1;

    stddevruns = (2 * posval * negval * (2 * posval * negval - posval - negval)) /
                 ((posval + negval) * (posval + negval) * (posval + negval - 1));

    runsz = (runs - expruns) / sqrt(stddevruns);

    /* Calculate Chi-Square distribution for local means test*/

    for (i = 0; i < local_means_count; i++) {
         double zstat = sqrt(M) * (local_means[i] - 127.5) / (sqrt((pow(256, 2) - 1) / 12.0));
         lm_chisq += (zstat * zstat);
      }

    /* Return results through arguments */

    *r_ent = ent;
    *r_chisq = chisq;
    *r_mean = datasum / totalc;
    *r_montepicalc = montepi;
    *r_scc = scc;
    *r_runs = runs;
    *r_runsz = runsz;
    *r_N = N;
    *r_lm_chisq = lm_chisq;
}

void rt_free() {
    if (local_means != NULL) {
        free(local_means);
        local_means = NULL;
    }
}
