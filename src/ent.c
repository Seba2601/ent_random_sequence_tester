/*
        ENT  --  Entropy calculation and analysis of putative
                 random sequences.

        Designed and implemented by John "Random" Walker in May 1985.

        Multiple analyses of random sequences added in December 1985.

        Bit stream analysis added in September 1997.

        Terse mode output, getopt() command line processing,
        optional stdin input, and HTML documentation added in
        October 1998.

        Documentation for the -t (terse output) option added
        in July 2006.

        Replaced table look-up for chi square to probability
        conversion with algorithmic computation in January 2008.

        For additional information and the latest version,
        see https://www.fourmilab.ch/random/

*/

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#else
#include <unistd.h>
#endif

#include "iso8859.h"
#include "randtest.h"

#define FALSE 0
#define TRUE  1

#ifdef M_PI
#define PI       M_PI
#else
#define PI       3.14159265358979323846
#endif

extern double poz(const double z);

extern double pochisq(const double ax, const int df);

/*  HELP  --  Print information on how to call  */

static void help(void)
{
        printf("ent --  Test randomness of file.  Call with");
        printf("\n            ent [options] [input-file]");
        printf("\n");
        printf("\n        Options:   -b   Treat input as a stream of bits");
        printf("\n                   -c   Print occurrence counts");
        printf("\n                   -f   Fold upper to lower case letters");
        printf("\n                   -t   Terse output in CSV format");
        printf("\n                   -p   Include p-values");
        printf("\n                   -u   Print this message\n");
        printf("\nVersion " VERSION);
        printf("\nBy John Walker");
        printf("\n   https://www.fourmilab.ch/\n");
}

/*  GETOPT  --  Dumb version of getopt for brain-dead Windows.  */

#ifdef _WIN32
static int optind = 1;

static int getopt(int argc, char *argv[], char *opts)
{
    static char *opp = NULL;
    int o;

    while (opp == NULL) {
        if ((optind >= argc) || (*argv[optind] != '-')) {
           return -1;
        }
        opp = argv[optind] + 1;
        optind++;
        if (*opp == 0) {
            opp = NULL;
        }
    }
    o = *opp++;
    if (*opp == 0) {
        opp = NULL;
    }
    return strchr(opts, o) == NULL ? '?' : o;
}
#endif

/*  Main program  */

int main(int argc, char *argv[])
{
        int i, oc, opt;
        long ccount[256];             /* Bins to count occurrences of values */
        long totalc = 0;              /* Total character count */
        char *samp;

        double ent;                   /*Variables for Entropy test*/
        double chisq, chip;           /*Variables for Chi-square test*/
        double mean, meanz, meanp;    /*Variables for Mean test*/
        double montepi, mcount,       /* Variables for Monte Carlo test*/
               montepiz, montepip;
        double scc, scct, sccp;       /* Variables for Serial Correlation test*/
        double runsz, runsp;          /* Variables for Runs test*/
        int runs;
        
        int N;                        /* Variables for Local Means test*/     
        double lm_chisq,locmeanp;            

        FILE *fp = stdin;
        int counts = FALSE,           /* Print character counts */
            fold = FALSE,             /* Fold upper to lower */
            binary = FALSE,           /* Treat input as a bitstream */
            terse = FALSE,            /* Terse (CSV format) output */ 
            csp = FALSE;              /* Include p-values */

        while ((opt = getopt(argc, argv, "bcftpuv?BCFTPUV")) != -1) {
            switch (toISOlower(opt)) {
                 case 'b':
                    binary = TRUE;
                    break;

                 case 'c':
                    counts = TRUE;
                    break;

                 case 'f':
                    fold = TRUE;
                    break;

                 case 't':
                    terse = TRUE;
                    break;

                 case 'p':
                    csp = TRUE;
                    break;

                 case '?':
                 case 'u':
                 case 'v':
                    help();
                    return 0;
            }
        }
        if (optind < argc) {
           if (optind != (argc - 1)) {
              printf("Duplicate file name.\n");
              help();
              return 2;
           }
           if ((fp = fopen(argv[optind], "rb")) == NULL) {
              printf("Cannot open file %s\n", argv[optind]);
              return 2;
           }
        }

#ifdef _WIN32

            /** Warning!  On systems which distinguish text mode and
                binary I/O (MS-DOS, Macintosh, etc.) the modes in the open
                statement for "fp" should have forced the input file into
                binary mode.  But what if we're reading from standard
                input?  Well, then we need to do a system-specific tweak
                to make sure it's in binary mode.  While we're at it,
                let's set the mode to binary regardless of however fopen
                set it.

                The following code, conditional on _WIN32, sets binary
                mode using the method prescribed by Microsoft Visual C 7.0
                ("Monkey C"); this may require modification if you're
                using a different compiler or release of Monkey C.      If
                you're porting this code to a different system which
                distinguishes text and binary files, you'll need to add
                the equivalent call for that system. */

            _setmode(_fileno(fp), _O_BINARY);
#endif

        samp = binary ? "bit" : "byte";
        memset(ccount, 0, sizeof ccount);

        /* Initialise for calculations */

        rt_init(binary);

        /* Scan input file and count character occurrences */

        while ((oc = fgetc(fp)) != EOF) {
           unsigned char ocb;

           if (fold && isISOalpha(oc) && isISOupper(oc)) {
              oc = toISOlower(oc);
           }
           ocb = (unsigned char) oc;
           totalc += binary ? 8 : 1;
           if (binary) {
            int b;
            unsigned char ob = ocb;

            for (b = 0; b < 8; b++) {
                ccount[ob & 1]++;
                ob >>= 1;
            }
           } else {
               ccount[ocb]++;         /* Update counter for this bin */
           }
           rt_add(&ocb, 1);
        }
        fclose(fp);

        /* Complete calculation */

        rt_end(&ent, &chisq, &mean, &montepi, &scc, &runs, &runsz, &N, &lm_chisq);

        /* Calculate p-value for Chi-square test*/

        chip = pochisq(chisq, (binary ? 1 : 255));

        meanz = (sqrt(totalc) * (mean - 127.5)) / (sqrt((pow(256, 2) - 1) / 12.0));

        meanp = 2 * (1 - poz(fabs(meanz)));

        mcount = floor(totalc / 6.0);

        montepiz = sqrt(mcount) * (montepi - PI) / (4 * sqrt(PI/4 * (1 - PI/4)));

        montepip = 2 * (1 - poz(fabs(montepiz)));

        scct = (scc * sqrt(totalc - 2)) / sqrt(1 - (scc * scc));

        sccp = 2 * (1 - poz(fabs(scct) * (1 - 1 / (4 * totalc)) * pow((1 + (pow(scct, 2)) / (2 * totalc)), -0.5)));

        runsp = 2 * (1 - poz(fabs(runsz)));

        locmeanp = pochisq(lm_chisq, N);

        /* Print terse output if requested */

        if (terse) {
            if (csp) {
                  printf("0,File-%ss,Entropy,Chi-square,Chi-square-p-val,Mean,Mean-p-val,Monte-Carlo-Pi,Monte-Carlo-Pi-p-val,Serial-Correlation,Serial-Correlation-p-val,Runs,Runs-p-val,Local-Means-X^2,Local-Means-p-val\n", binary ? "bit" : "byte");
                  printf("1,%ld,%.10f,%f,%f,%f,%f,%f,%f,%f,%f,%d,%f,%f,%f\n", totalc, ent, chisq, chip, mean, meanp, montepi, montepip, scc, sccp, runs, runsp, lm_chisq, locmeanp);
            } else {
                  printf("0,File-%ss,Entropy,Chi-square,Mean,Monte-Carlo-Pi,Serial-Correlation,Runs,Local-Means-X^2\n", binary ? "bit" : "byte");
                  printf("1,%ld,%.10f,%f,%f,%f,%f,%d,%f\n", totalc, ent, chisq, mean, montepi, scc, runs, lm_chisq);
            }  
        }

        /* Print bin counts if requested */

        if (counts) {
           if (terse) {
              printf("2,Value,Occurrences,Fraction\n");
           } else {
              printf("Value Char Occurrences Fraction\n");
           }
           for (i = 0; i < (binary ? 2 : 256); i++) {
              if (terse) {
                 printf("3,%d,%ld,%f\n", i,
                    ccount[i], ((double) ccount[i] / totalc));
              } else {
                 if (ccount[i] > 0) {
                    printf("%3d   %c   %10ld   %f\n", i,
                       /* The following expression shows ISO 8859-1
                          Latin1 characters and blanks out other codes.
                          The test for ISO space replaces the ISO
                          non-blanking space (0xA0) with a regular
                          ASCII space, guaranteeing it's rendered
                          properly even when the font doesn't contain
                          that character, which is the case with many
                          X fonts. */
                       (!isISOprint(i) || isISOspace(i)) ? ' ' : i,
                       ccount[i], ((double) ccount[i] / totalc));
                 }
              }
           }
           if (!terse) {
              printf("\nTotal:    %10ld   %f\n\n", totalc, 1.0);
           }
        }

        /* Print calculated results */

        if (!terse && !csp) {
           printf("StringENT | Results report\n");
           printf("--------------------------------------------------\n");
           printf("Entropy = %.10f bits per %s.", ent, samp);
           printf("\nOptimum compression would reduce the size\n");
           printf("of this %ld %s file by %d percent.\n", totalc, samp,
                    (short) ((100 * ((binary ? 1 : 8) - ent) /
                              (binary ? 1.0 : 8.0))));
           printf("--------------------------------------------------\n");
           printf(
              "Chi square distribution for %ld samples is %1.2f, and randomly\n",
              totalc, chisq);
           if (chip < 0.0001) {
              printf("would exceed this value less than 0.01 percent of the times.\n");
           } else if (chip > 0.9999) {
              printf("would exceed this value more than than 99.99 percent of the times.\n");
           } else {
              printf("would exceed this value %1.2f percent of the times.\n",
                 chip * 100);
           }
           printf("--------------------------------------------------\n");
           printf(
              "Arithmetic mean value of data %ss is %1.4f (%.1f = random).\n",
              samp, mean, binary ? 0.5 : 127.5);
           printf("--------------------------------------------------\n");
           printf("Monte Carlo value for Pi is %1.9f (error %1.2f percent).\n",
              montepi, 100.0 * (fabs(PI - montepi) / PI));
           printf("--------------------------------------------------\n");
           printf("Serial correlation coefficient is ");
           if (scc >= -99999) {
              printf("%1.6f (totally uncorrelated = 0.0).\n", scc);
           } else {
              printf("undefined (all values equal!).\n");
           }
           printf("--------------------------------------------------\n");
           printf("The number of runs test is %d runs.\n", runs);
           printf("--------------------------------------------------\n");
           printf("The local means test's X^2 statistic is %f for %d blocks.\n", lm_chisq, N);
           printf("--------------------------------------------------\n");
        } else if (!terse) {
           printf("StringENT | Results report\n");
           printf("--------------------------------------------------\n");
           printf("Entropy = %.10f bits per %s.", ent, samp);
           printf("\nOptimum compression would reduce the size\n");
           printf("of this %ld %s file by %d percent.\n", totalc, samp,
                    (short) ((100 * ((binary ? 1 : 8) - ent) /
                              (binary ? 1.0 : 8.0))));
           printf("--------------------------------------------------\n");
           printf(
              "Chi square distribution for %ld samples is %1.2f.\n",
              totalc, chisq);
           printf("p-value   %f\n", chip);
           printf("--------------------------------------------------\n");
           printf(
              "Arithmetic mean value of data %ss is %1.4f (%.1f = random).\n",
              samp, mean, binary ? 0.5 : 127.5);
           printf("p-value   %f\n", meanp);   
           printf("--------------------------------------------------\n");
           printf("Monte Carlo value for Pi is %1.9f (error %1.2f percent).\n",
              montepi, 100.0 * (fabs(PI - montepi) / PI));
           printf("p-value   %f\n", montepip);
           printf("--------------------------------------------------\n");
           printf("Serial correlation coefficient is ");
           if (scc >= -99999) {
              printf("%1.6f (totally uncorrelated = 0.0).\n", scc);
           } else {
              printf("undefined (all values equal!).\n");
           }
           printf("p-value   %f\n", sccp);
           printf("--------------------------------------------------\n");
           printf("The number of runs test is %d runs.\n", runs);
           printf("p-value   %f\n", runsp);
           printf("--------------------------------------------------\n");
           printf("The local means test's X^2 statistic is %f for %d blocks.\n", lm_chisq, N);
           printf("p-value   %f\n", locmeanp);
           printf("--------------------------------------------------\n");  
        }

        return 0;
}
