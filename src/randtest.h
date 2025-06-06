
/*  Random test function prototypes  */

extern void rt_init(int binmode);
extern void rt_add(void *buf, int bufl);
extern void rt_end(double *r_ent, double *r_chisq, double *r_mean,
                   double *r_montepicalc, double *r_scc, int *r_runs,
                   double *r_runsz, int *r_N, double *r_lm_chisq);
extern void rt_free(void);