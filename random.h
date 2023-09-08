/* Prototypes for all user accessible RANLIB routines */

#ifndef __NETLIB_RANDOM__
#define __NETLIB_RANDOM__

#ifdef __cplusplus
extern "C" {
#endif

void advnst (long k);
float genbet (float aa, float bb);
float genchi (float df);
float genexp (float av);
float genf (float dfn, float dfd);
float gengam (float a, float r);
void genmn (float* parm, float* x, float* work);
void genmul (long n, float* p, long ncat, long* ix);
float gennch (float df, float xnonc);
float gennf (float dfn, float dfd, float xnonc);
float gennor (float av, float sd);
void genprm (long* iarray, int larray);
float genunf (float low, float high);
void getsd (long* iseed1, long* iseed2);
void gscgn (long getset, long* g);
long ignbin (long n, float pp);
long ignnbn (long n, float p);
long ignlgi (void);
long ignpoi (float mu);
long ignuin (long low, long high);
void initgn (long isdtyp);
long mltmod (long a, long s, long m);
void phrtsd (char* phrase, long* seed1, long* seed2);
float ranf (void);
void setall (long iseed1, long iseed2);
void setant (long qvalue);
void setgmn (float* meanv, float* covm, long p, float* parm);
void setsd (long iseed1, long iseed2);
float sexpo (void);
float sgamma (float a);
float snorm (void);

#ifdef __cplusplus
}
#endif

#endif
