#include <stdio.h>
#include <math.h>

#ifndef PATPC_NOGSL
#include <gsl/gsl_rng.h>
#endif

int main(void){

#ifdef PATPC_NOGSL
 // do nothing if we are not allowed to use GSL
 fprintf(stderr, "Please install GSL and run 'make' to use this simulator!\n");
 return 1;
#else

 double time_sec;
 double period_sec=400;
 double phase;

 const gsl_rng_type * T;
 gsl_rng * r;

 int i, n = 10;

 gsl_rng_env_setup();

 T = gsl_rng_default;
 r = gsl_rng_alloc (T);

 for ( time_sec=0; time_sec<20000; time_sec+=period_sec ) {
  // background uniform in phase
  for (i = 0; i < 2*n; i++){
   phase = gsl_rng_uniform (r);
   printf ("%.5f\n", time_sec+phase*period_sec);
  }
  // excess in phase 0.0-0.25
  for (i = 0; i < n; i++){
   phase = gsl_rng_uniform (r) / 4.0;
   printf ("%.5f\n", time_sec+phase*period_sec);
  }
 }

 gsl_rng_free (r);

#endif

 return 0;
}
