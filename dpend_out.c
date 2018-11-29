/* dpend_out.c
 *
 * Used by dpend_mpfr.c to write out data.
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <mpfr.h>


void output_polar(FILE *file, mpfr_t* t, mpfr_t* th1, mpfr_t* w1, mpfr_t* th2, mpfr_t* w2) {
  mpfr_fprintf(file, "%0.32RNf %0.32RNf %0.32RNF %0.32RNF %0.32RNF\n", 
              *t, *th1, *w1, *th2, *w2);
}

void output_cartesian(FILE *file, int nbits, mpfr_t* t, mpfr_t* th1, mpfr_t* w1, mpfr_t* th2, mpfr_t* w2, double L1, double L2) {

  /* Convert to cartesian coordinates. */
  mpfr_t x1, y1, x2, y2, temp;
  mpfr_inits2(nbits, x1, y1, x2, y2, temp, NULL);

  // Calculate x1
  mpfr_set_d(x1, L1, MPFR_RNDN);
  mpfr_sin(temp, th1, MPFR_RNDN);
  mpfr_mul(x1, x1, temp, MPFR_RNDN);

  // Calculate y1
  mpfr_set_d(y1, -L1, MPFR_RNDN);
  mpfr_cos(temp, th1, MPFR_RNDN);
  mpfr_mul(y1, y1, temp, MPFR_RNDN);

  // Calculate x2
  mpfr_set_d(x2, L2, MPFR_RNDN);
  mpfr_sin(temp, th2, MPFR_RNDN);
  mpfr_mul(x2, x2, temp, MPFR_RNDN);
  mpfr_add(x2, x2, x1, MPFR_RNDN);

  // Calculate y2
  mpfr_set_d(y2, -L2, MPFR_RNDN);
  mpfr_cos(temp, th2, MPFR_RNDN);
  mpfr_mul(y2, y2, temp, MPFR_RNDN);
  mpfr_add(y2, y2, y1, MPFR_RNDN);

  mpfr_fprintf(file, "%0.32RNf %0.32RNf %0.32RNF %0.32RNF %0.32RNF\n", 
               t, x1, y1, x2, y2); 
}
