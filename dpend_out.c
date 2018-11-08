/* dpend_out.c
 *
 * Used by dpend_mpfr.c to write out data.
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <mpfr.h>


void output_polar(FILE *file, seconds_t t, y_t *yin) {
  mpfr_fprintf(file, "%0.32RNf %0.32RNf %0.32RNF %0.32RNF %0.32RNF\n", 
              t, yin->th1, yin->w1, yin->th2, yin->w2);
}

void output_cartesian(FILE *file, seconds_t t, y_t *yin, double L1, double L2) {

  /* Convert to cartesian coordinates. */
  mpfr_t x1, y1, x2, y2, temp;
  mpfr_inits2(nbits, x1, y1, x2, y2, temp, NULL);

  // Calculate x1
  mpfr_set_d(x1, L1, MPFR_RNDN);
  mpfr_sin(temp, yin->th1, MPFR_RNDN);
  mpfr_mul(x1, x1, temp, MPFR_RNDN);

  // Calculate y1
  mpfr_set_d(y1, -L1, MPFR_RNDN);
  mpfr_cos(temp, yin->th1, MPFR_RNDN);
  mpfr_mul(y1, y1, temp, MPFR_RNDN);

  // Calculate x2
  mpfr_set_d(x2, L2, MPFR_RNDN);
  mpfr_sin(temp, yin->th2, MPFR_RNDN);
  mpfr_mul(x2, x2, temp, MPFR_RNDN);
  mpfr_add(x2, x2, x1, MPFR_RNDN);

  // Calculate y2
  mpfr_set_d(y2, -L2, MPFR_RNDN);
  mpfr_cos(temp, yin->th2, MPFR_RNDN);
  mpfr_mul(y2, y2, temp, MPFR_RNDN);
  mpfr_add(y2, y2, y1, MPFR_RNDN);

  mpfr_fprintf(cartesian_output, "%0.32RNf %0.32RNf %0.32RNF %0.32RNF %0.32RNF\n", 
               t_curr, x1, y1, x2, y2); 
}
