/* dpend_out.c
 *
 * Used by dpend_mpfr.c to write out data.
 *
 */

#include "dpend.h"
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <mpfr.h>


void output_polar(mpfr_t* t, y_t* y) 
{
  mpfr_fprintf(polar_output, "%0.32RNf,%0.32RNf,%0.32RNF,%0.32RNF,%0.32RNF\n", 
              *t, y->th1, y->w1, y->th2, y->w2);
}

void output_lyapunov(mpfr_t* t, mpfr_t* exp) 
{
  mpfr_fprintf(lexp_output, "%0.32RNf,%0.32RNf\n", 
              *t, *exp);
}

void output_cartesian(mpfr_t* t, y_t* y) 
{

  /* Convert to cartesian coordinates. */
  mpfr_t x1, y1, x2, y2, temp;
  mpfr_inits2(nbits, x1, y1, x2, y2, temp, NULL);

  // Calculate x1
  mpfr_set_d(x1, L1, MPFR_RNDN);
  mpfr_sin(temp, y->th1, MPFR_RNDN);
  mpfr_mul(x1, x1, temp, MPFR_RNDN);

  // Calculate y1
  mpfr_set_d(y1, -L1, MPFR_RNDN);
  mpfr_cos(temp, y->th1, MPFR_RNDN);
  mpfr_mul(y1, y1, temp, MPFR_RNDN);

  // Calculate x2
  mpfr_set_d(x2, L2, MPFR_RNDN);
  mpfr_sin(temp, y->th2, MPFR_RNDN);
  mpfr_mul(x2, x2, temp, MPFR_RNDN);
  mpfr_add(x2, x2, x1, MPFR_RNDN);

  // Calculate y2
  mpfr_set_d(y2, -L2, MPFR_RNDN);
  mpfr_cos(temp, y->th2, MPFR_RNDN);
  mpfr_mul(y2, y2, temp, MPFR_RNDN);
  mpfr_add(y2, y2, y1, MPFR_RNDN);

  mpfr_fprintf(cartesian_output, "%0.32RNf,%0.32RNf,%0.32RNF,%0.32RNF,%0.32RNF\n", 
               *t, x1, y1, x2, y2); 
}

void output_energy(mpfr_t* t, y_t* y){
  
  mpfr_t p_energy, k_energy, t_energy, temp, temp2, cos_th1, cos_th2, del, cos_del;
  mpfr_inits2(nbits, p_energy, k_energy, t_energy, temp, temp2, cos_th1, cos_th2, del, cos_del, NULL);
  double half = 0.5;
  double neg = -1.0;
  double mass_sum = M1+M2;

  // PE
  mpfr_cos(cos_th1, y->th1, MPFR_RNDN);     
  mpfr_cos(cos_th2, y->th2, MPFR_RNDN);

  mpfr_mul_d(p_energy, cos_th1, L1, MPFR_RNDN); 
  mpfr_mul_d(p_energy, p_energy, mass_sum, MPFR_RNDN);
  mpfr_mul_d(temp, cos_th2, L2, MPFR_RNDN);
  mpfr_mul_d(temp, temp, M2, MPFR_RNDN);
  mpfr_add(p_energy, p_energy, temp, MPFR_RNDN);
  mpfr_mul_d(p_energy, p_energy, -G, MPFR_RNDN);

  //KE
  mpfr_mul(temp, y->w1, y->w1, MPFR_RNDN);
  mpfr_mul_d(temp, temp, L1, MPFR_RNDN);
  mpfr_mul_d(temp, temp, L1, MPFR_RNDN);
  mpfr_mul_d(temp, temp, M1, MPFR_RNDN);

  mpfr_mul(temp2, y->w2, y->w2, MPFR_RNDN);
  mpfr_mul_d(temp2, temp2, L2, MPFR_RNDN);
  mpfr_mul_d(temp2, temp2, L2, MPFR_RNDN);
  mpfr_mul_d(temp2, temp2, M2, MPFR_RNDN);

  mpfr_add(temp, temp, temp2, MPFR_RNDN);

  mpfr_mul(temp2, y->w1, y->w1, MPFR_RNDN);
  mpfr_mul_d(temp2, temp2, L1, MPFR_RNDN);
  mpfr_mul_d(temp2, temp2, L1, MPFR_RNDN);
  mpfr_mul_d(temp2, temp2, M2, MPFR_RNDN);

  mpfr_add(temp, temp, temp2, MPFR_RNDN);
  mpfr_mul_d(temp, temp, half, MPFR_RNDN);
  
  mpfr_mul(temp2, y->w1, y->w2, MPFR_RNDN);
  mpfr_mul_d(temp2, temp2, L1, MPFR_RNDN);
  mpfr_mul_d(temp2, temp2, L2, MPFR_RNDN);
  mpfr_mul_d(temp2, temp2, M2, MPFR_RNDN);
  mpfr_sub(del, y->th1, y->th2, MPFR_RNDN);
  mpfr_cos(cos_del, del, MPFR_RNDN); 
  mpfr_mul(temp2, temp2, cos_del, MPFR_RNDN);

  mpfr_add(k_energy, temp, temp2, MPFR_RNDN);
  mpfr_mul_d(temp2, p_energy, neg, MPFR_RNDN);
  mpfr_mul(k_energy, temp, temp2, MPFR_RNDN);

  //total energy
  mpfr_add(t_energy, p_energy, k_energy, MPFR_RNDN);

  mpfr_fprintf(energy_output, "%0.32RNf,%0.32RNf,%0.32RNF,%0.32RNF\n", 
               *t, p_energy, k_energy, t_energy); 

  
}


