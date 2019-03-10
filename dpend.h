#include <stdio.h>
#include <gmp.h>
#include<mpfr.h>


#define N 4    /* number of equations to be solved */
#define G 9.8  /* gravity (m/s^2) */
#define L1 1.0 /* length of pendulum 1 (m) */
#define L2 1.0 /* length of pendulum 2 (m) */
#define M1 1.0 /* mass of pendulum 1 (kg) */
#define M2 1.0 /* mass of pendulum 2 (kg) */

int nbits;  /* number of bits to use for mantissa */
mpfr_t h;   /* step size */

typedef struct {
  mpfr_t th1;     /* angle of pend 1 */
  mpfr_t w1;      /* angular velocity of pend 1 */
  mpfr_t th2;     /* angle of pend 2 */
  mpfr_t w2;      /* angular velocity of pend 2 */
} y_t;

FILE *polar_output;
FILE *cartesian_output;
FILE *energy_output;
FILE *lexp_output;

FILE *final_lexp_output;
FILE *all_ics;
