/* 
 * dpend_mpfr.c - Arbitrary precision solution to double pendulum ODEs
 * 		  using fourth order Runge-Kutta. 
 *
 * Parameters are passed in at the command line:
 * 
 * $./dpend_mpfr.c TH10 W10 TH20 W20 NSTEP
 *
 * TH10, TH20 - initial angles of the pendulums (degrees)
 * W10, W20 - initial angular velocities of the pendulums (degrees per second)
 * NSTEP - number of times steps (default step size is h=0.0001)
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>
#include "dpend.h"
#include "dpend_out.h"
#include "dpend_math.h"
#include <dirent.h>
#include <sys/stat.h>
#include <errno.h>
#include <math.h>

int main(int argc, char *argv[])
{
  int NSTEP = atoi(argv[5]);
  int mantissa_sz[5] = {113,11,24,53,64};
 // int eps[5] = {0, 2^(-11), 2^(-24), 2^(-53), 2^(-64)};
  h = pow( atoi(argv[6]), atoi(argv[7]) );  
  create_output_directory();
  output_initial_conditions(argv);

  

  y_t * y_actual;
  y_actual = malloc(NSTEP*sizeof(y_t));


  //mpfr_t radian_conv;

  /* Preform Runga Kutta method to solve the dpend system for each mantissa size. */
  for (int j = 0; j < 5; j++) {
    mpfr_t t_curr, t_next;
    mpfr_t mag, dot, r_error;
    y_t yin, yout;

    nbits = mantissa_sz[j];

    /* Create constant for converting angles to radians. */

    //mpfr_init2(h, 113);
    // mpfr_set_d(h, 0.001, MPFR_RNDN);
    mpfr_t radian_conv;
    mpfr_init2(radian_conv, nbits);
    //mpfr_const_pi(radian_conv, MPFR_RNDN);
    mpfr_set_d(radian_conv, 3.14159, MPFR_RNDN);  
    mpfr_div_ui(radian_conv, radian_conv, 180, MPFR_RNDN);

    mpfr_inits2(113, mag, dot, r_error, NULL);
    mpfr_set_d(r_error, 0.0, MPFR_RNDN);

    mpfr_inits2(nbits, yin.th1, yin.w1, yin.th2, yin.w2, yout.th1, yout.w1, yout.th2, yout.w2, NULL);

   // printf("\nA: underflow:  %d   overflow: %d   div0: %d   nanflag: %d   inexflag: %d   erangeflag: %d\n",
      //  mpfr_underflow_p(), mpfr_overflow_p(), mpfr_divby0_p(), mpfr_nanflag_p(), mpfr_inexflag_p(), mpfr_erangeflag_p());

    /* Create output files to hold results for calculations using nbits. */
    create_out_files(mantissa_sz[j]);

    /* Initilize and set time values. */
    mpfr_inits2(113, t_curr, t_next, NULL);
    mpfr_set_d(t_curr, 0.0, MPFR_RNDN);
    mpfr_set_d(t_next, 0.0, MPFR_RNDN);

    /* Set initial values converting angles to radians. */
    mpfr_mul_d(yin.th1, radian_conv, atof(argv[1]), MPFR_RNDN);  // th1[0] = TH1*PI/180.0
    mpfr_mul_d(yin.w1, radian_conv, atof(argv[2]) , MPFR_RNDN);    // w1[0] = W1*PI/180.0
    mpfr_mul_d(yin.th2, radian_conv, atof(argv[3]), MPFR_RNDN);  // th2[0] = TH2*PI/180.0
    mpfr_mul_d(yin.w2, radian_conv, atof(argv[4]) , MPFR_RNDN);    // w2[0] = W2*PI/180.0

   // printf("B: underflow:  %d   overflow: %d   div0: %d   nanflag: %d   inexflag: %d   erangeflag: %d\n",
     //   mpfr_underflow_p(), mpfr_overflow_p(), mpfr_divby0_p(), mpfr_nanflag_p(), mpfr_inexflag_p(), mpfr_erangeflag_p());

    magnitude(&yin, &mag);
    if (nbits != 113) {
      dot_product(&yin, &y_actual[0], &dot);
     // relative_error(&y_actual[0], &yout, &r_error);
      output_mag(&t_curr, &mag, &dot, &r_error);
    } else {
      y_actual[0] = yin;
      output_mag(&t_curr, &mag, NULL, NULL);
    }

    /* Output initial values. */
    output_polar(&t_curr, &yin);
    output_cartesian(&t_curr, &yin);
    //output_energy(&t_curr,  &yin);


    printf("Computing %d-bit result... ", nbits); fflush(stdout);
    /* Perform the integration. */
    for (int i = 1; i < NSTEP; i++) {
      mpfr_add_d(t_next, t_curr, h, MPFR_RNDN);		// update time
      runge_kutta(t_curr, &yin, &yout);      // preform runge kutta 
      
      /* Print output to files. */
      output_polar(&t_next, &yout);
     // output_cartesian(&t_next, &yout);

      //output_energy(&t_next, &yout);

      /* Set yin to yout. */
      mpfr_set(yin.th1, yout.th1, MPFR_RNDN);
      mpfr_set(yin.w1, yout.w1, MPFR_RNDN);
      mpfr_set(yin.th2, yout.th2, MPFR_RNDN);
      mpfr_set(yin.w2, yout.w2, MPFR_RNDN);
      mpfr_set(t_curr, t_next, MPFR_RNDN);

      if (nbits == 113) { y_actual[i] = yin; }

      magnitude(&yout, &mag);

      if (nbits != 113) { 
        dot_product(&yout, &y_actual[i], &dot);
        relative_error(&y_actual[i], &yout, &r_error);
        output_mag(&t_curr, &mag, &dot, &r_error);
      } else {
        output_mag(&t_curr, &mag, NULL, NULL);
      }

    }
//mpfr_clear_flags();
    printf("Done.\n");


    mpfr_clears(radian_conv, NULL);
    mpfr_clears(mag, dot, t_curr, t_next, yin.th1, yin.w1, yin.th2, yin.w2,
  		yout.th1, yout.w1, yout.th2, yout.w2, NULL);
    mpfr_free_cache();


    /* Close files. */
    fclose(polar_output);
    fclose(cartesian_output);
    fclose(mag_output);
    //fclose(energy_output);
  }


  
  free(y_actual);

  return 0;
}






















