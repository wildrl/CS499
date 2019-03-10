/* 
 * dpend_mpfr.c - Arbitrary precision solution to double pendulum ODEs
 * 		  using fourth order Runge-Kutta. 
 *
 * Parameters are passed in at the command line:
 * 
 * $./dpend_mpfr.c TMIN TMAX TH10 W10 TH20 W20 NSTEP BITS
 $./dpend_mpfr.c TMIN TMAX TH10 W10 TH20 W20 BITS
 *
 * TMIN, TMAX - start and end times (seconds)
 * TH10, TH20 - initial angles of the pendulums (degrees)
 * W10, W20 - initial angular velocities of the pendulums (degrees per second)
 * NSTEP - number of times steps
 * BITS - number of bits to use for the significand 
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

int main(int argc, char *argv[])
{
  unsigned int NSTEP;
  int mantissa_sz[5] = {11,24,53,64,113};

  /* Count number of files in mpfr_data in order to name the new file to be created. */
  int dir_count = 0;
  DIR *dirp;
  struct dirent *entry;
    
  dirp = opendir("./mpfr_data");
  while ((entry = readdir(dirp)) != NULL) {
    if (entry->d_type == DT_DIR) { /* If the entry is a regular file */
      dir_count++;
    }
  }
  closedir(dirp);

  /* Create ic_## directory to store output in. */
  char dir_name[30];
  snprintf(dir_name, 30, "./mpfr_data/ic_%d", dir_count);
  mkdir(dir_name, 0777);

  /* Record the initial conditions for this execution. */
  all_ics = fopen("./mpfr_data/all_ics.txt", "a");
  fprintf(all_ics, "ic_%d,%s,%s,%s,%s,%s,%s,%s\n", 
          dir_count, argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7]);
  fclose(all_ics);

  /* Open file to record final lyapunov exponents in for each run. */
  //final_lexp_output = fopen("./mpfr_data/final_lexps.csv", "a");

  // Preform Runga Kutta method to solve the dpend system for each mantissa size. */
  for (int j = 0; j < 5; j++) {

    nbits = mantissa_sz[j];

    /* Create output file for polar coordinates. */
    char polar_fn[50];
    snprintf(polar_fn, 50, "%s/polar%d.csv", dir_name, nbits);
    polar_output = fopen(polar_fn, "w");
    fprintf(polar_output, "time,th1,w1,th2,w2\n");

    /* Create output file for cartesian coordinates. */
    char cartesian_fn[50];
    snprintf(cartesian_fn, 50, "%s/cartesian%d.txt", dir_name, nbits);
    cartesian_output = fopen(cartesian_fn, "w");
    fprintf(cartesian_output, "time,x1,y1,x2,y2\n");

    /* Create output file for lypapunov exponent. *2/
    char lexp_fn[50];
    snprintf(lexp_fn, 50, "%s/lexp%d.csv", dir_name, nbits);
    lexp_output = fopen(lexp_fn, "w");
    fprintf(lexp_output, "time,lexp\n"); */

    /* Create output file for energy values. */
    //char energy_fn[50];
    //snprintf(energy_fn, 50, "%s/energy%d.csv", dir_name, nbits);
    //energy_output = fopen(energy_fn, "w");
    //fprintf(energy_output, "time,KE,PE,Total\n");

    mpfr_t h, TMIN, TMAX, t_curr, t_next, TH10, W10, TH20, W20;
    mpfr_inits2(nbits, TMIN, TMAX, TH10, W10, TH20, W20, NULL);
    mpfr_inits2(113, t_curr, t_next, h, NULL);

    y_t yin, yout;
    mpfr_inits2(nbits, yin.th1, yin.w1, yin.th2, yin.w2, yout.th1, yout.w1, yout.th2, yout.w2, NULL);

    /* obtain command line values */
    mpfr_set_flt(TMIN, atof(argv[1]), MPFR_RNDN);
    mpfr_set_flt(TMAX, atof(argv[2]), MPFR_RNDN);
    mpfr_set_flt(TH10, atof(argv[3]), MPFR_RNDN);
    mpfr_set_flt(W10, atof(argv[4]), MPFR_RNDN);
    mpfr_set_flt(TH20, atof(argv[5]), MPFR_RNDN);
    mpfr_set_flt(W20, atof(argv[6]), MPFR_RNDN);

    mpfr_set_d(h, 0.0001, MPFR_RNDN);
    NSTEP = (int) (atof(argv[2])-atof(argv[1]))/0.0001;

    /* Create constant for converting angles to radians. */
    mpfr_t radian_conv;
    mpfr_init2(radian_conv, nbits);
    mpfr_const_pi(radian_conv, MPFR_RNDN);
    mpfr_div_si(radian_conv, radian_conv, 180, MPFR_RNDN);

    /* Set initial values converting angles to radians. */
    mpfr_set(t_curr, TMIN, MPFR_RNDN);
    mpfr_set_d(t_next, 0.0, MPFR_RNDN);
    mpfr_mul(yin.th1, TH10, radian_conv, MPFR_RNDN);  // th1[0] = TH10*PI/180.0;
    mpfr_mul(yin.w1, W10, radian_conv, MPFR_RNDN);    // w1[0] = W10*PI/180.0;
    mpfr_mul(yin.th2, TH20, radian_conv, MPFR_RNDN);  // th2[0] = TH20*PI/180.0;
    mpfr_mul(yin.w2, W20, radian_conv, MPFR_RNDN);    // w2[0] = W20*PI/180.0; 

    /* Clean up. */
    mpfr_clears(TMIN, TMAX, TH10, W10, TH20, W20, radian_conv, NULL);

    /* Output initial values. */
    output_polar(&t_curr, &yin);
    output_cartesian(&t_curr, &yin);
    //output_energy(&t_curr,  &yin);

    /* Perform the integration. */
    for (int i = 0; i < 100000; i++) {

      mpfr_add(t_next, t_curr, h, MPFR_RNDN);		// update time
      runge_kutta(t_curr, &yin, &yout, h);      // preform runge kutta 

      /* Print output to files. */
      output_polar(&t_next, &yout);
      output_cartesian(&t_next, &yout);
      //output_energy(&t_next, &yout);

      /* Set yin to yout. */
      mpfr_set(yin.th1, yout.th1, MPFR_RNDN);
      mpfr_set(yin.w1, yout.w1, MPFR_RNDN);
      mpfr_set(yin.th2, yout.th2, MPFR_RNDN);
      mpfr_set(yin.w2, yout.w2, MPFR_RNDN);
      mpfr_set(t_curr, t_next, MPFR_RNDN);

    }
    printf("done with integration loop\n");

    mpfr_clears(h, t_curr, t_next, yin.th1, yin.w1, yin.th2, yin.w2,
  		yout.th1, yout.w1, yout.th2, yout.w2, NULL);
    mpfr_free_cache();

    /* Close files. */
    fclose(polar_output);
    fclose(cartesian_output);
    //fclose(energy_output);
  }

  return 0;
}






















