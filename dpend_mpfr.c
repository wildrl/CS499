/* 
 * dpend_mpfr.c - Arbitrary precision solution to double pendulum ODEs
 * 		  using fourth order Runge-Kutta. 
 *
 * Parameters are passed in at the command line:
 * 
 * $./dpend_mpfr.c TMIN TMAX TH10 W10 TH20 W20 NSTEP BITS
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
#include "dpend_out.h"
#include <dirent.h>
#include <sys/stat.h>
#include <errno.h>

/* hardwired parameters */

#define N 4 /* number of equations to be solved */
#define G 9.8 /* acc'n due to gravity, in m/s^2 */
#define L1 1.0 /* length of pendulum 1 in m */
#define L2 1.0 /* length of pendulum 2 in m */
#define M1 1.0 /* mass of pendulum 1 in kg */
#define M2 1.0 /* mass of pendulum 2 in kg */

int nbits;

typedef struct {
  mpfr_t th1;
  mpfr_t w1;
  mpfr_t th2;
  mpfr_t w2;
} y_t;

FILE *polar_output;
FILE *cartesian_output;
FILE *energy_output;
FILE *ly_exp_output;

FILE *final_exp_output;
FILE *ic_record;

void runge_kutta(mpfr_t t, y_t *yin, y_t *yout, mpfr_t h);
void derivs(y_t *yin, y_t *dydx);
void reset_yin_adj(mpfr_t *d0, mpfr_t *di, y_t *y0, y_t *y1_out, y_t *y1_in);
void lyapunov(mpfr_t *sum, mpfr_t *d0, mpfr_t *d1);
void calc_di(mpfr_t *di, y_t *y0, y_t *y1);

int main(int argc, char *argv[])
{
  unsigned int NSTEP;
  int mantissas[5] = {11,24,53,64,113};

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

  /* Create directory to store output in. */
  char dir_name[30];
  snprintf(dir_name, 30, "./mpfr_data/ic_%d", dir_count);
  mkdir(dir_name, 0777);


  final_exp_output = fopen("./mpfr_data/final_exp.csv", "a");

  ic_record = fopen("./mpfr_data/ic_record.txt", "a");
  fprintf(ic_record, "ic_%d    (%s, %s, %s, %s, %s, %s, %s)\n", dir_count, argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7]);


  for (int j = 0; j < 5; j++) {

    nbits = mantissas[j];

    /* Create output files. */
    char polar_fn[50];
    char ly_exp_fn[50];
    char cartesian_fn[50];
    //char energy_fn[30];
    snprintf(polar_fn, 50, "%s/polar%d.csv", dir_name, nbits);
    snprintf(ly_exp_fn, 50, "%s/ly_exp%d.csv", dir_name, nbits);
    snprintf(cartesian_fn, 50, "%s/cartesian%d.txt", dir_name, nbits);
    //snprintf(energy_fn, 30, "./mpfr_data/energy%s.csv", nbits);
    polar_output = fopen(polar_fn, "w");  
    ly_exp_output = fopen(ly_exp_fn, "w");
    //fseek(final_exp_output, 0, SEEK_END);
    cartesian_output = fopen(cartesian_fn, "w");
    //energy_output = fopen(energy_fn, "w");

    mpfr_t h, TMIN, TMAX, t_curr, t_next, TH10, W10, TH20, W20;
    mpfr_inits2(nbits, h, TMIN, TMAX, t_curr, t_next, TH10, W10, TH20, W20, NULL);

    y_t yin, yout;
    mpfr_inits2(nbits, yin.th1, yin.w1, yin.th2, yin.w2, yout.th1, yout.w1, yout.th2, yout.w2, NULL);

    /* obtain command line values */
    mpfr_set_flt(TMIN, atof(argv[1]), MPFR_RNDN);
    mpfr_set_flt(TMAX, atof(argv[2]), MPFR_RNDN);
    mpfr_set_flt(TH10, atof(argv[3]), MPFR_RNDN);
    mpfr_set_flt(W10, atof(argv[4]), MPFR_RNDN);
    mpfr_set_flt(TH20, atof(argv[5]), MPFR_RNDN);
    mpfr_set_flt(W20, atof(argv[6]), MPFR_RNDN);
    NSTEP = atoi(argv[7]);

    /* Calculate stepsize for integration */
    mpfr_sub(h, TMAX, TMIN, MPFR_RNDN);
    NSTEP--;
    mpfr_div_ui(h, h, NSTEP, MPFR_RNDN);  //h = (TMAX - TMIN)/(NSTEP - 1.0)
    NSTEP++;

    /* Create constant for converting angles to radians. */
    mpfr_t radian_conv;
    mpfr_init2(radian_conv, nbits);
    mpfr_const_pi(radian_conv, MPFR_RNDN);
    mpfr_div_si(radian_conv, radian_conv, 180, MPFR_RNDN);

    /* Set initial values converting angles to radians. */
    mpfr_set(t_curr, TMIN, MPFR_RNDN);
    mpfr_mul(yin.th1, TH10, radian_conv, MPFR_RNDN);  // th1[0] = TH10*PI/180.0;
    mpfr_mul(yin.w1, W10, radian_conv, MPFR_RNDN);    // w1[0] = W10*PI/180.0;
    mpfr_mul(yin.th2, TH20, radian_conv, MPFR_RNDN);  // th2[0] = TH20*PI/180.0;
    mpfr_mul(yin.w2, W20, radian_conv, MPFR_RNDN);    // w2[0] = W20*PI/180.0; 

    /* Clean up. */
    mpfr_clears(TMIN, TMAX, TH10, W10, TH20, W20, radian_conv, NULL);

    /* ********************** FOR LYAPUNOV EXP: ********************** */
    y_t yin_adj, yout_adj;
    mpfr_t d0, di, sum, delta, exp;
    mpfr_inits2(nbits, yin_adj.th1, yin_adj.w1, yin_adj.th2, yin_adj.w2, yout_adj.th1, yout_adj.w1, yout_adj.th2, yout_adj.w2, NULL);
    mpfr_inits2(53, d0, di, sum, delta, exp, NULL);
   
    mpfr_set_d(sum, 0.0, MPFR_RNDN);
    int n = -nbits;

    // calc d0; square root of machine precision
    mpfr_set_d(d0, 2.0, MPFR_RNDN);
    mpfr_pow_si(d0, d0, n, MPFR_RNDN);
    mpfr_sqrt(d0, d0, MPFR_RNDN);

    // pick a point that is d0 away from the initial condition
    mpfr_div_si(delta, d0, 2, MPFR_RNDN);

    mpfr_add(yin_adj.th1, yin.th1, delta, MPFR_RNDN);
    mpfr_add(yin_adj.w1, yin.w1, delta, MPFR_RNDN);
    mpfr_add(yin_adj.th2, yin.th2, delta, MPFR_RNDN);
    mpfr_add(yin_adj.w2, yin.w2, delta, MPFR_RNDN);
    
   
    double c;
    /* *************************************************************** */

    /* Output initial values. */
    output_polar(polar_output, &t_curr,  &yin.th1, &yin.w1, &yin.th2, &yin.w2);
    output_cartesian(cartesian_output, nbits, &t_curr, &yin.th1, &yin.w1, &yin.th2, &yin.w2, L1, L2);
    //output_energy(energy_output, nbits, &t_curr,  &yin.th1, &yin.w1, &yin.th2, &yin.w2, L1, L2, G);

    /* Perform the integration. */
    for (int i = 0; i < NSTEP; i++) {

      mpfr_add(t_next, t_curr, h, MPFR_RNDN);		// update time
      runge_kutta(t_curr, &yin, &yout, h);    		// preform runge kutta 
      runge_kutta(t_curr, &yin_adj, &yout_adj, h);	// preform runge kutta on adjusted IC

      calc_di(&di, &yout, &yout_adj);
      lyapunov(&sum, &d0, &di);
      reset_yin_adj(&d0, &di, &yout, &yout_adj, &yin_adj);

      /* Print output to files. */
      output_polar(polar_output, &t_next, &yout.th1, &yout.w1, &yout.th2, &yout.w2);
      output_cartesian(cartesian_output, nbits, &t_next, &yout.th1, &yout.w1, &yout.th2, &yout.w2, L1, L2);
      //output_energy(energy_output, nbits, &t_next, &yout.th1, &yout.w1, &yout.th2, &yout.w2, L1, L2, G);

      if (i != 0 && i%10 == 0) {
        c = 1.0/(double) i;
        mpfr_set_d(exp, c, MPFR_RNDN);
        mpfr_div(exp, exp, h, MPFR_RNDN);
        mpfr_mul(exp, exp, sum, MPFR_RNDN);

        output_lyapunov(ly_exp_output, &t_next, &exp);
      }

      /* Set yin to yout. */
      mpfr_set(yin.th1, yout.th1, MPFR_RNDN);
      mpfr_set(yin.w1, yout.w1, MPFR_RNDN);
      mpfr_set(yin.th2, yout.th2, MPFR_RNDN);
      mpfr_set(yin.w2, yout.w2, MPFR_RNDN);
      mpfr_set(t_curr, t_next, MPFR_RNDN);

      calc_di(&d0, &yin, &yin_adj);
      //if (mpfr_cmp(di, d0) != 0) {mpfr_printf("i=%d:  BAD IC d0 = %.24Rf  di = %.24Rf\n", i, d0, di); }
      //else { printf("GOOD IC i=%d\n",i); }
    }

    c = 1.0/((double) NSTEP);

    mpfr_set_d(exp, c, MPFR_RNDN);
    mpfr_div(exp, exp, h, MPFR_RNDN);
    mpfr_mul(exp, exp, sum, MPFR_RNDN);

    mpfr_printf("Lyapunov Exponent: %.24Rf\n", exp);
    mpfr_fprintf(final_exp_output, "%s,%s,%s,%s,%s,%s,%s,%0.32RNF\n", 
                argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], exp);

    /* Clean up. */
    mpfr_clears(yin_adj.th1, yin_adj.w1, yin_adj.th2, yin_adj.w2, 
  		yout_adj.th1, yout_adj.w1, yout_adj.th2, yout_adj.w2,
  		d0, di, sum, delta, exp, NULL);
    mpfr_clears(h, t_curr, t_next, yin.th1, yin.w1, yin.th2, yin.w2,
  		yout.th1, yout.w1, yout.th2, yout.w2, NULL);
    mpfr_free_cache();

    /* Close files. */
    fclose(polar_output);
    fclose(ly_exp_output);
    fclose(cartesian_output);
    //fclose(energy_output);
  }

  /* Close files. */
  fclose(final_exp_output);
  fclose(ic_record);

  return 0;
}






/* function to fill array of derivatives dydx at t */
void derivs(y_t *yin, y_t *dydx)
{
  mpfr_t num1, den1, num2, den2, temp1, temp2;
  mpfr_t const_mass_sum, del, cos_del, sin_del, sin_cos_del, sin_th1, sin_th2, w1_sqr, w2_sqr;

  // CALC: const_mass_sum = M1 + M2
  mpfr_init2(const_mass_sum, nbits);
  mpfr_set_d(const_mass_sum, M1, MPFR_RNDN);
  mpfr_add_d(const_mass_sum, const_mass_sum, M2, MPFR_RNDN);

  // INIT angle_t variables
  mpfr_inits2(nbits, del, cos_del, sin_del, sin_cos_del, sin_th1, sin_th2, NULL);

  // CALC: del = yin->th2 - yin->th1;
  mpfr_sub(del, yin->th2, yin->th1, MPFR_RNDN);

  // CALC: all sin and cos values
  mpfr_sin_cos(sin_del, cos_del, del, MPFR_RNDN);     // sin(del) & cos(del)
  mpfr_mul(sin_cos_del, sin_del, cos_del, MPFR_RNDN); // sin(del)*cos(del)
  mpfr_sin(sin_th1, yin->th1, MPFR_RNDN);             // sin(th_1)
  mpfr_sin(sin_th2, yin->th2, MPFR_RNDN);             // sin(th_2)

  // INIT velocity_t variables
  mpfr_inits2(nbits, w1_sqr, w2_sqr, NULL);

  // CALC: velocity_t variables;
  mpfr_sqr(w1_sqr, yin->w1, MPFR_RNDN);   // w1^2
  mpfr_sqr(w2_sqr, yin->w2, MPFR_RNDN);   // w2^2
  
  // INIT numerators, denomenators, and temps
  mpfr_inits2(nbits, num1, den1, num2, den2, temp1, temp2, NULL);

  // SET: dydx->th1 = yin->w1;
  mpfr_set(dydx->th1, yin->w1, MPFR_RNDN);

  // CALC: den1 = (M1+M2)*L1 - M2*L1*cos(del)*cos(del);
  mpfr_mul_d(den1, const_mass_sum, L1, MPFR_RNDN);
  mpfr_set_d(temp1, M2, MPFR_RNDN);
  mpfr_mul_d(temp1, temp1, L1, MPFR_RNDN);
  mpfr_mul(temp1, temp1, cos_del, MPFR_RNDN);
  mpfr_mul(temp1, temp1, cos_del, MPFR_RNDN);
  mpfr_sub(den1, den1, temp1, MPFR_RNDN);

  // CALC num1 = M2(L1*w1_sqr*sin_cos_del + G*sin_th2*cos_del + L2*w2_sqr*sin_del) 
  //            - (cons_mass_sum)*G*sin_th1
  mpfr_mul_d(num1, w1_sqr, L1, MPFR_RNDN);
  mpfr_mul(num1, num1, sin_cos_del, MPFR_RNDN); // num1 = L1*w1_sqr*sin_cos_del
  mpfr_mul_d(temp1, sin_th2, G, MPFR_RNDN);
  mpfr_mul(temp1, temp1, cos_del, MPFR_RNDN);   // temp1 = G*sin_th2*cos_del 
  mpfr_mul_d(temp2, w2_sqr, L2, MPFR_RNDN);
  mpfr_mul(temp2, temp2, sin_del, MPFR_RNDN);   // temp2 = L2*w2_sqr*sin_del
  mpfr_add(num1, num1, temp1, MPFR_RNDN);
  mpfr_add(num1, num1, temp2, MPFR_RNDN);
  mpfr_mul_d(num1, num1, M2, MPFR_RNDN);          // num1 = M2(num1 + temp1 + temp2)
  mpfr_mul_d(temp1, const_mass_sum, G, MPFR_RNDN);
  mpfr_mul(temp1, temp1, sin_th1, MPFR_RNDN);   // temp1 = (cons_mass_sum)*G*sin_th1)
  mpfr_sub(num1, num1, temp1, MPFR_RNDN); // num1 = num1 - temp1

  mpfr_div(dydx->w1, num1, den1, MPFR_RNDN);  // dydx->w1 = num1/den1

  // SET: dydx->th2 = yin->w2;
  mpfr_set(dydx->th2, yin->w2, MPFR_RNDN); 

  // CALC: den2 = (L2/L1)*den1;
  mpfr_set_d(den2, L2/L1, MPFR_RNDN);
  mpfr_mul(den2, den2, den1, MPFR_RNDN);

  // CALC: num2 = (M1+M2) (G*sin_th1*cos_del - L1*w1_sqr*sin_del - G*sin_th2)
  //            - M2*L2*w2_sqr*sin_cos_del
  mpfr_mul_d(num2, sin_th1, G, MPFR_RNDN);
  mpfr_mul(num2, num2, cos_del, MPFR_RNDN);   // num2 = G*sin_th1*cos_del
  mpfr_mul_d(temp1, w1_sqr, L1, MPFR_RNDN);
  mpfr_mul(temp1, temp1, sin_del, MPFR_RNDN); // temp1 = L1*w1_sqr*sin_del 
  mpfr_mul_d(temp2, sin_th2, G, MPFR_RNDN);     // temp2 = G*sin_th2
  mpfr_sub(num2, num2, temp1, MPFR_RNDN);
  mpfr_sub(num2, num2, temp2, MPFR_RNDN);
  mpfr_mul(num2, num2, const_mass_sum, MPFR_RNDN);  // num2 = (M1+M2)(num2 - temp1 - temp2)
  mpfr_set_d(temp1, M2 * L2, MPFR_RNDN);
  mpfr_mul(temp1, temp1, w2_sqr, MPFR_RNDN);
  mpfr_mul(temp1, temp1, sin_cos_del, MPFR_RNDN);   // temp1 = M2*L2*w2_sqr*sin_cos_del
  mpfr_sub(num2, num2, temp1, MPFR_RNDN);

  mpfr_div(dydx->w2, num2, den2, MPFR_RNDN);   // dydx->w2 = num2/den2 

  /*  Clean up. */
  mpfr_clears(num1, num2, den1, den2, temp1, temp2, const_mass_sum, del, cos_del, sin_del, 
              sin_cos_del, sin_th1, sin_th2, w1_sqr, w2_sqr, NULL);
  mpfr_free_cache();

  return;

}

void runge_kutta(mpfr_t t, y_t *yin, y_t *yout, mpfr_t h)
{
  /* fourth order Runge-Kutta - see e.g. Numerical Recipes */
 
  y_t dydx, dydxt, yt, k1, k2, k3, k4;
  mpfr_t temp1;
  mpfr_t THREE, SIX, ONE_HALF;

  mpfr_inits2(nbits, THREE, SIX, ONE_HALF, NULL);
  mpfr_set_d(THREE, 3.0, MPFR_RNDN);
  mpfr_set_d(SIX, 6.0, MPFR_RNDN);
  mpfr_set_d(ONE_HALF, 0.5, MPFR_RNDN);
  
  // INIT y_t struct contents
  mpfr_inits2(nbits, dydx.th1, dydx.w1, dydx.th2, dydx.w2, NULL);
  mpfr_inits2(nbits, dydxt.th1, dydxt.w1, dydxt.th2, dydxt.w2, NULL);
  mpfr_inits2(nbits, yt.th1, yt.w1, yt.th2, yt.w2, NULL);
  mpfr_inits2(nbits, k1.th1, k1.w1, k1.th2, k1.w2, NULL);
  mpfr_inits2(nbits, k2.th1, k2.w1, k2.th2, k2.w2, NULL);
  mpfr_inits2(nbits, k3.th1, k3.w1, k3.th2, k3.w2, NULL);
  mpfr_inits2(nbits, k4.th1, k4.w1, k4.th2, k4.w2, NULL);

  // INIT temps
  mpfr_init2(temp1, nbits);
  

  derivs(yin, &dydx); /* first step */

  mpfr_mul(k1.th1, h, dydx.th1, MPFR_RNDN);       // k1.th1 = h*dydx.th1;
  mpfr_mul(yt.th1, ONE_HALF, k1.th1, MPFR_RNDN);  
  mpfr_add(yt.th1, yt.th1, yin->th1, MPFR_RNDN);  // yt.th1 = yin->th1 + 0.5*k1.th1;

  mpfr_mul(k1.w1, h, dydx.w1, MPFR_RNDN);         // k1.w1 = h*dydx.w1;
  mpfr_mul(yt.w1, ONE_HALF, k1.w1, MPFR_RNDN);  
  mpfr_add(yt.w1, yt.w1, yin->w1, MPFR_RNDN);     // yt.w1 = yin->w1 + 0.5*k1.w1;

  mpfr_mul(k1.th2, h, dydx.th2, MPFR_RNDN);       // k1.th2 = h*dydx.th2;
  mpfr_mul(yt.th2, ONE_HALF, k1.th2, MPFR_RNDN);  
  mpfr_add(yt.th2, yt.th2, yin->th2, MPFR_RNDN);  // yt.th2 = yin->th2 + 0.5*k1.th2;

  mpfr_mul(k1.w2, h, dydx.w2, MPFR_RNDN);         // k1.w2 = h*dydx.w2;
  mpfr_mul(yt.w2, ONE_HALF, k1.w2, MPFR_RNDN);  
  mpfr_add(yt.w2, yt.w2, yin->w2, MPFR_RNDN);     // yt.w2 = yin->w2 + 0.5*k1.w2;
  

  derivs(&yt, &dydxt); /* second step */ 

  mpfr_mul(k2.th1, h, dydxt.th1, MPFR_RNDN);      // k2.th1 = h*dydxt.th1;
  mpfr_mul(yt.th1, ONE_HALF, k2.th1, MPFR_RNDN);  
  mpfr_add(yt.th1, yt.th1, yin->th1, MPFR_RNDN);  // yt.th1 = yin->th1 + 0.5*k2.th1;

  mpfr_mul(k2.w1, h, dydxt.w1, MPFR_RNDN);        // k2.w1 = h*dydxt.w1;
  mpfr_mul(yt.w1, ONE_HALF, k2.w1, MPFR_RNDN);  
  mpfr_add(yt.w1, yt.w1, yin->w1, MPFR_RNDN);     // yt.w1 = yin->w1 + 0.5*k2.w1;

  mpfr_mul(k2.th2, h, dydxt.th2, MPFR_RNDN);      // k2.th2 = h*dydxt.th2;
  mpfr_mul(yt.th2, ONE_HALF, k2.th2, MPFR_RNDN);  
  mpfr_add(yt.th2, yt.th2, yin->th2, MPFR_RNDN);  // yt.th2 = yin->th2 + 0.5*k2.th2;

  mpfr_mul(k2.w2, h, dydxt.w2, MPFR_RNDN);        // k2.w2 = h*dydxt.w2;
  mpfr_mul(yt.w2, ONE_HALF, k2.w2, MPFR_RNDN);  
  mpfr_add(yt.w2, yt.w2, yin->w2, MPFR_RNDN);     // yt.w2 = yin->w2 + 0.5*k2.w2;
  

  derivs(&yt, &dydxt); /* third step */

  mpfr_mul(k3.th1, h, dydxt.th1, MPFR_RNDN);      // k3.th1 = h*dydxt.th1;
  mpfr_add(yt.th1, yt.th1, yin->th1, MPFR_RNDN);  // yt.th1 = yin->th1 + k3.th1;

  mpfr_mul(k3.w1, h, dydxt.w1, MPFR_RNDN);        //  k3.w1 = h*dydxt.w1;
  mpfr_add(yt.w1, yt.w1, yin->w1, MPFR_RNDN);     // yt.w1 = yin->w1 + k3.w1;

  mpfr_mul(k3.th2, h, dydxt.th2, MPFR_RNDN);      // k3.th2 = h*dydxt.th2;
  mpfr_add(yt.th2, yt.th2, yin->th2, MPFR_RNDN);  // yt.th2 = yin->th2 + k3.th2;

  mpfr_mul(k3.w2, h, dydxt.w2, MPFR_RNDN);        // k3.w2 = h*dydxt.w2;
  mpfr_add(yt.w2, yt.w2, yin->w2, MPFR_RNDN);     // yt.w2 = yin->w2 + k3.w2;


  derivs(&yt, &dydxt); /* fourth step */

  mpfr_mul(k4.th1, h, dydxt.th1, MPFR_RNDN);  // k4.th1 = h*dydxt.th1;
  // CALC: yout->th1 = yin->th1 + k1.th1/6.0 + k2.th1/3.0 + k3.th1/3.0 + k4.th1/6.0;
  mpfr_set(yout->th1, yin->th1, MPFR_RNDN);   // yout->th1 = yin->th1
  mpfr_add(temp1, k1.th1, k4.th1, MPFR_RNDN);
  mpfr_div(temp1, temp1, SIX, MPFR_RNDN);     // temp1 = (k1.th1 + k4.th1)/6.0
  mpfr_add(yout->th1, yout->th1, temp1, MPFR_RNDN);
  mpfr_add(temp1, k2.th1, k3.th1, MPFR_RNDN);
  mpfr_div(temp1, temp1, THREE, MPFR_RNDN);     // temp1 = (k1.th1 + k4.th1)/3.0
  mpfr_add(yout->th1, yout->th1, temp1, MPFR_RNDN);

  mpfr_mul(k4.w1, h, dydxt.w1, MPFR_RNDN);  // k4.w1 = h*dydxt.w1;
  // CALC: yout->w1 = yin->w1 + k1.w1/6.0 + k2.w1/3.0 + k3.w1/3.0 + k4.w1/6.0;
  mpfr_set(yout->w1, yin->w1, MPFR_RNDN);   // yout->w1 = yin->w1
  mpfr_add(temp1, k1.w1, k4.w1, MPFR_RNDN);
  mpfr_div(temp1, temp1, SIX, MPFR_RNDN);   // temp1 = (k1.w1 + k4.w1)/6.0
  mpfr_add(yout->w1, yout->w1, temp1, MPFR_RNDN);
  mpfr_add(temp1, k2.w1, k3.w1, MPFR_RNDN);
  mpfr_div(temp1, temp1, THREE, MPFR_RNDN);   // temp1 = (k1.w1 + k4.w1)/3.0
  mpfr_add(yout->w1, yout->w1, temp1, MPFR_RNDN);

  mpfr_mul(k4.th2, h, dydxt.th2, MPFR_RNDN);  // k4.th2 = h*dydxt.th2;
  // CALC: yout->th2 = yin->th2 + k1.th2/6.0 + k2.th2/3.0 + k3.th2/3.0 + k4.th2/6.0;
  mpfr_set(yout->th2, yin->th2, MPFR_RNDN);   // yout->th2 = yin->th2
  mpfr_add(temp1, k1.th2, k4.th2, MPFR_RNDN);
  mpfr_div(temp1, temp1, SIX, MPFR_RNDN);     // temp1 = (k1.th2 + k4.th2)/6.0
  mpfr_add(yout->th2, yout->th2, temp1, MPFR_RNDN);
  mpfr_add(temp1, k2.th2, k3.th2, MPFR_RNDN);
  mpfr_div(temp1, temp1, THREE, MPFR_RNDN);     // temp1 = (k1.th2 + k4.th2)/3.0
  mpfr_add(yout->th2, yout->th2, temp1, MPFR_RNDN);

  mpfr_mul(k4.w2, h, dydxt.w2, MPFR_RNDN);  // k4.w2 = h*dydxt.w2;
  // CALC: yout->w2 = yin->w2 + k1.w2/6.0 + k2.w2/3.0 + k3.w2/3.0 + k4.w2/6.0;
  mpfr_set(yout->w2, yin->w2, MPFR_RNDN);   // yout->w2 = yin->w2
  mpfr_add(temp1, k1.w2, k4.w2, MPFR_RNDN);
  mpfr_div(temp1, temp1, SIX, MPFR_RNDN);     // temp1 = (k1.w2 + k4.w2)/6.0
  mpfr_add(yout->w2, yout->w2, temp1, MPFR_RNDN);
  mpfr_add(temp1, k2.w2, k3.w2, MPFR_RNDN);
  mpfr_div(temp1, temp1, THREE, MPFR_RNDN);     // temp1 = (k1.w2 + k4.w2)/3.0
  mpfr_add(yout->w2, yout->w2, temp1, MPFR_RNDN);

 /* Clean up. */
  mpfr_clears(dydx.th1, dydx.w1, dydx.th2, dydx.w2, 
              dydxt.th1, dydxt.w1, dydxt.th2, dydxt.w2, 
              yt.th1, yt.w1, yt.th2, yt.w2, 
              k1.th1, k1.w1, k1.th2, k1.w2,
              k2.th1, k2.w1, k2.th2, k2.w2,
              k3.th1, k3.w1, k3.th2, k3.w2, 
              k4.th1, k4.w1, k4.th2, k4.w2, 
              temp1, THREE, SIX, ONE_HALF, NULL);
  mpfr_free_cache();

  return;
}

void lyapunov(mpfr_t *sum, mpfr_t *d0, mpfr_t *di) {
  mpfr_t temp;
  mpfr_init2(temp, 53);

  mpfr_div(temp, *di, *d0, MPFR_RNDN);
  mpfr_log(temp, temp, MPFR_RNDN);

  mpfr_add(*sum, *sum, temp, MPFR_RNDN);

  mpfr_clears(temp, NULL);
}


void calc_di(mpfr_t *di, y_t *y0, y_t *y1) {
  mpfr_t temp;
  mpfr_init_set_d(temp, 0.0, 53);
  
  mpfr_sub(temp, y0->th1, y1->th1, MPFR_RNDN);
  mpfr_sqr(temp, temp, MPFR_RNDN);
  mpfr_set(*di, temp, MPFR_RNDN);

  mpfr_sub(temp, y0->w1, y1->w1, MPFR_RNDN);
  mpfr_sqr(temp, temp, MPFR_RNDN);
  mpfr_add(*di, *di, temp, MPFR_RNDN);

  mpfr_sub(temp, y0->th2, y1->th2, MPFR_RNDN);
  mpfr_sqr(temp, temp, MPFR_RNDN);
  mpfr_add(*di, *di, temp, MPFR_RNDN);

  mpfr_sub(temp, y0->w2, y1->w2, MPFR_RNDN);
  mpfr_sqr(temp, temp, MPFR_RNDN);
  mpfr_add(*di, *di, temp, MPFR_RNDN);

  mpfr_sqrt(*di, *di, MPFR_RNDN);

  mpfr_clear(temp);
}

void reset_yin_adj(mpfr_t *d0, mpfr_t *di, y_t *y0, y_t *y1_out, y_t *y1_in) {
  mpfr_t d;
  mpfr_init_set_d(d, 0.0, 53);
//mpfr_inits2(nbits, temp, d, NULL);

  mpfr_div(d, *d0, *di, MPFR_RNDN);
  
  mpfr_sub(y1_in->th1, y1_out->th1, y0->th1, MPFR_RNDN);
  mpfr_mul(y1_in->th1, y1_in->th1, d, MPFR_RNDN);
  mpfr_add(y1_in->th1, y1_in->th1, y0->th1, MPFR_RNDN);

  mpfr_sub(y1_in->w1, y1_out->w1, y0->w1, MPFR_RNDN);
  mpfr_mul(y1_in->w1, y1_in->w1, d, MPFR_RNDN);
  mpfr_add(y1_in->w1, y1_in->w1, y0->w1, MPFR_RNDN);

  mpfr_sub(y1_in->th2, y1_out->th2, y0->th2, MPFR_RNDN);
  mpfr_mul(y1_in->th2, y1_in->th2, d, MPFR_RNDN);
  mpfr_add(y1_in->th2, y1_in->th2, y0->th2, MPFR_RNDN);

  mpfr_sub(y1_in->w2, y1_out->w2, y0->w2, MPFR_RNDN);
  mpfr_mul(y1_in->w2, y1_in->w2, d, MPFR_RNDN);
  mpfr_add(y1_in->w2, y1_in->w2, y0->w2, MPFR_RNDN);

  mpfr_clears(d, NULL);
}























