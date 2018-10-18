/* 
 * E dpend_mpfr.c
 *
 * Example code to solve double pendulum ODEs using fourth order 
 * Runge-Kutta. 
 *
 * Parameters are passed in at the command line:
 * 
 * $./solve_dpend TMIN TMAX TH10 W10 TH20 W20 NSTEP > pendulum.txt 
 *
 * where TMIN and TMAX are the starting and ending times (in seconds), 
 * TH10 and TH20 are the initial angles (degrees), and W10 and W20 
 * are the initial angular velocities (degrees per second), and 
 * NSTEP is the number of integrations steps. This example illustrates 
 * using redirection to write the results to file in a file
 * pendulum.txt. Note that there is no checking for accuracy, so the
 * user needs to choose a suitable NSTEP. Also angles written to file
 * are in radians.
 *
 * As an example, the data for the first animated gif on the web page 
 * may be generated with   
 *
 * $./solve_dpend 0.0 10.0 90.0 0.00 -10.0 0.0 1000 > outfile.txt
 *
 * (only every fifth frame is used in the animation).
 * 
 * M.S. Wheatland, 2004
 *
 * ------------------------------------------------------------------
 * Edited by Becky Wild to support generic floating-point types.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>

/* hardwired parameters */

#define PI 3.14159265
#define N 4 /* number of equations to be solved */
#define G 9.8 /* acc'n due to gravity, in m/s^2 */
#define L1 1.0 /* length of pendulum 1 in m */
#define L2 1.0 /* length of pendulum 2 in m */
#define M1 1.0 /* mass of pendulum 1 in kg */
#define M2 1.0 /* mass of pendulum 2 in kg */

size_t nbits;

typedef mpfr_t seconds_t;  // seconds_t because linux has its own time_t
typedef mpfr_t angle_t;
typedef mpfr_t velocity_t;

typedef struct {
  angle_t th1;
  velocity_t w1;
  angle_t th2;
  velocity_t w2;
} y_t;

void runge_kutta(seconds_t t, y_t *yin, y_t *yout, seconds_t h);
void derivs(y_t *yin, y_t *dydx);

int main(int argc, char *argv[])
{
  unsigned int i = 0, NSTEP;

  /* Get percision from args before any initilizing. */
  nbits = atoi(argv[8]);

  seconds_t h, TMIN, TMAX, t_curr, t_next;
  mpfr_inits2(nbits, h, TMIN, TMAX, t_curr, t_next, NULL);

  angle_t TH10, TH20;
  mpfr_inits2(nbits, TH10, TH20, NULL);

  velocity_t W10, W20;
  mpfr_inits2(nbits, W10, W20, NULL);

  y_t yin, yout;
  mpfr_inits2(nbits, yin.th1, yin.w1, yin.th2, yin.w2, yout.th1, yout.w1, yout.th2, yout.w2, NULL);

  mpfr_t x1, y1, x2, y2, temp;
  mpfr_inits2(nbits, x1, y1, x2, y2, temp, NULL);

  /* obtain command line values */
  mpfr_set_flt(TMIN, atof(argv[1]), MPFR_RNDN);
  mpfr_set_flt(TMAX, atof(argv[2]), MPFR_RNDN);
  mpfr_set_flt(TH10, atof(argv[3]), MPFR_RNDN);
  mpfr_set_flt(W10, atof(argv[4]), MPFR_RNDN);
  mpfr_set_flt(TH20, atof(argv[5]), MPFR_RNDN);
  mpfr_set_flt(W20, atof(argv[6]), MPFR_RNDN);
  NSTEP = atoi(argv[7]);

  if (NSTEP < 2) {
    printf("Number of steps must be greater than 1");
    return 0;
  }

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

  /* Get initial coordinates. */
  mpfr_set_d(x1, L1, MPFR_RNDN);        // calc x1
  mpfr_sin(temp, yin.th1, MPFR_RNDN);
  mpfr_mul(x1, x1, temp, MPFR_RNDN);

  mpfr_set_d(y1, -L1, MPFR_RNDN);       // calc y1
  mpfr_cos(temp, yin.th1, MPFR_RNDN);
  mpfr_mul(y1, y1, temp, MPFR_RNDN);

  mpfr_set_d(x2, L2, MPFR_RNDN);        // calc x2
  mpfr_sin(temp, yin.th2, MPFR_RNDN);
  mpfr_mul(x2, x2, temp, MPFR_RNDN);
  mpfr_add(x2, x2, x1, MPFR_RNDN);

  mpfr_set_d(y2, -L2, MPFR_RNDN);        // calc y2
  mpfr_cos(temp, yin.th2, MPFR_RNDN);
  mpfr_mul(y2, y2, temp, MPFR_RNDN);
  mpfr_add(y2, y2, y1, MPFR_RNDN);


  /* Print initial values. */
 // mpfr_printf("%0.6RNf %0.6RNf %0.6RNF %0.6RNF %0.6RNF\n", t_curr, yin.th1, yin.w1, yin.th2, yin.w2);
  mpfr_printf("%0.6RNf %0.6RNf %0.6RNF %0.6RNF %0.6RNF\n", t_curr, x1, y1, x2, y2);

  /* perform the integration */
  for (i = 0; i < NSTEP - 1; i++)
  { 
    mpfr_add(t_next, t_curr, h, MPFR_RNDN); // update time
    runge_kutta(t_curr, &yin, &yout, h);    // preform runge kutta

    /* Get coordinates. */
    mpfr_set_d(x1, L1, MPFR_RNDN);        // calc x1
    mpfr_sin(temp, yout.th1, MPFR_RNDN);
    mpfr_mul(x1, x1, temp, MPFR_RNDN);

    mpfr_set_d(y1, -L1, MPFR_RNDN);       // calc y1
    mpfr_cos(temp, yout.th1, MPFR_RNDN);
    mpfr_mul(y1, y1, temp, MPFR_RNDN);

    mpfr_set_d(x2, L2, MPFR_RNDN);        // calc x2
    mpfr_sin(temp, yout.th2, MPFR_RNDN);
    mpfr_mul(x2, x2, temp, MPFR_RNDN);
    mpfr_add(x2, x2, x1, MPFR_RNDN);

    mpfr_set_d(y2, -L2, MPFR_RNDN);        // calc y2
    mpfr_cos(temp, yout.th2, MPFR_RNDN);
    mpfr_mul(y2, y2, temp, MPFR_RNDN);
    mpfr_add(y2, y2, y1, MPFR_RNDN);

    /* Print "t th1 w1 th2 w2" */
    //mpfr_printf("%0.6RNf %0.6RNf %0.6RNF %0.6RNF %0.6RNF\n", 
    //            t_next, yout.th1, yout.w1, yout.th2, yout.w2);  
    mpfr_printf("%0.6RNf %0.6RNf %0.6RNF %0.6RNF %0.6RNF\n", t_curr, x1, y1, x2, y2);  

    /* Set yin to yout. */
    mpfr_set(yin.th1, yout.th1, MPFR_RNDN);
    mpfr_set(yin.w1, yout.w1, MPFR_RNDN);
    mpfr_set(yin.th2, yout.th2, MPFR_RNDN);
    mpfr_set(yin.w2, yout.w2, MPFR_RNDN);
  
    mpfr_set(t_curr, t_next, MPFR_RNDN);
  }

  /* Clean up. */
  mpfr_clears(h, t_curr, t_next, yin.th1, yin.w1, yin.th2, yin.w2, yout.th1, yout.w1, yout.th2, yout.w2, NULL);
  mpfr_free_cache();

  return 0;
}

/* function to fill array of derivatives dydx at t */
void derivs(y_t *yin, y_t *dydx)
{
  mpfr_t num1, den1, num2, den2, temp1, temp2;
  mpfr_t const_mass_sum;
  angle_t del, cos_del, sin_del, sin_cos_del, sin_th1, sin_th2;
  velocity_t w1_sqr, w2_sqr;

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

void runge_kutta(seconds_t t, y_t *yin, y_t *yout, seconds_t h)
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



