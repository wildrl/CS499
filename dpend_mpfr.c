/* solve_dpend.c
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
  float in_TMIN, in_TMAX, in_TH10, in_W10, in_TH20, in_W20;
  unsigned int i = 0, in_NSTEP, NSTEP;

  seconds_t h, TMIN, TMAX;
  mpfr_inits2(200, h, TMIN, TMAX);

  angle_t TH10, TH20;
  mpfr_inits2(200, TH10, TH20);

  velocity_t W10, W20;
  mpfr_inits2(200, W10, W20);

  y_t yin, yout;
  mpfr_inits2(200, yin.th1, yin.w1, yin.th2, yin.w2, yout.th1, yout.w1, yout.th2, yout.w2);

  seconds_t *t;
  angle_t *th1, *th2;
  velocity_t *w1, *w2;

  /* obtain command line values */

  mpfr_set(TMIN, atof(argv[1]), MPFR_RNDN);
  mpfr_set(TMAX, atof(argv[2]), MPFR_RNDN);
  mpfr_set(TH10, atof(argv[3]), MPFR_RNDN);
  mpfr_set(W10, atof(argv[4]), MPFR_RNDN);
  mpfr_set(TH20, atof(argv[5]), MPFR_RNDN);
  mpfr_set(W20, atof(argv[6]), MPFR_RNDN);

  NSTEP = atoi(argv[7]);

  /* allocate memory for arrays of values of time, angles 1 and 2,
     and angular velocities 1 and 2 respectively */ 

  t = (seconds_t *) malloc(NSTEP*sizeof(seconds_t)); 
  th1 = (angle_t *) malloc(NSTEP*sizeof(angle_t)); 
  w1 = (velocity_t *) malloc(NSTEP*sizeof(velocity_t));
  th2 = (angle_t *) malloc(NSTEP*sizeof(angle_t));
  w2 = (velocity_t *) malloc(NSTEP*sizeof(velocity_t));

  /* stepsize for integration */

  // CALC h: h = (TMAX - TMIN)/(NSTEP - 1.0);
  mpfr_set(h, TMAX, MPFR_RNDN);
  mpfr_sub(h, h, TMIN, MPFR_RNDN);
  NSTEP--;
  mpfr_div_ui(h, h, NSTEP, MPFR_RNDN);
  NSTEP++;
 
  /* Define array of t values */
  seconds_t h_temp;
  mpfr_init2(h_temp, 200);

  for (i = 0; i < NSTEP; i++) {
    // CALC t: t[i] = TMIN + h*i;
    mpfr_set(h_temp, h, MPFR_RNDN);             // h_temp = h
    mpfr_mul_ui(h_temp, h_temp, i, MPFR_RNDN);  // h_temp *= i
    mpfr_set(t[i], TMIN, MPFR_RNDN);            // t = TMIN
    mpfr_add(t[i], t[i], h_temp, MPFR_RNDN);          // t += h_temp
  }

  /* initial values - convert all angles to radians */

  // CALC radian_conv: PI/180
  mpfr_t radian_conv;
  mpfr_init2(radian_conv, 200);
  mpfr_set(radian_conv, 180, MPFR_RNDN);
  mpfr_mul(radian_conv, radian_conv, mpfr_const_pi, MPFR_RNDN);

  mpfr_set(th1[0], TH10, MPFR_RNDN);  // th1[0] = TH10*PI/180.0;
  mpfr_set(w1[0], W10, MPFR_RNDN);    // w1[0] = W10*PI/180.0;
  mpfr_set(th2[0], TH20, MPFR_RNDN);  // th2[0] = TH20*PI/180.0;
  mpfr_set(w2[0], W20, MPFR_RNDN);    // w2[0] = W20*PI/180.0; 

  /* perform the integration */

  mpfr_out_str (stdout, 10, 0, t[0], MPFR_RNDD);    putchar(' ');
  mpfr_out_str (stdout, 10, 0, th1[0], MPFR_RNDD);  putchar(' ');
  mpfr_out_str (stdout, 10, 0, w1[0], MPFR_RNDD);   putchar(' ');
  mpfr_out_str (stdout, 10, 0, th2[0], MPFR_RNDD);  putchar(' ');
  mpfr_out_str (stdout, 10, 0, w2[0], MPFR_RNDD);   putchar(' ');   putchar ('\n');

  for (i = 0; i < NSTEP - 1; i++)
  { 
    mpfr_set(yin.th1, th1[i], MPFR_RNDN);   // yin.th1 = th1[i];
    mpfr_set(yin.w1, w1[i], MPFR_RNDN);     // yin.w1 = w1[i];
    mpfr_set(yin.th2, th2[i], MPFR_RNDN);   // yin.th2 = th2[i];
    mpfr_set(yin.w2, w2[i], MPFR_RNDN);     // yin.w2 = w2[i];

    runge_kutta(t[i], &yin, &yout, h);

    mpfr_set(th1[i+1], yout.th1, MPFR_RNDN); // th1[i+1] = yout.th1;
    mpfr_set(w1[i+1], yout.w1, MPFR_RNDN);   // w1[i+1] = yout.w1;
    mpfr_set(th2[i+1], yout.th2, MPFR_RNDN); // th2[i+1] = yout.th2;
    mpfr_set(w2[i+1], yout.w2, MPFR_RNDN);   // w2[i+1] = yout.w2;

    mpfr_out_str (stdout, 10, 0, t[i+1], MPFR_RNDD);    putchar(' ');
    mpfr_out_str (stdout, 10, 0, th1[i+1], MPFR_RNDD);  putchar(' ');
    mpfr_out_str (stdout, 10, 0, w1[i+1], MPFR_RNDD);   putchar(' ');
    mpfr_out_str (stdout, 10, 0, th2[i+1], MPFR_RNDD);  putchar(' ');
    mpfr_out_str (stdout, 10, 0, w2[i+1], MPFR_RNDD);   putchar(' ');   putchar ('\n');
  }

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
  mpfr_init2(const_mass_sum, 200);
  mpfr_set(const_mass_sum, M1, MPFR_RNDN);
  mpfr_add(const_mass_sum, const_mass_sum, M2, MPFR_RNDN);

  // INIT angle_t variables
  mpfr_inits2(200, del, cos_del, sin_del, sin_cos_del, sin_th1, sin_th2);

  // CALC: del = yin->th2 - yin->th1;
  mpfr_sub(del, yin->th2, yin->th1, MPFR_RNDN);

  // CALC: all sin and cos values
  mpfr_sin_cos(sin_del, cos_del, del, MPFR_RNDN);     // sin(del) & cos(del)
  mpfr_mul(sin_cos_del, sin_del, cos_del, MPFR_RNDN); // sin(del)*cos(del)
  mpfr_sin(sin_th1, yin->th1, MPFR_RNDN);             // sin(th_1)
  mpfr_sin(sin_th2, yin->th2, MPFR_RNDN);             // sin(th_2)

  // INIT velocity_t variables
  mpfr_inits2(200, w1_sqr, w2_sqr);

  // CALC: velocity_t variables;
  mpfr_sqr(w1_sqr, yin->w1, MPFR_RNDN);   // w1^2
  mpfr_sqr(w2_sqr, yin->w2, MPFR_RNDN);   // w2^2
  
  // INIT numerators, denomenators, and temps
  mpfr_inits2(200, num1, den1, num2, den2, temp1, temp2);

  // SET: dydx->th1 = yin->w1;
  set(dydx->th1, yin->w1, MPFR_RNDN); 

  // CALC: den1 = (M1+M2)*L1 - M2*L1*cos(del)*cos(del);
  mpfr_mul(den1, const_mass_sum, L1, MPFR_RNDN);
  mpfr_set(temp1, M2, MPFR_RNDN);
  mpfr_mul(temp1, temp1, L1, MPFR_RNDN);
  mpfr_mul(temp1, temp1, cos_del, MPFR_RNDN);
  mpfr_mul(temp1, temp1, cos_del, MPFR_RNDN);
  mpfr_sub(den1, den1, temp1, MPFR_RNDN);

  // CALC num1 = M2(L1*w1_sqr*sin_cos_del + G*sin_th2*cos_del + L2*w2_sqr*sin_del) 
  //            - (cons_mass_sum)*G*sin_th1
  mpfr_mul(num1, L1, w1_sqr, MPFR_RNDN);
  mpfr_mul(num1, num1, sin_cos_del, MPFR_RNDN); // num1 = L1*w1_sqr*sin_cos_del
  mpfr_mul(temp1, G, sin_th2, MPFR_RNDN);
  mpfr_mul(temp1, temp1, cos_del, MPFR_RNDN);   // temp1 = G*sin_th2*cos_del 
  mpfr_mul(temp2, L2, w2_sqr, MPFR_RNDN);
  mpfr_mul(temp2, temp2, sin_del, MPFR_RNDN);   // temp2 = L2*w2_sqr*sin_del
  mpfr_add(num1, num1, temp1, MPFR_RNDN);
  mpfr_add(num1, num1, temp2, MPFR_RNDN);
  mpfr_mul(num1, num1, M2, MPFR_RNDN);          // num1 = M2(num1 + temp1 + temp2)
  mpfr_mul(temp1, const_mass_sum, G, MPFR_RNDN);
  mpfr_mul(temp1, temp1, sin_th1, MPFR_RNDN);   // temp1 = (cons_mass_sum)*G*sin_th1)
  mpfr_sub(num1, num1, temp1, MPFR_RNDN); // num1 = num1 - temp1

  mpfr_div(dydx->w1, num1, den1, MPFR_RNDN);  // dydx->w1 = num1/den1

  // SET: dydx->th2 = yin->w2;
  set(dydx->th2, yin->w2, MPFR_RNDN); 

  // CALC: den2 = (L2/L1)*den1;
  mpfr_div(den2, L2, L1, MPFR_RNDN);
  mpfr_mul(den2, den2, den1, MPFR_RNDN);

  // CALC: num2 = (M1+M2) (G*sin_th1*cos_del - L1*w1_sqr*sin_del - G*sin_th2)
  //            - M2*L2*w2_sqr*sin_cos_del
  mpfr_mul(num2, G, sin_th1, MPFR_RNDN);
  mpfr_mul(num2, num2, cos_del, MPFR_RNDN);   // num2 = G*sin_th1*cos_del
  mpfr_mul(temp1, L1, w1_sqr, MPFR_RNDN);
  mpfr_mul(temp1, temp1, sin_del, MPFR_RNDN); // temp1 = L1*w1_sqr*sin_del 
  mpfr_mul(temp2, G, sin_th2, MPFR_RNDN);     // temp2 = G*sin_th2
  mpfr_sub(num2, num2, temp1, MPFR_RNDN);
  mpfr_sub(num2, num2, temp2, MPFR_RNDN);
  mpfr_mul(num2, num2, const_mass_sum, MPFR_RNDN);  // num2 = (M1+M2)(num2 - temp1 - temp2)
  mpfr_mul(temp1, M2, L2, MPFR_RNDN);
  mpfr_mul(temp1, temp1, w2_sqr, MPFR_RNDN);
  mpfr_mul(temp1, temp1, sin_cos_del, MPFR_RNDN);   // temp1 = M2*L2*w2_sqr*sin_cos_del
  mpfr_sub(num2, num2, temp1, MPFR_RNDN);

  mpfr_div(dydx->w2, num2, den2, MPFR_RNDN);   // dydx->w2 = num2/den2


  return;

}

void runge_kutta(seconds_t t, y_t *yin, y_t *yout, seconds_t h)
{
  /* fourth order Runge-Kutta - see e.g. Numerical Recipes */
 
  y_t dydx, dydxt, yt, k1, k2, k3, k4;
  mpfr_t temp1;
  mpfr_t THREE, SIX, ONE_HALF;

  mpfr_inits2(200, THREE, SIX, ONE_HALF);
  mpfr_set(THREE, 3.0, MPFR_RNDN);
  mpfr_set(SIX, 6.0, MPFR_RNDN);
  mpfr_set(ONE_HALF, 0.5, MPFR_RNDN);
  
  // INIT y_t struct contents
  mpfr_inits2(200, dydx.th1, dydx.w2, dydx.th2, dydx.w2);
  mpfr_inits2(200, dydxt.th1, dydxt.w2, dydxt.th2, dydxt.w2);
  mpfr_inits2(200, yt.th1, yt.w2, yt.th2, yt.w2);
  mpfr_inits2(200, k1.th1, k1.w2, k1.th2, k1.w2);
  mpfr_inits2(200, k2.th1, k2.w2, k2.th2, k2.w2);
  mpfr_inits2(200, k3.th1, k3.w2, k3.th2, k3.w2);

  // INIT temps
  mpfr_init2(temp1, 200);
  
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

  mpfr_mul(k4.w1, h, dydxt.w1, MPFR_RNDN);  // k4.w2 = h*dydxt.w2;
  // CALC: yout->w2 = yin->w2 + k1.w2/6.0 + k2.w2/3.0 + k3.w2/3.0 + k4.w2/6.0;
  mpfr_set(yout->w2, yin->w2, MPFR_RNDN);   // yout->w2 = yin->w2
  mpfr_add(temp1, k1.w2, k4.w2, MPFR_RNDN);
  mpfr_div(temp1, temp1, SIX, MPFR_RNDN);     // temp1 = (k1.w2 + k4.w2)/6.0
  mpfr_add(yout->th2, yout->th2, temp1, MPFR_RNDN);
  mpfr_add(temp1, k2.w2, k3.w2, MPFR_RNDN);
  mpfr_div(temp1, temp1, THREE, MPFR_RNDN);     // temp1 = (k1.w2 + k4.w2)/3.0
  mpfr_add(yout->w2, yout->w2, temp1, MPFR_RNDN);
  
 
  return;
}

