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

/* hardwired parameters */

#define PI 3.14159265
#define N 4 /* number of equations to be solved */
#define G 9.8 /* acc'n due to gravity, in m/s^2 */
#define L1 1.0 /* length of pendulum 1 in m */
#define L2 1.0 /* length of pendulum 2 in m */
#define M1 1.0 /* mass of pendulum 1 in kg */
#define M2 1.0 /* mass of pendulum 2 in kg */

typedef double tyme_t;  // tyme_t because linux has its own time_t
typedef double angle_t;
typedef double velocity_t;

typedef struct {
  angle_t th1;
  velocity_t w1;
  angle_t th2;
  velocity_t w2;
} y_t;

void runge_kutta(tyme_t t, y_t *yin, y_t *yout, tyme_t h);
void derivs(y_t *yin, y_t *dydx);

int main(int argc, char *argv[])
{

  int i = 0, NSTEP;
  tyme_t h;
  tyme_t TMIN;
  tyme_t TMAX;
  angle_t TH10, TH20;
  velocity_t W10, W20;

  y_t yin, yout;
  tyme_t *t;
  angle_t *th1, *th2;
  velocity_t *w1, *w2;

  /* obtain command line values */

  TMIN = atof(argv[1]);
  TMAX = atof(argv[2]);
  TH10 = atof(argv[3]);
  W10 = atof(argv[4]);
  TH20 = atof(argv[5]);
  W20 = atof(argv[6]);
  NSTEP = atoi(argv[7]);

  /* allocate memory for arrays of values of time, angles 1 and 2,
     and angular velocities 1 and 2 respectively */ 

  t = (tyme_t *) malloc(NSTEP*sizeof(tyme_t)); 
  th1 = (angle_t *) malloc(NSTEP*sizeof(angle_t)); 
  w1 = (velocity_t *) malloc(NSTEP*sizeof(velocity_t));
  th2 = (angle_t *) malloc(NSTEP*sizeof(angle_t));
  w2 = (velocity_t *) malloc(NSTEP*sizeof(velocity_t));

  /* stepsize for integration */

  h = (TMAX - TMIN)/(NSTEP - 1.0);
 
  /* Define array of t values */

  for (i = 0; i < NSTEP; i++)
    t[i] = TMIN + h*i;

  /* initial values - convert all angles to radians */

  th1[0] = TH10*PI/180.0;
  w1[0] = W10*PI/180.0;
  th2[0] = TH20*PI/180.0;
  w2[0] = W20*PI/180.0; 

  /* perform the integration */

  printf("%f %f %f %f %f\n", t[0], th1[0], w1[0], th2[0], w2[0]);
  for (i = 0; i < NSTEP - 1; i++)
  { 
    yin.th1 = th1[i];
    yin.w1 = w1[i];
    yin.th2 = th2[i];
    yin.w2 = w2[i];
    runge_kutta(t[i], &yin, &yout, h);
    th1[i+1] = yout.th1;
    w1[i+1] = yout.w1;
    th2[i+1] = yout.th2;
    w2[i+1] = yout.w2;
  
    printf("%f %f %f %f %f\n", t[i+1], th1[i+1], w1[i+1], th2[i+1],
      w2[i+1]);
  }

  return 0;
 
}

void derivs(y_t *yin, y_t *dydx)
{

  /* function to fill array of derivatives dydx at t */

  float den1, den2, del;

  dydx->th1 = yin->w1; 
  
  del = yin->th2 - yin->th1;
  den1 = (M1+M2)*L1 - M2*L1*cos(del)*cos(del);
  dydx->w1 = (M2*L1*yin->w1*yin->w1*sin(del)*cos(del)
    + M2*G*sin(yin->th2)*cos(del) + M2*L2*yin->w2*yin->w2*sin(del)
    - (M1+M2)*G*sin(yin->th1))/den1;

  dydx->th2 = yin->w2;

  den2 = (L2/L1)*den1;
  dydx->w2 = (-M2*L2*yin->w2*yin->w2*sin(del)*cos(del)
    + (M1+M2)*G*sin(yin->th1)*cos(del) 
    - (M1+M2)*L1*yin->w1*yin->w1*sin(del)
    - (M1+M2)*G*sin(yin->th2))/den2;

  return;

}

void runge_kutta(tyme_t t, y_t *yin, y_t *yout, tyme_t h)
{
  /* fourth order Runge-Kutta - see e.g. Numerical Recipes */
 
  int i;
  y_t dydx, dydxt, yt, k1, k2, k3, k4; 
  
  derivs(yin, &dydx); /* first step */

  k1.th1 = h*dydx.th1;
  yt.th1 = yin->th1 + 0.5*k1.th1;
  k1.w1 = h*dydx.w1;
  yt.w1 = yin->w1 + 0.5*k1.w1;
  k1.th2 = h*dydx.th2;
  yt.th2 = yin->th2 + 0.5*k1.th2;
  k1.w2 = h*dydx.w2;
  yt.w2 = yin->w2 + 0.5*k1.w2;
  

  derivs(&yt, &dydxt); /* second step */ 

  k2.th1 = h*dydxt.th1;
  yt.th1 = yin->th1 + 0.5*k2.th1;
  k2.w1 = h*dydxt.w1;
  yt.w1 = yin->w1 + 0.5*k2.w1;
  k2.th2 = h*dydxt.th2;
  yt.th2 = yin->th2 + 0.5*k2.th2;
  k2.w2 = h*dydxt.w2;
  yt.w2 = yin->w2 + 0.5*k2.w2;
  

  derivs(&yt, &dydxt); /* third step */

  k3.th1 = h*dydxt.th1;
  yt.th1 = yin->th1 + k3.th1;
  k3.w1 = h*dydxt.w1;
  yt.w1 = yin->w1 + k3.w1;
  k3.th2 = h*dydxt.th2;
  yt.th2 = yin->th2 + k3.th2;
  k3.w2 = h*dydxt.w2;
  yt.w2 = yin->w2 + k3.w2;


  derivs(&yt, &dydxt); /* fourth step */

  k4.th1 = h*dydxt.th1;
  yout->th1 = yin->th1 + k1.th1/6.0 + k2.th1/3.0 + k3.th1/3.0 + k4.th1/6.0;
  k4.w1 = h*dydxt.w1;
  yout->w1 = yin->w1 + k1.w1/6.0 + k2.w1/3.0 + k3.w1/3.0 + k4.w1/6.0;
  k4.th2 = h*dydxt.th2;
  yout->th2 = yin->th2 + k1.th2/6.0 + k2.th2/3.0 + k3.th2/3.0 + k4.th2/6.0;
  k4.w2 = h*dydxt.w2;
  yout->w2 = yin->w2 + k1.w2/6.0 + k2.w2/3.0 + k3.w2/3.0 + k4.w2/6.0;
  
 
  return;

}
