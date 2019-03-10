#include "dpend.h"

/* 
 * Calculates the derivitive of each component in yin.
 * Result is stored in dydx.
 */
void derivs(y_t *yin, y_t *dydx) {

  /* Declare variables. */
  mpfr_t num1, den1;        // numerator and denomenator value for w1'
  mpfr_t num2, den2;        // numerator and denomenator value for w2'
  mpfr_t mass_sum;          // sum of pendulum masses, M1 + M2
  mpfr_t del;               // difference of pendulum angles, th2 - th1
  mpfr_t cos_del, sin_del;  // cos(del), sin(del)
  mpfr_t sin_th1, sin_th2;  // sin(th1), sin(th2)
  mpfr_t sin_cos_del;       // sin(del) * cos(del)
  mpfr_t w1_sqr, w2_sqr;    // w1^2, w2^2
  mpfr_t aux;               // auxilary variables used for calculations


  /* Initilize variable precision. */
  mpfr_inits2(nbits, mass_sum, w1_sqr, w2_sqr, num1, den1, num2, den2, 
    aux, del, cos_del, sin_del, sin_cos_del, sin_th1, sin_th2, NULL);


  /* Initilize variables. */
  mpfr_set_d(mass_sum, M1 + M2, MPFR_RNDN);           // mass_sum = M1 + M2
  mpfr_sub(del, yin->th2, yin->th1, MPFR_RNDN);       // del = th2 - th1;

  mpfr_sin_cos(sin_del, cos_del, del, MPFR_RNDN);     // sin_del = sin(del), cos_del = cos(del)
  mpfr_mul(sin_cos_del, sin_del, cos_del, MPFR_RNDN); // sin_cos_del = sin(del) * cos(del)
  mpfr_sin(sin_th1, yin->th1, MPFR_RNDN);             // sin_th1 = sin(th1)
  mpfr_sin(sin_th2, yin->th2, MPFR_RNDN);             // sin_th2 = sin(th2)

  mpfr_sqr(w1_sqr, yin->w1, MPFR_RNDN);               // w1_sqr = w1^2
  mpfr_sqr(w2_sqr, yin->w2, MPFR_RNDN);               // w2_sqr = w2^2


  /* Calculate th1' and th2'. */
  mpfr_set(dydx->th1, yin->w1, MPFR_RNDN);  // th1' = w1
  mpfr_set(dydx->th2, yin->w2, MPFR_RNDN);  // th2' = w2


  /* Calculate w1'. */
  // den1 = (M1+M2)*L1 - M2*L1*cos(del)^2;
  mpfr_mul_d(den1, mass_sum, L1, MPFR_RNDN);
  mpfr_set_d(aux, M2, MPFR_RNDN);
  mpfr_mul_d(aux, aux, L1, MPFR_RNDN);
  mpfr_mul(aux, aux, cos_del, MPFR_RNDN);
  mpfr_mul(aux, aux, cos_del, MPFR_RNDN);
  mpfr_sub(den1, den1, aux, MPFR_RNDN);

  // num1 = M2( L1*w1_sqr*sin_cos_del + G*sin_th2*cos_del + L2*w2_sqr*sin_del )  - (mass_sum)*G*sin_th1
  mpfr_mul_d(num1, w1_sqr, L1, MPFR_RNDN);
  mpfr_mul(num1, num1, sin_cos_del, MPFR_RNDN); // num1 = L1*w1_sqr*sin_cos_del
  mpfr_mul_d(aux, sin_th2, G, MPFR_RNDN);
  mpfr_mul(aux, aux, cos_del, MPFR_RNDN);       // aux = G*sin_th2*cos_del 
  mpfr_add(num1, num1, aux, MPFR_RNDN);         // num1 += aux
  mpfr_mul_d(aux, w2_sqr, L2, MPFR_RNDN);
  mpfr_mul(aux, aux, sin_del, MPFR_RNDN);       // aux = L2*w2_sqr*sin_del
  mpfr_add(num1, num1, aux, MPFR_RNDN);         // num1 += aux
  mpfr_mul_d(num1, num1, M2, MPFR_RNDN);        // num1 *= M2
  mpfr_mul_d(aux, mass_sum, G, MPFR_RNDN);
  mpfr_mul(aux, aux, sin_th1, MPFR_RNDN);       // aux = (mass_sum)*G*sin_th1)
  mpfr_sub(num1, num1, aux, MPFR_RNDN);         // num1 = num1 - aux

  mpfr_div(dydx->w1, num1, den1, MPFR_RNDN);    // w1' = num1/den1

  
  /* Calculate w2'. */
  // den2 = (M1+M2)*L2 - M2*L2*cos(del)^2;
  mpfr_mul_d(den2, mass_sum, L2, MPFR_RNDN);
  mpfr_set_d(aux, M2, MPFR_RNDN);
  mpfr_mul_d(aux, aux, L2, MPFR_RNDN);
  mpfr_mul(aux, aux, cos_del, MPFR_RNDN);
  mpfr_mul(aux, aux, cos_del, MPFR_RNDN);
  mpfr_sub(den2, den2, aux, MPFR_RNDN);

  // num2 = (M1+M2) (G*sin_th1*cos_del - L1*w1_sqr*sin_del - G*sin_th2) - M2*L2*w2_sqr*sin_cos_del
  mpfr_mul_d(num2, sin_th1, G, MPFR_RNDN);
  mpfr_mul(num2, num2, cos_del, MPFR_RNDN);   // num2 = G*sin_th1*cos_del
  mpfr_mul_d(aux, w1_sqr, L1, MPFR_RNDN);
  mpfr_mul(aux, aux, sin_del, MPFR_RNDN);     // aux = L1*w1_sqr*sin_del 
  mpfr_sub(num2, num2, aux, MPFR_RNDN);       // num2 -= aux;
  mpfr_mul_d(aux, sin_th2, G, MPFR_RNDN);     // aux = G*sin_th2
  mpfr_sub(num2, num2, aux, MPFR_RNDN);       // num2 -= aux
  mpfr_mul(num2, num2, mass_sum, MPFR_RNDN);  // num2 *= (M1+M2)
  mpfr_set_d(aux, M2 * L2, MPFR_RNDN);
  mpfr_mul(aux, aux, w2_sqr, MPFR_RNDN);
  mpfr_mul(aux, aux, sin_cos_del, MPFR_RNDN); // aux = M2*L2*w2_sqr*sin_cos_del
  mpfr_sub(num2, num2, aux, MPFR_RNDN);       // num2 *= aux

  mpfr_div(dydx->w2, num2, den2, MPFR_RNDN);  // w2' = num2/den2 


  /*  Clean up. */
  mpfr_clears(num1, num2, den1, den2, aux, mass_sum, del, cos_del, sin_del, 
              sin_cos_del, sin_th1, sin_th2, w1_sqr, w2_sqr, NULL);
  mpfr_free_cache();
}

/* 
 * Preforms runge-kutta method to integrate yin at time t.
 * Result is stored in yout. 
 */
void runge_kutta(mpfr_t t, y_t *yin, y_t *yout, mpfr_t h) 
{
 
  y_t dydx, dydxt, yt, k1, k2, k3, k4;
  mpfr_t aux;
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
  mpfr_init2(aux, nbits);
  
  /* First step. */
  derivs(yin, &dydx);

  mpfr_mul(k1.th1, h, dydx.th1, MPFR_RNDN);       // k1.th1 = h*dydx.th1
  mpfr_mul_d(yt.th1, k1.th1, 0.5, MPFR_RNDN);  // yt.th1 = yin->th1 + 0.5*k1.th1
  mpfr_add(yt.th1, yt.th1, yin->th1, MPFR_RNDN);

  mpfr_mul(k1.w1, h, dydx.w1, MPFR_RNDN);         // k1.w1 = h*dydx.w1
  mpfr_mul_d(yt.w1, k1.w1, 0.5, MPFR_RNDN);    // yt.w1 = yin->w1 + 0.5*k1.w1
  mpfr_add(yt.w1, yt.w1, yin->w1, MPFR_RNDN);

  mpfr_mul(k1.th2, h, dydx.th2, MPFR_RNDN);       // k1.th2 = h*dydx.th2
  mpfr_mul_d(yt.th2, k1.th2, 0.5, MPFR_RNDN);  // yt.th2 = yin->th2 + 0.5*k1.th2
  mpfr_add(yt.th2, yt.th2, yin->th2, MPFR_RNDN);

  mpfr_mul(k1.w2, h, dydx.w2, MPFR_RNDN);         // k1.w2 = h*dydx.w2
  mpfr_mul_d(yt.w2, k1.w2, 0.5, MPFR_RNDN);    // yt.w2 = yin->w2 + 0.5*k1.w2
  mpfr_add(yt.w2, yt.w2, yin->w2, MPFR_RNDN);
  

  /* Second step. */
  derivs(&yt, &dydxt);

  mpfr_mul(k2.th1, h, dydxt.th1, MPFR_RNDN);      // k2.th1 = h*dydxt.th1
  mpfr_mul_d(yt.th1, k2.th1, 0.5, MPFR_RNDN);  // yt.th1 = yin->th1 + 0.5*k2.th1
  mpfr_add(yt.th1, yt.th1, yin->th1, MPFR_RNDN);

  mpfr_mul(k2.w1, h, dydxt.w1, MPFR_RNDN);        // k2.w1 = h*dydxt.w1
  mpfr_mul_d(yt.w1, k2.w1, 0.5, MPFR_RNDN);    // yt.w1 = yin->w1 + 0.5*k2.w1
  mpfr_add(yt.w1, yt.w1, yin->w1, MPFR_RNDN);

  mpfr_mul(k2.th2, h, dydxt.th2, MPFR_RNDN);      // k2.th2 = h*dydxt.th2
  mpfr_mul_d(yt.th2, k2.th2, 0.5, MPFR_RNDN);  // yt.th2 = yin->th2 + 0.5*k2.th2
  mpfr_add(yt.th2, yt.th2, yin->th2, MPFR_RNDN);

  mpfr_mul(k2.w2, h, dydxt.w2, MPFR_RNDN);        // k2.w2 = h*dydxt.w2
  mpfr_mul_d(yt.w2, k2.w2, 0.5, MPFR_RNDN);    // yt.w2 = yin->w2 + 0.5*k2.w2  
  mpfr_add(yt.w2, yt.w2, yin->w2, MPFR_RNDN);
  

  /* Third step. */
  derivs(&yt, &dydxt);

  mpfr_mul(k3.th1, h, dydxt.th1, MPFR_RNDN);      // k3.th1 = h*dydxt.th1;
  mpfr_add(yt.th1, k3.th1, yin->th1, MPFR_RNDN);  // yt.th1 = yin->th1 + k3.th1;

  mpfr_mul(k3.w1, h, dydxt.w1, MPFR_RNDN);        // k3.w1 = h*dydxt.w1;
  mpfr_add(yt.w1, k3.w1, yin->w1, MPFR_RNDN);     // yt.w1 = yin->w1 + k3.w1;

  mpfr_mul(k3.th2, h, dydxt.th2, MPFR_RNDN);      // k3.th2 = h*dydxt.th2;
  mpfr_add(yt.th2, k3.th2, yin->th2, MPFR_RNDN);  // yt.th2 = yin->th2 + k3.th2;

  mpfr_mul(k3.w2, h, dydxt.w2, MPFR_RNDN);        // k3.w2 = h*dydxt.w2;
  mpfr_add(yt.w2, k3.w2, yin->w2, MPFR_RNDN);     // yt.w2 = yin->w2 + k3.w2;


  /* Fourth step. */
  derivs(&yt, &dydxt);

  mpfr_mul(k4.th1, h, dydxt.th1, MPFR_RNDN);      // k4.th1 = h*dydxt.th1;
  // yout->th1 = yin->th1 + (k1.th1 + k4.th1)/6.0 + (k2.th1 + k3.th1)/3.0
  mpfr_set(yout->th1, yin->th1, MPFR_RNDN);
  mpfr_add(aux, k1.th1, k4.th1, MPFR_RNDN);
  mpfr_div_d(aux, aux, 6.0, MPFR_RNDN);
  mpfr_add(yout->th1, yout->th1, aux, MPFR_RNDN);
  mpfr_add(aux, k2.th1, k3.th1, MPFR_RNDN);
  mpfr_div_d(aux, aux, 3.0, MPFR_RNDN);
  mpfr_add(yout->th1, yout->th1, aux, MPFR_RNDN);


  mpfr_mul(k4.w1, h, dydxt.w1, MPFR_RNDN);        // k4.w1 = h*dydxt.w1;
  // yout->w1 = yin->w1 + (k1.w1 + k4.w1)/6.0 + (k2.w1 + k3.w1)/3.0
  mpfr_set(yout->w1, yin->w1, MPFR_RNDN);
  mpfr_add(aux, k1.w1, k4.w1, MPFR_RNDN);
  mpfr_div_d(aux, aux, 6.0, MPFR_RNDN);
  mpfr_add(yout->w1, yout->w1, aux, MPFR_RNDN);
  mpfr_add(aux, k2.w1, k3.w1, MPFR_RNDN);
  mpfr_div_d(aux, aux, 3.0, MPFR_RNDN);
  mpfr_add(yout->w1, yout->w1, aux, MPFR_RNDN);


  mpfr_mul(k4.th2, h, dydxt.th2, MPFR_RNDN);      // k4.th2 = h*dydxt.th2;
  // yout->th2 = yin->th2 + (k1.th2 + k4.th2)/6.0 + (k2.th2 + k3.th2)/3.0
  mpfr_set(yout->th2, yin->th2, MPFR_RNDN);
  mpfr_add(aux, k1.th2, k4.th2, MPFR_RNDN);
  mpfr_div_d(aux, aux, 6.0, MPFR_RNDN);
  mpfr_add(yout->th2, yout->th2, aux, MPFR_RNDN);
  mpfr_add(aux, k2.th2, k3.th2, MPFR_RNDN);
  mpfr_div_d(aux, aux, 3.0, MPFR_RNDN);
  mpfr_add(yout->th2, yout->th2, aux, MPFR_RNDN);


  mpfr_mul(k4.w2, h, dydxt.w2, MPFR_RNDN);        // k4.w2 = h*dydxt.w2;
  // yout->w2 = yin->w2 + (k1.w2 + k4.w2)/6.0 + (k2.w2 + k3.w2)/3.0
  mpfr_set(yout->w2, yin->w2, MPFR_RNDN);
  mpfr_add(aux, k1.w2, k4.w2, MPFR_RNDN);
  mpfr_div_d(aux, aux, 6.0, MPFR_RNDN);
  mpfr_add(yout->w2, yout->w2, aux, MPFR_RNDN);
  mpfr_add(aux, k2.w2, k3.w2, MPFR_RNDN);
  mpfr_div_d(aux, aux, 3.0, MPFR_RNDN);
  mpfr_add(yout->w2, yout->w2, aux, MPFR_RNDN);


 /* Clean up. */
  mpfr_clears(dydx.th1, dydx.w1, dydx.th2, dydx.w2, 
              dydxt.th1, dydxt.w1, dydxt.th2, dydxt.w2, 
              yt.th1, yt.w1, yt.th2, yt.w2, 
              k1.th1, k1.w1, k1.th2, k1.w2,
              k2.th1, k2.w1, k2.th2, k2.w2,
              k3.th1, k3.w1, k3.th2, k3.w2, 
              k4.th1, k4.w1, k4.th2, k4.w2, 
              aux, THREE, SIX, ONE_HALF, NULL);
  mpfr_free_cache();
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

void reset_yin_s(mpfr_t *d0, mpfr_t *di, y_t *y0, y_t *y1_out, y_t *y1_in) {
  mpfr_t d;
  mpfr_init_set_d(d, 0.0, 53);

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







