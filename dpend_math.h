#include <stdio.h>
#include <stdlib.h>
//#include "dpend.h"
#include <math.h>
#include <gmp.h>
#include <mpfr.h>

void runge_kutta(mpfr_t t, y_t *yin, y_t *yout, mpfr_t h);
void derivs(y_t *yin, y_t *dydx);

void reset_yin_s(mpfr_t *d0, mpfr_t *di, y_t *y0, y_t *y1_out, y_t *y1_in);
void lyapunov(mpfr_t *sum, mpfr_t *d0, mpfr_t *d1);
void calc_di(mpfr_t *di, y_t *y0, y_t *y1);
