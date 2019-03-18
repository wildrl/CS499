#include <stdio.h>
#include <stdlib.h>
//#include "dpend.h"
#include <math.h>
#include <gmp.h>
#include <mpfr.h>

void runge_kutta(mpfr_t t, y_t *yin, y_t *yout, mpfr_t h);
void derivs(y_t *yin, y_t *dydx);

void magnitude (y_t *y, mpfr_t *magnitude);
void dot_product(y_t *y1, y_t *y2, mpfr_t *dot);
