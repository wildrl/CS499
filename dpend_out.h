#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <mpfr.h>

void output_polar(FILE *file, mpfr_t* t, mpfr_t* th1, mpfr_t* w1, mpfr_t* th2, mpfr_t* w2);
void output_cartesian(FILE *file, int nbits, mpfr_t* t, mpfr_t* th1, mpfr_t* w1, mpfr_t* th2, mpfr_t* w2, double L1, double L2);
output_energy(FILE *file, int nbits, mpfr_t* t, mpfr_t* th1, mpfr_t* w1, mpfr_t* th2, mpfr_t* w2, double L1, double L2, double G);
