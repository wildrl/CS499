//#include "dpend.h"
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <mpfr.h>

/*void output_polar(FILE *file, mpfr_t* t, mpfr_t* th1, mpfr_t* w1, mpfr_t* th2, mpfr_t* w2);
void output_lyapunov(FILE *file, mpfr_t* t, mpfr_t* exp);
void output_cartesian(FILE *file, int nbits, mpfr_t* t, mpfr_t* th1, mpfr_t* w1, mpfr_t* th2, mpfr_t* w2, double L1, double L2);
void output_energy(FILE *file, int nbits, mpfr_t* t, mpfr_t* th1, mpfr_t* w1, mpfr_t* th2, mpfr_t* w2, double L1, double L2, double G);*/

void output_polar(mpfr_t* t, y_t* y);
void output_lyapunov(mpfr_t* t, mpfr_t* exp);
void output_cartesian(mpfr_t* t, y_t* y);
void output_energy(mpfr_t* t, y_t* y);
void output_mag(mpfr_t* t, mpfr_t *mag, mpfr_t *dot);

void create_out_files();
void create_output_directory();
void output_initial_condition();
