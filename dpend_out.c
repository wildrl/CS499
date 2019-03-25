/* dpend_out.c
 *
 * Used by dpend_mpfr.c to write out data.
 *
 */

#include "dpend.h"
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <mpfr.h>
#include <dirent.h>
#include <sys/stat.h>

void output_initial_conditions(char* argv[]) {
  all_ics = fopen("./mpfr_data/all_ics.txt", "a");
  fprintf(all_ics, "ic_%d,%s,%s,%s,%s,%s,%f\n", 
          dir_count, argv[1], argv[2], argv[3], argv[4], argv[5], h);
  fclose(all_ics);
}

void create_output_directory() {
  /* Count number of files in mpfr_data in order to name the new file to be created. */
  dir_count = 0;
  DIR *dirp;
  struct dirent *entry;
    
  dirp = opendir("./mpfr_data");
  while ((entry = readdir(dirp)) != NULL) {
    if (entry->d_type == DT_DIR) { /* If the entry is a regular file */
      dir_count++;
    }
  }
  closedir(dirp);

  /* Create ic_# directory to store output in. */
  snprintf(dir_name, 30, "./mpfr_data/ic_%d", dir_count);
  mkdir(dir_name, 0777);
}

void create_out_files(int num) {
    char polar_fn[50], cartesian_fn[50], mag_fn[50]; //energy_fn[50];

    /* Create output file for polar coordinates. */
    snprintf(polar_fn, 50, "%s/polar%d.csv", dir_name, num);
    polar_output = fopen(polar_fn, "w");
    fprintf(polar_output, "time,th1,w1,th2,w2\n");

    /* Create output file for cartesian coordinates. */
    snprintf(cartesian_fn, 50, "%s/cartesian%d.txt", dir_name, num);
    cartesian_output = fopen(cartesian_fn, "w");
    fprintf(cartesian_output, "time,x1,y1,x2,y2\n");

    /* Create output file for magnitude. */
    snprintf(mag_fn, 50, "%s/mag%d.csv", dir_name, num);
    mag_output = fopen(mag_fn, "w");
    fprintf(mag_output, "time,magnitude,dot_product\n");

    /* Create output file for energy values. */
    //snprintf(energy_fn, 50, "%s/energy%d.csv", dir_name, num);
    //energy_output = fopen(energy_fn, "w");
    //fprintf(energy_output, "time,KE,PE,Total\n");
}

void output_mag(mpfr_t* t, mpfr_t *mag, mpfr_t *dot, mpfr_t *r_error) {

  if (nbits != 113) {
    mpfr_fprintf(mag_output, "%0.32RNF,%0.32RNF,%0.32RNF,%0.32RNF\n", 
              *t, mag, dot, r_error);
  } else {
    mpfr_fprintf(mag_output, "%0.32RNF,%0.32RNF,0\n", 
              *t, mag);
  }

}

void output_polar(mpfr_t* t, y_t* y) 
{
  mpfr_fprintf(polar_output, "%0.32RNf,%0.32RNf,%0.32RNF,%0.32RNF,%0.32RNF\n", 
              *t, y->th1, y->w1, y->th2, y->w2);
}

void output_cartesian(mpfr_t* t, y_t* y) 
{

  /* Convert to cartesian coordinates. */
  mpfr_t x1, y1, x2, y2, temp;
  mpfr_inits2(nbits, x1, y1, x2, y2, temp, NULL);

  // Calculate x1
  mpfr_set_d(x1, L1, MPFR_RNDN);
  mpfr_sin(temp, y->th1, MPFR_RNDN);
  mpfr_mul(x1, x1, temp, MPFR_RNDN);

  // Calculate y1
  mpfr_set_d(y1, -L1, MPFR_RNDN);
  mpfr_cos(temp, y->th1, MPFR_RNDN);
  mpfr_mul(y1, y1, temp, MPFR_RNDN);

  // Calculate x2
  mpfr_set_d(x2, L2, MPFR_RNDN);
  mpfr_sin(temp, y->th2, MPFR_RNDN);
  mpfr_mul(x2, x2, temp, MPFR_RNDN);
  mpfr_add(x2, x2, x1, MPFR_RNDN);

  // Calculate y2
  mpfr_set_d(y2, -L2, MPFR_RNDN);
  mpfr_cos(temp, y->th2, MPFR_RNDN);
  mpfr_mul(y2, y2, temp, MPFR_RNDN);
  mpfr_add(y2, y2, y1, MPFR_RNDN);

  mpfr_fprintf(cartesian_output, "%0.32RNf,%0.32RNf,%0.32RNF,%0.32RNF,%0.32RNF\n", 
               *t, x1, y1, x2, y2); 
}

void output_energy(mpfr_t* t, y_t* y){
  
  mpfr_t p_energy, k_energy, t_energy, temp, temp2, cos_th1, cos_th2, del, cos_del;
  mpfr_inits2(nbits, p_energy, k_energy, t_energy, temp, temp2, cos_th1, cos_th2, del, cos_del, NULL);
  double half = 0.5;
  double neg = -1.0;
  double mass_sum = M1+M2;

  // PE
  mpfr_cos(cos_th1, y->th1, MPFR_RNDN);     
  mpfr_cos(cos_th2, y->th2, MPFR_RNDN);

  mpfr_mul_d(p_energy, cos_th1, L1, MPFR_RNDN); 
  mpfr_mul_d(p_energy, p_energy, mass_sum, MPFR_RNDN);
  mpfr_mul_d(temp, cos_th2, L2, MPFR_RNDN);
  mpfr_mul_d(temp, temp, M2, MPFR_RNDN);
  mpfr_add(p_energy, p_energy, temp, MPFR_RNDN);
  mpfr_mul_d(p_energy, p_energy, -G, MPFR_RNDN);

  //KE
  mpfr_mul(temp, y->w1, y->w1, MPFR_RNDN);
  mpfr_mul_d(temp, temp, L1, MPFR_RNDN);
  mpfr_mul_d(temp, temp, L1, MPFR_RNDN);
  mpfr_mul_d(temp, temp, M1, MPFR_RNDN);

  mpfr_mul(temp2, y->w2, y->w2, MPFR_RNDN);
  mpfr_mul_d(temp2, temp2, L2, MPFR_RNDN);
  mpfr_mul_d(temp2, temp2, L2, MPFR_RNDN);
  mpfr_mul_d(temp2, temp2, M2, MPFR_RNDN);

  mpfr_add(temp, temp, temp2, MPFR_RNDN);

  mpfr_mul(temp2, y->w1, y->w1, MPFR_RNDN);
  mpfr_mul_d(temp2, temp2, L1, MPFR_RNDN);
  mpfr_mul_d(temp2, temp2, L1, MPFR_RNDN);
  mpfr_mul_d(temp2, temp2, M2, MPFR_RNDN);

  mpfr_add(temp, temp, temp2, MPFR_RNDN);
  mpfr_mul_d(temp, temp, half, MPFR_RNDN);
  
  mpfr_mul(temp2, y->w1, y->w2, MPFR_RNDN);
  mpfr_mul_d(temp2, temp2, L1, MPFR_RNDN);
  mpfr_mul_d(temp2, temp2, L2, MPFR_RNDN);
  mpfr_mul_d(temp2, temp2, M2, MPFR_RNDN);
  mpfr_sub(del, y->th1, y->th2, MPFR_RNDN);
  mpfr_cos(cos_del, del, MPFR_RNDN); 
  mpfr_mul(temp2, temp2, cos_del, MPFR_RNDN);

  mpfr_add(k_energy, temp, temp2, MPFR_RNDN);
  mpfr_mul_d(temp2, p_energy, neg, MPFR_RNDN);
  mpfr_mul(k_energy, temp, temp2, MPFR_RNDN);

  //total energy
  mpfr_add(t_energy, p_energy, k_energy, MPFR_RNDN);

  mpfr_fprintf(energy_output, "%0.32RNf,%0.32RNf,%0.32RNF,%0.32RNF\n", 
               *t, p_energy, k_energy, t_energy); 

  
}


