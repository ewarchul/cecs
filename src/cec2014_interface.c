#include "cec2014_interface.h"
#include <string.h>

void cec2014_func(double *x, double *f, int nx, int mx, int func_num) {
  int cf_num = 10, i, j;
  if (ini_flag == 1) {
    if ((n_flag != nx) || (func_flag != func_num)) {
      ini_flag = 0;
    }
  }

  if (ini_flag == 0) {
    FILE *fpt;
    char FileName[256];
    free(M);
    free(OShift);
    free(y);
    free(z);
    free(x_bound);
    y = (double *)malloc(sizeof(double) * nx);
    z = (double *)malloc(sizeof(double) * nx);
    x_bound = (double *)malloc(sizeof(double) * nx);
    for (i = 0; i < nx; i++)
      x_bound[i] = 100.0;

    if (!(nx == 2 || nx == 10 || nx == 20 || nx == 30 || nx == 50 ||
          nx == 100)) {
      printf("\nError: Test functions are only defined for "
             "D=2,10,20,30,50,100.\n");
    }
    if (nx == 2 && ((func_num >= 17 && func_num <= 22) ||
                    (func_num >= 29 && func_num <= 30))) {
      printf("\nError: hf01,hf02,hf03,hf04,hf05,hf06,cf07&cf08 are NOT defined "
             "for D=2.\n");
    }

    /* Load Matrix M*/
    sprintf(FileName, "%s/M_%d_D%d.txt", extdata, func_num, nx);
    fpt = fopen(FileName, "r");
    if (fpt == NULL) {
      printf("\n Error: Cannot open input file for reading \n");
    }
    if (func_num < 23) {
      M = (double *)malloc(nx * nx * sizeof(double));
      if (M == NULL)
        printf("\nError: there is insufficient memory available!\n");
      for (i = 0; i < nx * nx; i++) {
        if (fscanf(fpt, "%lf", &M[i]) != 1) {
          printf("\nError\n");
        }
      }
    } else {
      M = (double *)malloc(cf_num * nx * nx * sizeof(double));
      if (M == NULL)
        printf("\nError: there is insufficient memory available!\n");
      for (i = 0; i < cf_num * nx * nx; i++) {
        if (fscanf(fpt, "%lf", &M[i]) != 1) { 
          printf("\nError\n");
        }
      }
    }
    fclose(fpt);

    /* Load shift_data */

    sprintf(FileName, "%s/shift_data_%d.txt", extdata, func_num);
    fpt = fopen(FileName, "r");
    if (fpt == NULL) {
      printf("\n Error: Cannot open input file for reading \n");
    }

    if (func_num < 23) {
      OShift = (double *)malloc(nx * sizeof(double));
      if (OShift == NULL)
        printf("\nError: there is insufficient memory available!\n");
      for (i = 0; i < nx; i++) {
        if (fscanf(fpt, "%lf", &OShift[i]) != 1) {
          printf("\nError\n");
        }
      }
    } else {
      OShift = (double *)malloc(nx * cf_num * sizeof(double));
      if (OShift == NULL)
        printf("\nError: there is insufficient memory available!\n");
      for (i = 0; i < cf_num - 1; i++) {
        for (j = 0; j < nx; j++) {
          if (fscanf(fpt, "%lf", &OShift[i * nx + j]) != 1) {
          printf("\nError\n");
          }
        }
        if (fscanf(fpt, "%*[^\n]%*c") != 1) {
          printf("\nError\n");
        }
      }
      for (j = 0; j < nx; j++) {
        if (fscanf(fpt, "%lf", &OShift[(cf_num - 1) * nx + j]) != 1) {
          printf("\nError\n");
        }
      }
    }
    fclose(fpt);

    /* Load Shuffle_data */

    if (func_num >= 17 && func_num <= 22) {

      sprintf(FileName, "%s/shuffle_data_%d_D%d.txt", extdata, func_num, nx);
      fpt = fopen(FileName, "r");
      if (fpt == NULL) {
        printf("\n Error: Cannot open input file for reading \n");
      }
      SS = (int *)malloc(nx * sizeof(int));
      if (SS == NULL)
        printf("\nError: there is insufficient memory available!\n");
      for (i = 0; i < nx; i++) {
        if (fscanf(fpt, "%d", &SS[i]) != 1) {
          printf("\nError\n");
        }
      }
      fclose(fpt);
    } else if (func_num == 29 || func_num == 30) {

      sprintf(FileName, "%s/shuffle_data_%d_D%d.txt", extdata, func_num, nx);
      fpt = fopen(FileName, "r");
      if (fpt == NULL) {
        printf("\n Error: Cannot open input file for reading \n");
      }
      SS = (int *)malloc(nx * cf_num * sizeof(int));
      if (SS == NULL)
        printf("\nError: there is insufficient memory available!\n");
      for (i = 0; i < nx * cf_num; i++) {
        if (fscanf(fpt, "%d", &SS[i]) != 1) {
          printf("\nError\n");
        }
      }
      fclose(fpt);
    }

    n_flag = nx;
    func_flag = func_num;
    ini_flag = 1;
  }

  for (i = 0; i < mx; i++) {
    switch (func_num) {
    case 1:
      cec2014_ellips_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      f[i] += 100.0;
      break;
    case 2:
      cec2014_bent_cigar_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      f[i] += 200.0;
      break;
    case 3:
      cec2014_discus_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      f[i] += 300.0;
      break;
    case 4:
      cec2014_rosenbrock_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      f[i] += 400.0;
      break;
    case 5:
      cec2014_ackley_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      f[i] += 500.0;
      break;
    case 6:
      cec2014_weierstrass_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      f[i] += 600.0;
      break;
    case 7:
      cec2014_griewank_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      f[i] += 700.0;
      break;
    case 8:
      cec2014_rastrigin_func(&x[i * nx], &f[i], nx, OShift, M, 1, 0);
      f[i] += 800.0;
      break;
    case 9:
      cec2014_rastrigin_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      f[i] += 900.0;
      break;
    case 10:
      cec2014_schwefel_func(&x[i * nx], &f[i], nx, OShift, M, 1, 0);
      f[i] += 1000.0;
      break;
    case 11:
      cec2014_schwefel_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      f[i] += 1100.0;
      break;
    case 12:
      cec2014_katsuura_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      f[i] += 1200.0;
      break;
    case 13:
      cec2014_happycat_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      f[i] += 1300.0;
      break;
    case 14:
      cec2014_hgbat_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      f[i] += 1400.0;
      break;
    case 15:
      cec2014_grie_rosen_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      f[i] += 1500.0;
      break;
    case 16:
      cec2014_escaffer6_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      f[i] += 1600.0;
      break;
    case 17:
      cec2014_hf01(&x[i * nx], &f[i], nx, OShift, M, SS, 1, 1);
      f[i] += 1700.0;
      break;
    case 18:
      cec2014_hf02(&x[i * nx], &f[i], nx, OShift, M, SS, 1, 1);
      f[i] += 1800.0;
      break;
    case 19:
      cec2014_hf03(&x[i * nx], &f[i], nx, OShift, M, SS, 1, 1);
      f[i] += 1900.0;
      break;
    case 20:
      cec2014_hf04(&x[i * nx], &f[i], nx, OShift, M, SS, 1, 1);
      f[i] += 2000.0;
      break;
    case 21:
      cec2014_hf05(&x[i * nx], &f[i], nx, OShift, M, SS, 1, 1);
      f[i] += 2100.0;
      break;
    case 22:
      cec2014_hf06(&x[i * nx], &f[i], nx, OShift, M, SS, 1, 1);
      f[i] += 2200.0;
      break;
    case 23:
      cec2014_cf01(&x[i * nx], &f[i], nx, OShift, M, 1);
      f[i] += 2300.0;
      break;
    case 24:
      cec2014_cf02(&x[i * nx], &f[i], nx, OShift, M, 1);
      f[i] += 2400.0;
      break;
    case 25:
      cec2014_cf03(&x[i * nx], &f[i], nx, OShift, M, 1);
      f[i] += 2500.0;
      break;
    case 26:
      cec2014_cf04(&x[i * nx], &f[i], nx, OShift, M, 1);
      f[i] += 2600.0;
      break;
    case 27:
      cec2014_cf05(&x[i * nx], &f[i], nx, OShift, M, 1);
      f[i] += 2700.0;
      break;
    case 28:
      cec2014_cf06(&x[i * nx], &f[i], nx, OShift, M, 1);
      f[i] += 2800.0;
      break;
    case 29:
      cec2014_cf07(&x[i * nx], &f[i], nx, OShift, M, SS, 1);
      f[i] += 2900.0;
      break;
    case 30:
      cec2014_cf08(&x[i * nx], &f[i], nx, OShift, M, SS, 1);
      f[i] += 3000.0;
      break;
    default:
      printf("\nError: There are only 30 test functions in this test suite!\n");
      f[i] = 0.0;
      break;
    }
  }
}