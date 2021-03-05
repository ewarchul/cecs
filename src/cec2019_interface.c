#include "cec2019/cec2019_interface.h"

void cec2019_func(double *x, double *f, int nx, int mx, int func_num) {
  int cf_num = 10, i, j, ret;
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

    if (!(nx == 2 || nx == 10 || nx == 9 || nx == 16 || nx == 18)) {
      printf("Error: Test functions are only defined for D=10, 9, 16, 18 \n F1 "
             "is defined on D=9 \n F2 is defined on D=16 \n F3 is defined on "
             "D=18 \n F4-F10 are defined on D=10.");
    }

    /* Load Matrix M*/
    if (func_num > 3) {
      sprintf(FileName, "%s/M_%d_D%d.txt", extdata, func_num, nx);
      fpt = fopen(FileName, "r");
      if (fpt == NULL) {
        perror("Error: Cannot open input file for reading");
      }
      M = (double *)malloc(nx * nx * sizeof(double));
      if (M == NULL)
        perror("Error: there is insufficient memory available!");
      for (i = 0; i < nx * nx; i++) {
        ret = fscanf(fpt, "%lf", &M[i]);
        if (ret != 1) {
          perror("Error reading from the input file");
        }
      }
      fclose(fpt);
    }

    /* Load shift_data */
    if (func_num > 3) {
      sprintf(FileName, "%s/shift_data_%d.txt", extdata, func_num);
      fpt = fopen(FileName, "r");
      if (fpt == NULL) {
        perror("Error: Cannot open input file for reading");
      }

      OShift = (double *)malloc(nx * sizeof(double));
      if (OShift == NULL)
        printf("\nError: there is insufficient memory available!\n");
      for (i = 0; i < nx; i++) {
        ret = fscanf(fpt, "%lf", &OShift[i]);
        if (ret != 1) {
          perror("Error reading from the input file");
        }
      }

      fclose(fpt);
    }

    n_flag = nx;
    func_flag = func_num;
    ini_flag = 1;
    // printf("Function has been initialized!\n");
  }

  for (i = 0; i < mx; i++) {
    switch (func_num) {
    case 1:
      cec2019_Chebyshev(&x[i * nx], nx, &f[i]);
      f[i] += 1.0;
      break;

    case 2:
      cec2019_Hilbert(&x[i * nx], nx, &f[i]);
      f[i] += 1.0;
      break;

    case 3:
      cec2019_Lennard_Jones(&x[i * nx], nx, &f[i]);
      f[i] += 1.0;
      break;

    case 4:
      cec2019_rastrigin_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      f[i] += 1.0;
      break;

    case 5:
      cec2019_griewank_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      f[i] += 1.0;
      break;

    case 6:
      cec2019_weierstrass_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      f[i] += 1.0;
      break;

    case 7:
      cec2019_schwefel_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      f[i] += 1.0;
      break;

    case 8:
      cec2019_escaffer6_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      f[i] += 1.0;
      break;

    case 9:
      cec2019_happycat_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      f[i] += 1.0;
      break;

    case 10:
      cec2019_ackley_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      f[i] += 1.0;
      break;
    default:
      printf("\nError: There are only 10 test functions in this test suite!\n");
      f[i] = 0.0;
      break;
    }
  }
}
