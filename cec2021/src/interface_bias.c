#include "interface_bias.h"

void cec21_bias_func(double *x, double *f, int nx, int mx, int func_num) {

  int cf_num = 10, i, j;

  if (func_num < 1 || func_num > 10) {
    printf("\nError: Test function %d is not defined.\n", func_num);
  }

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

    if (!(nx == 2 || nx == 10 || nx == 20)) {
      printf("\nError: Test functions are only defined for D=2,10,20.\n");
    }
    if (nx == 2 && (func_num == 5 || func_num == 6 || func_num == 7)) {
      printf("\nError:  NOT defined for D=2.\n");
    }

    /* Load Matrix M*/
    sprintf(FileName, "%s/M_%d_D%d_nr.txt", extdata, func_num, nx);
    fpt = fopen(FileName, "r");
    if (fpt == NULL) {
      printf("\n Error: Cannot open M_%d_D%d_nr.txt for reading \n", func_num,
             nx);
    }
    if (func_num < 7) {
      M = (double *)malloc(nx * nx * sizeof(double));
      if (M == NULL)
        printf("\nError: there is insufficient memory available!\n");
      for (i = 0; i < nx * nx; i++) {
        fscanf(fpt, "%lf", &M[i]);
      }
    } else {
      M = (double *)malloc(cf_num * nx * nx * sizeof(double));
      if (M == NULL)
        printf("\nError: there is insufficient memory available!\n");
      for (i = 0; i < cf_num * nx * nx; i++) {
        fscanf(fpt, "%lf", &M[i]);
      }
    }
    fclose(fpt);

    /* Load shift_data */
    sprintf(FileName, "%s/shift_data_%d_ns.txt", extdata, func_num);
    fpt = fopen(FileName, "r");
    if (fpt == NULL) {
      printf("\n Error: Cannot open shift_data_%d_ns.txt for reading \n",
             func_num);
    }

    if (func_num < 7) {
      OShift = (double *)malloc(nx * sizeof(double));
      if (OShift == NULL)
        printf("\nError: there is insufficient memory available!\n");
      for (i = 0; i < nx; i++) {
        fscanf(fpt, "%lf", &OShift[i]);
      }
    } else {
      // OShift=(double *)malloc(nx*sizeof(double));
      OShift = (double *)malloc(nx * cf_num * sizeof(double));
      if (OShift == NULL)
        printf("\nError: there is insufficient memory available!\n");
      for (i = 0; i < cf_num - 1; i++) {
        for (j = 0; j < nx; j++) {
          fscanf(fpt, "%lf", &OShift[i * nx + j]);
        }
        fscanf(fpt, "%*[^\n]%*c");
      }
      for (j = 0; j < nx; j++) {
        fscanf(fpt, "%lf", &OShift[nx * (cf_num - 1) + j]);
      }
    }
    fclose(fpt);

    /* Load Shuffle_data */

    if (func_num >= 5 && func_num <= 7) {
      sprintf(FileName, "%s/shuffle_data_%d_D%d.txt", extdata, func_num, nx);
      fpt = fopen(FileName, "r");
      if (fpt == NULL) {
        printf("\n Error: Cannot open shuffle_data_%d_D%d.txt for reading \n",
               func_num, nx);
      }
      SS = (int *)malloc(nx * sizeof(int));
      if (SS == NULL)
        printf("\nError: there is insufficient memory available!\n");
      for (i = 0; i < nx; i++) {
        fscanf(fpt, "%d", &SS[i]);
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
      bent_cigar_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      f[i] += 100.0;
      break;
    case 2:
      schwefel_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      f[i] += 1100.0;
      break;
    case 3:
      bi_rastrigin_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      f[i] += 700.0;
      break;
    case 4:
      grie_rosen_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      f[i] += 1900.0;
      break;
    case 5:
      hf01(&x[i * nx], &f[i], nx, OShift, M, SS, 1, 1);
      f[i] = f[i] + 1700.0;
      break;
    case 6:
      hf06(&x[i * nx], &f[i], nx, OShift, M, SS, 1, 1);
      f[i] += 1600.0;
      break;
    case 7:
      hf05(&x[i * nx], &f[i], nx, OShift, M, SS, 1, 1);
      f[i] += 2100.0;
      break;
    case 8:
      cf02(&x[i * nx], &f[i], nx, OShift, M, 1);
      f[i] += 2200.0;
      break;
    case 9:
      cf04(&x[i * nx], &f[i], nx, OShift, M, 1);
      f[i] += 2400.0;
      break;
    case 10:
      cf05(&x[i * nx], &f[i], nx, OShift, M, 1);
      f[i] += 2500.0;
      break;
    default:
      printf("\nError: There are only 10 test functions in this test suite!\n");
      f[i] = 0.0;
      break;
    }
  }
}