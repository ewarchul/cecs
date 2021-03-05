#include "cec2015/cec2015_interface.h"

void cec2015_func(double *x, double *f, int nx, int mx, int func_num) {
  // It's used to skip the rest of line untile get to a newline;
  int cf_nums[] = {0, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 5, 5, 5, 7, 10};
  int bShuffle[] = {0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0};

  char tmpchar;
  int cf_num, i, j;
  int ret;
  double *bias;

  cf_num = cf_nums[func_num];
  if (ini_flag == 1) {
    if ((n_flag != nx) || (func_flag != func_num)) {
      ini_flag = 0;
    }
  }

  if (ini_flag == 0) {
    FILE *fpt;
    char FileName[256];
    free(M);
    free(bias);
    free(OShift);
    free(y);
    free(z);
    free(x_bound);
    y = (double *)malloc(sizeof(double) * nx);
    z = (double *)malloc(sizeof(double) * nx);
    x_bound = (double *)malloc(sizeof(double) * nx);
    for (i = 0; i < nx; i++)
      x_bound[i] = 100.0;

    if (!(nx == 2 || nx == 10 || nx == 30 || nx == 50 || nx == 100)) {
      perror(
          "\nError: Test functions are only defined for D=2,10,30,50,100.\n");
    }
    if (nx == 2 && ((func_num >= 6 && func_num <= 8) || (func_num == 10) ||
                    (func_num == 13))) {
      perror("\nError: hf01,hf02,hf03,cf02&cf05 are NOT defined for D=2.\n");
    }

    /* Load Matrix M*/
    sprintf(FileName, "%s/M_%d_D%d.txt", extdata, func_num, nx);
    fpt = fopen(FileName, "r");
    if (fpt == NULL) {
      perror("\n Error : Cannot open input file for reading \n");
    }

    M = (double *)malloc(cf_num * nx * nx * sizeof(double));
    if (M == NULL) {
      perror("\nError: there is insufficient memory available!\n");
    }

    for (i = 0; i < cf_num * nx * nx; i++) {
      ret = fscanf(fpt, "%lf", &M[i]);
      if (ret != 1) {
        perror("\n Error: Cannot open input file for reading \n");
      }
    }

    fclose(fpt);

    /* Load Bias_value bias*/
    if (cf_num > 1) {
      sprintf(FileName, "%s/bias_%d.txt", extdata, func_num);
      fpt = fopen(FileName, "r");
      if (fpt == NULL) {
        perror("\n Error: Cannot open input file for reading \n");
      }

      bias = (double *)malloc(cf_num * sizeof(double));
      if (bias == NULL) {
        perror("\nError: there is insufficient memory available!\n");
      }
      for (i = 0; i < cf_num; i++) {
        ret = fscanf(fpt, "%lf", &bias[i]);
        if (ret != 1) {
          perror("\n Error: Cannot open input file for reading \n");
        }
      }
      fclose(fpt);
    }

    /* Load shift_data */
    sprintf(FileName, "%s/shift_data_%d.txt", extdata, func_num);
    fpt = fopen(FileName, "r");
    if (fpt == NULL) {
      perror("\n Error: Cannot open input file for reading \n");
    }
    OShift = (double *)malloc(cf_num * nx * sizeof(double));
    if (OShift == NULL) {
      perror("\nError: there is insufficient memory available!\n");
    }

    for (i = 0; i < nx * cf_nums[func_num]; i++) {
      ret = fscanf(fpt, "%lf", &OShift[i]);
      if (ret != 1) {
        perror("\nError: there is insufficient memory available!\n");
      }
      // if finish reading a row, skip the rest of the line to the next line
      // bug fixed 2014/12/26 by Q. Chen
      if (cf_nums[func_num] > 1 && ((i + 1) % nx) == 0) {
        // printf("%d\t%d\t%d\n", cf_nums[func_num], nx, i);
        ret = fscanf(fpt, "%c", &tmpchar);
        if (ret != 1) {
          perror("\nError: there is insufficient memory available!\n");
        }

        while (tmpchar != '\n') {
          // printf("%c", tmpchar);
          ret = fscanf(fpt, "%c", &tmpchar);
          if (ret != 1) {
            perror("\nError: there is insufficient memory available!\n");
          }
        }
      }
    }
    fclose(fpt);

    /* Load Shuffle_data */
    if (bShuffle[func_num] == 1) {
      sprintf(FileName, "%s/shuffle_data_%d_D%d.txt", extdata, func_num, nx);
      fpt = fopen(FileName, "r");
      if (fpt == NULL) {
        perror("\n Error10: Cannot open input file for reading \n");
      }
      SS = (int *)malloc(nx * cf_num * sizeof(int));
      if (SS == NULL) {
        perror("\nError11: there is insufficient memory available!\n");
      }
      for (i = 0; i < nx * cf_num; i++) {
        ret = fscanf(fpt, "%d", &SS[i]);
        if (ret != 1) {
          perror("\n Error12: Cannot open input file for reading \n");
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
      cec2015_ellips_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      f[i] += 100.0;
      break;
    case 2:
      cec2015_bent_cigar_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      f[i] += 200.0;
      break;
    case 3:
      cec2015_ackley_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      f[i] += 300.0;
      break;
    case 4:
      cec2015_rastrigin_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      f[i] += 400.0;
      break;
    case 5:
      cec2015_schwefel_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      f[i] += 500.0;
      break;
    case 6:
      cec2015_hf01(&x[i * nx], &f[i], nx, OShift, M, SS, 1, 1);
      f[i] += 600.0;
      break;
    case 7:
      cec2015_hf02(&x[i * nx], &f[i], nx, OShift, M, SS, 1, 1);
      f[i] += 700.0;
      break;
    case 8:
      cec2015_hf03(&x[i * nx], &f[i], nx, OShift, M, SS, 1, 1);
      f[i] += 800.0;
      break;
    case 9:
      cec2015_cf01(&x[i * nx], &f[i], nx, OShift, M, bias, 1);
      f[i] += 900.0;
      break;
    case 10:
      cec2015_cf02(&x[i * nx], &f[i], nx, OShift, M, SS, bias, 1);
      f[i] += 1000.0;
      break;
    case 11:
      cec2015_cf03(&x[i * nx], &f[i], nx, OShift, M, bias, 1);
      f[i] += 1100.0;
      break;
    case 12:
      cec2015_cf04(&x[i * nx], &f[i], nx, OShift, M, bias, 1);
      f[i] += 1200.0;
      break;
    case 13:
      cec2015_cf05(&x[i * nx], &f[i], nx, OShift, M, SS, bias, 1);
      f[i] += 1300.0;
      break;
    case 14:
      cec2015_cf06(&x[i * nx], &f[i], nx, OShift, M, bias, 1);
      f[i] += 1400.0;
      break;
    case 15:
      cec2015_cf07(&x[i * nx], &f[i], nx, OShift, M, bias, 1);
      f[i] += 1500.0;
      break;
    default:
      perror("\nError: There are only 15 test functions in this test suite!\n");
      f[i] = 0.0;
      break;
    }
  }
}
