#include "../include/cec2015_interface.h"

/**
 * An interface for CEC 2017 functions.
 *
 * @param x
 * @param f
 * @param nx
 * @param mx
 * @param func_num
 */

void cec2015_func(char *datapath, double *x, double *f, int nx, int mx,
                  int func_num) {
  if (!(nx == 2 || nx == 10 || nx == 30 || nx == 50 || nx == 100)) {
    perror("Error: Test functions are only defined for D = 2, 10, 20, 30, 50, "
           "100.");
  }
  if (nx == 2 && ((func_num >= 6 && func_num <= 8) || (func_num == 10) ||
                  (func_num == 13))) {
    perror("Error: hf0{1..3}, cf0{2..5} are NOT defined for D=2.");
  }

  int cf_nums[] = {0, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 5, 5, 5, 7, 10};
  int bShuffle[] = {0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0};
  int cf_num = cf_nums[func_num];
  char FileName[256];
  /*
   * Load M matrix:
   */
  sprintf(FileName, "%s/M_%d_D%d.txt", datapath, func_num, nx);
  FILE *fptMData = fopen(FileName, "r");
  if (fptMData == NULL) {
    perror("Error: Cannot open input file for reading");
  }
  int MatrixSize = cf_num * nx * nx;
  double *M = malloc(MatrixSize * sizeof(double));
  if (M == NULL) {
    perror("Error: there is insufficient memory available!");
  } else {
    for (int i = 0; i < MatrixSize; ++i) {
      if (fscanf(fptMData, "%lf", &M[i]) == -1) {
        perror("Cannot read M matrix data.");
      }
    }
  }
  fclose(fptMData);
  /*
   * Load bias data:
   */
  double *bias = NULL;
  if (cf_num > 1) {
    sprintf(FileName, "%s/bias_%d.txt", datapath, func_num);
    FILE *fptBiasData = fopen(FileName, "r");
    if (fptBiasData == NULL) {
      perror("Error: Cannot open input file for reading");
    }
    bias = malloc(cf_num * sizeof(double));
    if (M == NULL) {
      perror("Error: there is insufficient memory available!");
    } else {
      for (int i = 0; i < cf_num; ++i) {
        if (fscanf(fptBiasData, "%lf", &bias[i]) == -1) {
          perror("Cannot read M matrix data.");
        }
      }
    }
    fclose(fptBiasData);
  }
  /*
   * Load shift data:
   */
  sprintf(FileName, "%s/shift_data_%d.txt", datapath, func_num);
  FILE *fptOShiftData = fopen(FileName, "r");
  if (fptOShiftData == NULL) {
    perror("Error: Cannot open input file for reading");
  }
  double *OShift = (double *)malloc(cf_num * nx * sizeof(double));
  if (OShift == NULL) {
    perror("Error: there is insufficient memory available!");
  }
  char tmpchar;
  for (int i = 0; i < nx * cf_nums[func_num]; ++i) {
    if (fscanf(fptOShiftData, "%lf", &OShift[i]) == -1) {
      perror("Cannot read shift data.");
    }
    if (cf_nums[func_num] > 1 && ((i + 1) % nx) == 0) {
      fscanf(fptOShiftData, "%c", &tmpchar);
      while (tmpchar != '\n') {
        fscanf(fptOShiftData, "%c", &tmpchar);
      }
    }
  }
  fclose(fptOShiftData);

  /*
   * Load shuffle data:
   */
  int *SS = NULL;
  int shuffleFlag = bShuffle[func_num] == 1 ? 1 : 0;
  if (shuffleFlag) {
    sprintf(FileName, "%s/shuffle_data_%d_D%d.txt", datapath, func_num, nx);
    FILE *fptShuffleData = fopen(FileName, "r");
    if (fptShuffleData == NULL) {
      perror("Error: Cannot open input file for reading");
    }
    int ShuffleSize = nx * cf_num;
    SS = malloc(ShuffleSize * sizeof(int));
    for (int i = 0; i < ShuffleSize; ++i) {
      if (fscanf(fptShuffleData, "%d", &SS[i]) == -1) {
        perror("Cannot read shuffle data");
      }
    }
    fclose(fptShuffleData);
  }
  for (int i = 0; i < mx; i++) {
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
      perror("Error: There are only 15 test functions in this test suite! [CEC2015-LB]");
      f[i] = 0.0;
      break;
    }
  }
  free(M);
  free(OShift);
  if (shuffleFlag) {
    free(SS);
  }
  if (cf_num > 1) {
    free(bias);
  }
}
