#include "../include/cec2021_interface.h"
#include <stdio.h>
#include <string.h>

double getFunctionBias(const int biasFlag, const int fnNumber) {
  double bias = 0.0;
  double fnBiasDict[10] = {100.0,  1100.0, 700.0,  1900.0, 1700.0,
                           1600.0, 2100.0, 2200.0, 2400.0, 2500.0};
  if (biasFlag) {
    bias = fnBiasDict[fnNumber - 1];
  } else {
    bias = 0.0;
  }
  return bias;
}

void cec2021_func(char *datapath, double *x, double *f, int nx, int mx,
                  int func_num, char *suite) {
  if (!(nx == 10 || nx == 20)) {
    perror("Error: Test functions are only defined for D = 10, 20.");
  }
  if (func_num < 1 || func_num > 10) {
    perror("Error: Test function is not defined");
  }

  int cf_num = 10;
  char FileName[256];

  if (!strcmp(suite, "basic") || !strcmp(suite, "bias") ||
      !strcmp(suite, "bias_shift") || !strcmp(suite, "shift")) {
    sprintf(FileName, "%s/M_%d_D%d_nr.txt", datapath, func_num, nx);
  } else {
    sprintf(FileName, "%s/M_%d_D%d.txt", datapath, func_num, nx);
  }
  /*
   * Load M matrix:
   */
  FILE *fptMData = fopen(FileName, "r");
  if (fptMData == NULL) {
    perror("Error: Cannot open input file for reading");
  }
  int MatrixSize = func_num < 7 ? nx * nx : cf_num * nx * nx;
  double *M = malloc(MatrixSize * sizeof(double));
  if (M == NULL) {
    perror("Error: there is insufficient memory available!");
  } else {
    for (int i = 0; i < MatrixSize; ++i) {
      if (fscanf(fptMData, "%lf", &M[i]) == -1) {
        break;
      }
    }
  }
  fclose(fptMData);
  /*
   * Load shift data:
   */
  if (!strcmp(suite, "basic") || !strcmp(suite, "rot") || !strcmp(suite, "bias")) {
    sprintf(FileName, "%s/shift_data_%d_ns.txt", datapath, func_num);
  } else {
    sprintf(FileName, "%s/shift_data_%d.txt", datapath, func_num);
  }
  FILE *fptOShiftData = fopen(FileName, "r");
  if (fptOShiftData == NULL) {
    perror("Error: Cannot open input file for reading");
  }
  int OShiftSize = func_num < 7 ? nx : cf_num * nx;
  double *OShift = malloc(OShiftSize * sizeof(double));
  if (OShift == NULL) {
    perror("Error: there is insufficient memory available!");
  }
  if (func_num < 7) {
    for (int i = 0; i < OShiftSize; ++i) {
      if (fscanf(fptOShiftData, "%lf", &OShift[i]) == -1) {
        break;
      }
    }
  } else {
    for (int i = 0; i < cf_num - 1; i++) {
      for (int j = 0; j < nx; j++) {
        int count = fscanf(fptOShiftData, "%lf", &OShift[i * nx + j]);
        if (count == -1) {
          break;
        }
      }
      int count = fscanf(fptOShiftData, "%*[^\n]%*c");
      if (count == -1) {
        break;
      }
    }
    for (int j = 0; j < nx; j++) {
      if (fscanf(fptOShiftData, "%lf", &OShift[nx * (cf_num - 1) + j]) == -1) {
        break;
      }
    }
  }
  fclose(fptOShiftData);
  /*
   * Load shuffle data:
   */
  int *SS = NULL;
  int const shuffleFlag = (func_num >= 5 && func_num <= 7) ? 1 : 0;
  if (shuffleFlag) {
    sprintf(FileName, "%s/shuffle_data_%d_D%d.txt", datapath, func_num, nx);
    FILE *fptShuffleData = fopen(FileName, "r");
    if (fptShuffleData == NULL) {
      perror("Error: Cannot open input file for reading");
    }
    int ShuffleSize = nx;
    SS = malloc(ShuffleSize * sizeof(int));
    for (int i = 0; i < ShuffleSize; ++i) {
      if (fscanf(fptShuffleData, "%d", &SS[i]) == -1) {
        perror("Cannot read shuffle data");
      }
    }
    fclose(fptShuffleData);
  }

  int shiftFlag = 0;
  if (!strcmp(suite, "bias_shift") || !strcmp(suite, "shift") ||
      !strcmp(suite, "bias_shift_rot") || !strcmp(suite, "shift_rot")) {
    shiftFlag = 1;
  }

  int biasFlag = 0;
  if (!strcmp(suite, "bias") || !strcmp(suite, "bias_shift") ||
      !strcmp(suite, "bias_rot") || !strcmp(suite, "bias_shift_rot")) {
    biasFlag = 1;
  }

  for (int i = 0; i < mx; i++) {
    switch (func_num) {
    case 1:
      cec2021_bent_cigar_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      break;
    case 2:
      cec2021_schwefel_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      break;
    case 3:
      cec2021_bi_rastrigin_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      break;
    case 4:
      cec2021_grie_rosen_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      break;
    case 5:
      cec2021_hf01(&x[i * nx], &f[i], nx, OShift, M, SS, 1, 1);
      break;
    case 6:
      cec2021_hf02(&x[i * nx], &f[i], nx, OShift, M, SS, 1, 1);
      break;
    case 7:
      cec2021_hf03(&x[i * nx], &f[i], nx, OShift, M, SS, 1, 1);
      break;
    case 8:
      shiftFlag ? cec2021_cf01_s(&x[i * nx], &f[i], nx, OShift, M, 1)
                : cec2021_cf01(&x[i * nx], &f[i], nx, OShift, M, 1);
      break;
    case 9:
      shiftFlag ? cec2021_cf02_s(&x[i * nx], &f[i], nx, OShift, M, 1)
                : cec2021_cf02(&x[i * nx], &f[i], nx, OShift, M, 1);
      break;
    case 10:
      shiftFlag ? cec2021_cf03_s(&x[i * nx], &f[i], nx, OShift, M, 1)
                : cec2021_cf03(&x[i * nx], &f[i], nx, OShift, M, 1);
      break;
    default:
      perror("Error: There are only 10 test functions in this test suite!");
      f[i] = 0.0;
      break;
    }
  }
  f[0] += getFunctionBias(biasFlag, func_num);
  free(M);
  free(OShift);
  if (shuffleFlag) {
    free(SS);
  }
}
