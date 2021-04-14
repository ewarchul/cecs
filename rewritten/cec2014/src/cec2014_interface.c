#include "../include/cec2014_interface.h"

/**
 * An interface for CEC 2017 functions.
 *
 * @param x
 * @param f
 * @param nx
 * @param mx
 * @param func_num
 */

void cec2014_func(char *datapath, double *x, double *f, int nx, int mx,
                  int func_num) {
  if (!(nx == 2 || nx == 10 || nx == 20 || nx == 30 || nx == 50 || nx == 100)) {
    perror("Error: Test functions are only defined for D = 2, 10, 20, 30, 50, "
           "100.");
  }
  if (nx == 2 && ((func_num >= 17 && func_num <= 22) ||
                  (func_num >= 29 && func_num <= 30))) {
    perror("Error: hf0{1..6}, cf0{7..8} are NOT defined for D=2.");
  }
  int cf_num = 10;
  char FileName[256];
  /*
   * Load M matrix:
   */
  sprintf(FileName, "%s/M_%d_D%d.txt", datapath, func_num, nx);
  FILE *fptMData = fopen(FileName, "r");
  if (fptMData == NULL) {
    perror("Error: Cannot open input file for reading");
  }
  int MatrixSize = func_num < 23 ? nx * nx : cf_num * nx * nx;
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
   * Load shift data:
   */
  sprintf(FileName, "%s/shift_data_%d.txt", datapath, func_num);
  FILE *fptOShiftData = fopen(FileName, "r");
  if (fptOShiftData == NULL) {
    perror("Error: Cannot open input file for reading");
  }
  int OShiftSize = func_num < 20 ? nx : cf_num * nx;
  double *OShift = malloc(OShiftSize * sizeof(double));
  if (OShift == NULL) {
    perror("Error: there is insufficient memory available!");
  }
  if (func_num < 23) {
    for (int i = 0; i < OShiftSize; ++i) {
      if (fscanf(fptOShiftData, "%lf", &OShift[i]) == -1) {
        perror("Cannot read shift data.");
      }
    }
  } else {
    for (int i = 0; i < cf_num - 1; i++) {
      for (int j = 0; j < nx; j++) {
        int count = fscanf(fptOShiftData, "%lf", &OShift[i * nx + j]);
        if (count == -1) {
          perror("Cannot read shift data.");
        }
      }
      int count = fscanf(fptOShiftData, "%*[^\n]%*c");
      if (count == -1) {
        perror("Cannot read shift data.");
      }
    }
    for (int j = 0; j < nx; j++) {
      if (fscanf(fptOShiftData, "%lf", &OShift[(cf_num - 1) * nx + j]) == -1) {
        perror("Cannot read shift data.");
      }
    }
  }
  fclose(fptOShiftData);
  /*
   * Load shuffle data:
   */
  int *SS = NULL;
  int shuffleFlag =
      ((func_num >= 17 && func_num <= 22) || (func_num == 29 || func_num == 30))
          ? 1
          : 0;
  if (shuffleFlag) {
    sprintf(FileName, "%s/shuffle_data_%d_D%d.txt", datapath, func_num, nx);
    FILE *fptShuffleData = fopen(FileName, "r");
    if (fptShuffleData == NULL) {
      perror("Error: Cannot open input file for reading");
    }
    int ShuffleSize = (func_num >= 11 && func_num <= 20) ? nx : cf_num * nx;
    SS = malloc(ShuffleSize * sizeof(int));
    for (int i = 0; i < ShuffleSize; ++i) {
      if (fscanf(fptShuffleData, "%d", &SS[i]) == -1) {
        perror("Cannot read shuffle data");
      }
    }
  }
  for (int i = 0; i < mx; i++) {
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
      perror("Error: There are only 30 test functions in this test suite [CEC2014]!");
      f[i] = 0.0;
      break;
    }
  }
  free(M);
  free(OShift);
  if (shuffleFlag) {
    free(SS);
  }
}
