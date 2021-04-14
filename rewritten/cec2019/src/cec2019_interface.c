#include "../include/cec2019_interface.h"

/**
 * An interface for CEC 2019 functions.
 *
 * @param x
 * @param f
 * @param nx
 * @param mx
 * @param func_num
 */

void cec2019_func(char *datapath, double *x, double *f, int nx, int mx,
                  int func_num) {
  if (!(nx == 2 || nx == 9 || nx == 10 || nx == 16 || nx == 18)) {
    perror("Error: Test functions are only defined for D=10, 9, 16, 18 \n\
          F1 is defined on D=9 \n F2 is defined on D=16 \n\
          F3 is defined on D=18 \n F4-F10 are defined on D=10.");
  }
  char FileName[256];
  /*
   * Load M matrix:
   */
  int externalDataFlag = func_num > 3 ? 1 : 0;
  double *M = NULL;
  double *OShift = NULL;
  if (externalDataFlag) {
    sprintf(FileName, "%s/M_%d_D%d.txt", datapath, func_num, nx);
    FILE *fptMData = fopen(FileName, "r");
    if (fptMData == NULL) {
      perror("Error: Cannot open input file for reading");
    }
    int MatrixSize = nx * nx;
    M = malloc(MatrixSize * sizeof(double));
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
    int OShiftSize = nx;
    OShift = malloc(OShiftSize * sizeof(double));
    if (OShift == NULL) {
      perror("Error: there is insufficient memory available!");
    }
    for (int i = 0; i < OShiftSize; ++i) {
      if (fscanf(fptOShiftData, "%lf", &OShift[i]) == -1) {
        perror("Cannot read shift  data.");
      }
    }
    fclose(fptOShiftData);
  }
  for (int i = 0; i < mx; i++) {
    switch (func_num) {
    case 1:
      Chebyshev(&x[i * nx], nx, &f[i]);
      f[i] += 1.0;
      break;
    case 2:
      Hilbert(&x[i * nx], nx, &f[i]);
      f[i] += 1.0;
      break;
    case 3:
      Lennard_Jones(&x[i * nx], nx, &f[i]);
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
      perror("Error: There are only 10 test functions in this test suite! [CEC2015-LB]");
      f[i] = 0.0;
      break;
    }
  }
  if (externalDataFlag) {
    free(M);
    free(OShift);
  }
}
