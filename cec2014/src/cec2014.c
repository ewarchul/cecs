#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include "cec14_test_func.h"

void cec14_test_func(double *x, double *f, int nx, int mx, int func_num);

double *OShift, *M, *y, *z, *x_bound;
int ini_flag = 0, n_flag, func_flag, *SS;
char *extdata;

void cec2014(char **extdatadir, int *i, double *X, int *row, int *col,
             double *f) {
  int r, c;
  double *x;

  extdata = *extdatadir;

  x = (double *)malloc(*col * sizeof(double));

  for (r = 0; r < *row; r++) {
     R_CheckUserInterrupt();

    for (c = 0; c < *col; c++) {
      x[c] = X[r + *row * c];
    }
    cec14_test_func(x, &f[r], *col, 1, *i);
  }

  free(x);
}