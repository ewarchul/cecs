/*
  CEC21 Test Function Suite for Single Objective Optimization
  Copyright 2020 Eryk Warchulski ewarchul@gmail.com

  Based on: https://github.com/P-N-Suganthan/2021-SO-BCO
*/

#include "cec2021.h"

void cec2021(char **extdatadir, int *i, double *X, int *row, int *col,
             double *f, char *suite) {
  int r, c;
  double *x;

  extdata = *extdatadir;

  x = (double *)malloc(*col * sizeof(double));

  for (r = 0; r < *row; r++) {
     R_CheckUserInterrupt();

    for (c = 0; c < *col; c++) {
      x[c] = X[r + *row * c];
    }
    cec2021_func(x, &f[r], *col, 1, *i, suite);
  }

  free(x);
}
