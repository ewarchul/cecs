/*
  CEC15 Test Function Suite for Single Objective Optimization
  Copyright 2021 Eryk Warchulski ewarchul@gmail.com

  Based on: https://github.com/P-N-Suganthan/CEC2014
*/

#include "cec2015/cec2015.h"

void cec2015(char **extdatadir, int *i, double *X, int *row, int *col,
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
    cec2015_func(x, &f[r], *col, 1, *i);
  }

  free(x);
}
