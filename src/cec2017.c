/*
  CEC17 Test Function Suite for Single Objective Optimization
  Copyright 2017 Dariusz Jagodziński <d.jagodzinski@elka.pw.edu.pl>

  Based on:
  http://web.mysites.ntu.edu.sg/epnsugan/PublicSite/Shared%20Documents/CEC-2017/Bound-Constrained/code.rar
*/

#include "cec2017.h"



void cec2017(char **extdatadir, int *i, double *X, int *row, int *col,
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
    cec2017_func(x, &f[r], *col, 1, *i);
  }

  free(x);
}
