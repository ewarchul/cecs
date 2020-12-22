#include <stdio.h>
#include <stdlib.h>
#include <R.h>

#include "interface.h"

double *OShift, *M, *y, *z, *x_bound;
int ini_flag = 0, n_flag, func_flag, *SS;
char *extdata;

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
    cec21_func(x, &f[r], *col, 1, *i, suite);
  }

  free(x);
}
