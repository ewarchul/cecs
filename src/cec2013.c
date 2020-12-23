/* File cec2013.c
Part of the cec2013 R package, http://www.rforge.net/cec2013/ ;
                               http://cran.r-project.org/web/packages/cec2013
Copyright 2013 Yasser Gonzalez-Fernandez & Mauricio Zambrano-Bigiarini
Distributed under GPL 3 or later
*/

#include "cec2013.h"


void cec2013(char **extdatadir, int *i, double *X, int *row, int *col,
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
    cec2013_func(x, &f[r], *col, 1, *i);
  }

  free(x);
}
