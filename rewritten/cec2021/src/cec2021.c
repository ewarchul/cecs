#include "../include/cec2021.h"
#include "../include/cec2021_interface.h"

void cec2021(char **extdatadir, int *i, double *X, int *row, int *col,
             double *f, char *suite) {

  double *x = malloc(*col * sizeof(double));

  for (int r = 0; r < *row; r++) {
    for (int c = 0; c < *col; c++) {
      x[c] = X[r + *row * c];
    }
    cec2021_func(*extdatadir, x, &f[r], *col, 1, *i, suite);
  }
  free(x);
}
