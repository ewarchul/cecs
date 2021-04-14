#include "../include/cec2019.h"
#include "../include/cec2019_interface.h"

void cec2019(char **extdatadir, int *i, double *X, int *row, int *col,
             double *f) {

  double *x = malloc(*col * sizeof(double));

  for (int r = 0; r < *row; r++) {
    for (int c = 0; c < *col; c++) {
      x[c] = X[r + *row * c];
    }
    cec2019_func(*extdatadir, x, &f[r], *col, 1, *i);
  }
  free(x);
}
