#include "cec.h"

void cec(char **extdatadir, int *cecVersion, int *i, double *X, int *row,
         int *col, double *f, char *suite) {

  double *x = malloc(*col * sizeof(double));

  for (int r = 0; r < *row; r++) {
    for (int c = 0; c < *col; c++) {
      x[c] = X[r + *row * c];
    }
    switch (*cecVersion) {
    case 2014:
      cec2014_interface(*extdatadir, x, &f[r], *col, 1, *i);
      break;
    case 2015:
      cec2015_interface(*extdatadir, x, &f[r], *col, 1, *i);
      break;
    case 2017:
      cec2017_interface(*extdatadir, x, &f[r], *col, 1, *i);
      break;
    case 2019:
      cec2019_interface(*extdatadir, x, &f[r], *col, 1, *i);
      break;
    case 2021:
      cec2021_interface(*extdatadir, x, &f[r], *col, 1, *i, suite);
      break;
    default:
      perror("Not implemented");
      break;
    }
  }
  free(x);
}
