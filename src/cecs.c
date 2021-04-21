#include "cecs.h"

CecData cd = {
  .prevDimension = 0,
  .prevFunction = 0,
  .dataLoaded = 0,
};

void cecs(char **extdatadir, char **suite, char *cec, int *problem,
          double *input, int *row, int *col, double *output) {

  double *x = malloc(*col * sizeof(double));

  for (int r = 0; r < *row; r++) {
    R_CheckUserInterrupt();
    for (int c = 0; c < *col; c++) {
      x[c] = input[r + *row * c];
      switch (*cec) {
      case 14:
        cec2014_interface(*extdatadir, x, output, *col, *row, *problem);
        break;
      case 15:
        cec2015_interface(*extdatadir, x, output, *col, *row, *problem);
        break;
      case 17:
        cec2017_interface(*extdatadir, x, output, *col, *row, *problem);
        break;
      case 19:
        cec2019_interface(*extdatadir, x, output, *col, *row, *problem);
        break;
      case 21:
        cec2021_interface(*extdatadir, x, output, *col, *row, *problem, *suite);
        break;
      }
    }
  }
  free(x);
}
