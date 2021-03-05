#include "cecs.h"

void cecs(char **extdatadir, char **suite, char *cec, int *i, double *X,
          int *row, int *col, double *f) {
  switch (*cec) {
  case 13:
    cec2013(extdatadir, i, X, row, col, f);
    break;
  case 14:
    cec2014(extdatadir, i, X, row, col, f);
    break;
  case 15:
    cec2015(extdatadir, i, X, row, col, f);
    break;
  case 17:
    cec2017(extdatadir, i, X, row, col, f);
    break;
  case 21:
    cec2021(extdatadir, suite, i, X, row, col, f);
    break;
  }
}
