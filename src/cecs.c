/*
  CEC Test Function Suite for Single Objective Optimization
  Copyright 2020 Eryk Warchulski ewarchul@gmail.com
*/

#include "cecs.h"

void cecs(char **extdatadir, char **cec, int *i, double *X, int *row, int *col,
             double *f, char **suite) {
    char *cecx = *cec;
    if (!strcmp(cecx, "cec2013")) {

        cec2013(extdatadir, i, X, row, col, f);

    } else if (!strcmp(cecx, "cec2014") != 0) {

        cec2014(extdatadir, i, X, row, col, f);

    } else if (!strcmp(cecx, "cec2017")) {

        cec2017(extdatadir, i, X, row, col, f);

    } else if (!strcmp(cecx, "cec2021")) {

        cec2021(extdatadir, i, X, row, col, f, suite);

    } else {

        return;
    }
}
