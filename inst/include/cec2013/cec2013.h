#ifndef CEC2013_H_
#define CEC2013_H_

#include <R.h>
#include <stdio.h>
#include <stdlib.h>

#include "cec2013_interface.h"

void cec2013(char **extdatadir, int *i, double *X, int *row, int *col,
             double *f);

#endif //CEC2013_H_
