#ifndef CEC2008_H_
#define CEC2008_H_

#include <R.h>
#include <stdio.h>
#include <stdlib.h>

#include "cec2008_interface.h"

void cec2008(char **extdatadir, int *i, double *X, int *row, int *col,
             double *f);

#endif //CEC2008_H_
