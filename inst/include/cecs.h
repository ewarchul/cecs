#ifndef CECS_H
#define CECS_H

#include "interfaces.h"
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <R.h>

void cecs(char **extdatadir, char **suite, char *cec, int *i, double *X,
          int *row, int *col, double *f);


#endif // CECS_H
