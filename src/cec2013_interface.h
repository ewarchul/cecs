#ifndef CEC2013_INTERFACE_H_
#define CEC2013_INTERFACE_H_

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cec2013_functions.h"

extern double *OShift, *M, *y, *z, *x_bound;
extern int ini_flag, n_flag, func_flag, *SS;
extern char *extdata;

void cec2013_func(double *, double *, int, int, int);

#endif // CEC2013_INTERFACE_H_