#ifndef INTERFACE_H_
#define INTERFACE_H_

#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <string.h>

#include "cec2021_functions.h"
#include "interface_basic.h"
#include "interface_shift.h"
#include "interface_rot.h"
#include "interface_bias.h"
#include "interface_bias_rot.h"
#include "interface_bias_shift.h"
#include "interface_shift_rot.h"
#include "interface_bias_shift_rot.h"

extern double *OShift, *M, *y, *z, *x_bound;
extern int ini_flag, n_flag, func_flag, *SS;
extern char *extdata;

void cec21_func(double *, double *, int, int, int, char *);

#endif // INTERFACE_H_
