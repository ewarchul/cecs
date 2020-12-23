#ifndef CEC2021_INTERFACE_H_
#define CEC2021_INTERFACE_H_

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cec2021_functions.h"

extern double *OShift, *M, *y, *z, *x_bound;
extern int ini_flag, n_flag, func_flag, *SS;
extern char *extdata;

void cec2021_func(double *, double *, int, int, int, char *);
void cec2021_basic_func(double *, double *, int, int, int);
void cec2021_shift_func(double *, double *, int, int, int);
void cec2021_bias_func(double *, double *, int, int, int);
void cec2021_rot_func(double *, double *, int, int, int);
void cec2021_bias_rot_func(double *, double *, int, int, int);
void cec2021_bias_shift_func(double *, double *, int, int, int);
void cec2021_bias_shift_rot_func(double *, double *, int, int, int);
void cec2021_shift_rot_func(double *, double *, int, int, int);

#endif // CEC2021_INTERFACE_H_
