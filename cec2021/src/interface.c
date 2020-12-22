#include "interface.h"
#include <string.h>

extern double *OShift, *M, *y, *z, *x_bound;
extern int ini_flag, n_flag, func_flag, *SS;
extern char *extdata;

void cec21_func(double *x, double *f, int nx, int mx, int func_num,
                char *suite) {
  if (strcmp(suite, "basic")) {

    cec21_basic_func(x, f, nx, mx, func_num);

  } else if (strcmp(suite, "shift")) {

    cec21_shift_func(x, f, nx, mx, func_num);

  } else if (strcmp(suite, "rot")) {

    cec21_rot_func(x, f, nx, mx, func_num);

  } else if (strcmp(suite, "bias")) {

    cec21_bias_func(x, f, nx, mx, func_num);

  } else if (strcmp(suite, "shift_rot")) {

    cec21_shift_rot_func(x, f, nx, mx, func_num);

  } else if (strcmp(suite, "bias_rot")) {

    cec21_bias_rot_func(x, f, nx, mx, func_num);

  } else if (strcmp(suite, "bias_shift")) {

    cec21_bias_shift_func(x, f, nx, mx, func_num);

  } else if (strcmp(suite, "bias_shift_rot")) {

    cec21_shift_rot_func(x, f, nx, mx, func_num);

  } else {
    return;
  }
}