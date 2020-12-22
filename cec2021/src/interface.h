#ifndef INTERFACE_H_
#define INTERFACE_H_

#include "interface_basic.h"
#include "interface_shift.h"
#include "interface_rot.h"
#include "interface_bias.h"
#include "interface_bias_rot.h"
#include "interface_bias_shift.h"
#include "interface_shift_rot.h"
#include "interface_bias_shift_rot.h"

void cec21_func(double *, double *, int, int, int, char *);

#endif // INTERFACE_H_