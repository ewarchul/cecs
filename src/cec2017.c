/*
  CEC17 Test Function Suite for Single Objective Optimization
  Copyright 2017 Dariusz Jagodziński <d.jagodzinski@elka.pw.edu.pl>

  Based on:
  http://web.mysites.ntu.edu.sg/epnsugan/PublicSite/Shared%20Documents/CEC-2017/Bound-Constrained/code.rar
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>

#include <R.h>

void cec2017_func(double *, double *,int,int,int);

double *OShift,*M,*y,*z,*x_bound;
int ini_flag=0,n_flag,func_flag,*SS;
char *extdata;

void cec2017(char **extdatadir, int *i, double *X, int *row, int *col, double *f)
{
  int r, c;
  double *x;

  extdata = *extdatadir;

  x = (double *) malloc(*col * sizeof(double));

  for (r = 0; r < *row; r++) {
    R_CheckUserInterrupt();

    for (c = 0; c < *col; c++) {
      x[c] = X[r + *row * c];
    }
    cec17_test_func(x, &f[r] , *col, 1, *i);
  }

  free(x);

}
