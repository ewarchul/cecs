/* File cec2013.c
Part of the cec2013 R package, http://www.rforge.net/cec2013/ ; 
                               http://cran.r-project.org/web/packages/cec2013
Copyright 2013 Yasser Gonzalez-Fernandez & Mauricio Zambrano-Bigiarini
Distributed under GPL 3 or later
*/

#include <stdio.h>
#include <stdlib.h>

#include <R.h>

void test_func(double *x, double *f, int nx, int mx, int func_num);

double *OShift, *M, *y, *z, *x_bound;
int ini_flag = 0, n_flag, func_flag;
char *extdata;

void cec2013(char **extdatadir, int *i, double *X, int *row, int *col, double *f)
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
        test_func(x, &f[r] , *col, 1, *i);
    }

    free(x);
}
