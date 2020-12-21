/*
  CEC20 Test Function Suite for Single Objective Bound Constrained Numerical Optimization 
 

*/

The test functions are as follows for different case
1) cec21_bias_shift_rot_func(double *, double *,int,int,int);
2) cec21_bias_shift_func(double *, double *,int,int,int);
3) cec21_bias_rot_func(double *, double *,int,int,int);
4) cec21_shift_rot_func(double *, double *,int,int,int);
5) cec21_rot_func(double *, double *,int,int,int);
6) cec21_shift_func(double *, double *,int,int,int);
7) cec21_bias_func(double *, double *,int,int,int);
8) cec21_basic_func(double *, double *,int,int,int);
Example:
test_func(x, f, dimension,population_size,func_num);

where test_func according to case.

main.cpp is an example function about how to use test_func.cpp


#include <WINDOWS.H>    
#include <stdio.h>
#include <math.h>
#include <malloc.h>
void test_func(double *, double *,int,int,int);
double *OShift,*M,*y,*z,*x_bound;
int ini_flag=0,n_flag,func_flag;
void main()
{
...
}

lshade algorithm has been used for example.

For Linux Users:
Please  change %xx in fscanf and fprintf and do use "WINDOWS.H". 
