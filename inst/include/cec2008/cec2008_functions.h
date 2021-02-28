#ifndef CEC2008_FUNCTIONS_H_
#define CEC2008_FUNCTIONS_H_

#include <math.h>
#include <stdlib.h>

#include "../consts.h"
#include "../globals.h"

#define abss(a) (a < 0 ? (-a) : a)
#define max(a, b) (a > b ? a : b)
#define min(a, b) (a > b ? b : a)

void cec2008_Shifted_Sphere(double *, double *, int);
//double cec2008_Schwefel_Problem(int dim, double *x);
//double cec2008_Shifted_Rosenbrock(int dim, double *x);
//double cec2008_Shifted_Rastrigin(int dim, double *x);
//double cec2008_Shifted_Griewank(int dim, double *x);
//double cec2008_Shifted_Ackley(int dim, double *x);
//double cec2008_Shifted_Weierstrass(int dim, double *x);

#endif // CEC2008_FUNCTIONS_H_
