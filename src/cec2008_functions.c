#include "cec2008/cec2008_functions.h"
#include "cec2008/cec2008_data.h"

void cec2008_Shifted_Sphere(double *x, double *f, int nx) {
  int i;
  double z;
  f[0] = 0.0;
  for (i = 0; i < nx; i++) {
    z = x[i] - cec2008_sphere[i];
    f[0] += z * z;
  }
}

/*double cec2008_Schwefel_Problem(int dim, double *x) {*/
  /*int i;*/
  /*double z;*/
  /*double F = abss(x[0] - schwefel[0]);*/
  /*for (i = 1; i < dim; i++) {*/
    /*z = x[i] - schwefel[i];*/
    /*F = max(F, abss(z));*/
  /*}*/
  /*return F + f_bias[1];*/
/*}*/

/*double cec2008_Shifted_Rosenbrock(int dim, double *x) {*/
  /*int i;*/
  /*double z[dim];*/
  /*double F = 0;*/

  /*for (i = 0; i < dim; i++)*/
    /*z[i] = x[i] - rosenbrock[i] + 1;*/

  /*for (i = 0; i < dim - 1; i++) {*/
    /*F = F + 100 * (pow((pow(z[i], 2) - z[i + 1]), 2)) + pow((z[i] - 1), 2);*/
  /*}*/
  /*return F + f_bias[2];*/
/*}*/

/*double cec2008_Shifted_Rastrigin(int dim, double *x) {*/
  /*int i;*/
  /*double z;*/
  /*double F = 0;*/
  /*for (i = 0; i < dim; i++) {*/
    /*z = x[i] - rastrigin[i];*/
    /*F = F + (pow(z, 2) - 10 * cos(2 * M_PI * z) + 10);*/
  /*}*/
  /*return F + f_bias[3];*/
/*}*/

/*double cec2008_Shifted_Griewank(int dim, double *x) {*/
  /*int i;*/
  /*double z;*/
  /*double F1 = 0;*/
  /*double F2 = 1;*/
  /*for (i = 0; i < dim; i++) {*/
    /*z = x[i] - griewank[i];*/
    /*F1 = F1 + (pow(z, 2) / 4000);*/
    /*F2 = F2 * (cos(z / sqrt(i + 1)));*/
  /*}*/
  /*return (F1 - F2 + 1 + f_bias[4]);*/
/*}*/

/*double cec2008_Shifted_Ackley(int dim, double *x) {*/
  /*int i;*/
  /*double z;*/
  /*double Sum1 = 0;*/
  /*double Sum2 = 0;*/
  /*double F = 0;*/
  /*for (i = 0; i < dim; i++) {*/
    /*z = x[i] - ackley[i];*/
    /*Sum1 = Sum1 + pow(z, 2);*/
    /*Sum2 = Sum2 + cos(2 * M_PI * z);*/
  /*}*/
  /*F = -20 * exp(-0.2 * sqrt(Sum1 / dim)) - exp(Sum2 / dim) + 20 + E + f_bias[5];*/

  /*return F;*/
/*}*/
