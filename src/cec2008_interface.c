#include "cec2008/cec2008_interface.h"

void cec2008_func(double *x, double *f, int nx, int func_num) {

  if (nx > 1000) {
    nx = 1000;
  }

  switch (func_num) {
  case 1:
    cec2008_Shifted_Sphere(x, f, nx);
    break;
  /*case 2:*/
    /*cec2008_Schwefel_Problem(nx, x);*/
    /*break;*/
  /*case 3:*/
    /*cec2008_Shifted_Rosenbrock(nx, x);*/
    /*break;*/
  /*case 4:*/
    /*cec2008_Shifted_Rastrigin(nx, x);*/
    /*break;*/
  /*case 5:*/
    /*cec2008_Shifted_Griewank(nx, x);*/
    /*break;*/
  /*case 6:*/
    /*cec2008_Shifted_Ackley(nx, x);*/
    /*break;*/
  default:
    perror("Function number out of range (1-6)");
    exit(0);
    break;
  }
}
