#ifndef CEC2019_FUNCTIONS_H_
#define CEC2019_FUNCTIONS_H_

#include "cec2019/cec2019_functions.h"

void cec2019_schaffer_F7_func(double *x, double *f, int nx, double *Os,
                              double *Mr, int s_flag,
                              int r_flag) /* Schwefel's 1.2  */
{
  int i;
  double tmp;
  f[0] = 0.0;
  cec2019_sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag); /* shift and rotate */
  for (i = 0; i < nx - 1; i++) {
    z[i] = pow(y[i] * y[i] + y[i + 1] * y[i + 1], 0.5);
    tmp = sin(50.0 * pow(z[i], 0.2));
    f[0] += pow(z[i], 0.5) + pow(z[i], 0.5) * tmp * tmp;
  }
  f[0] = f[0] * f[0] / (nx - 1) / (nx - 1);
}

void cec2019_griewank_func(double *x, double *f, int nx, double *Os, double *Mr,
                           int s_flag, int r_flag) /* Griewank's  */
{
  int i;
  double s, p;
  s = 0.0;
  p = 1.0;

  cec2019_sr_func(x, z, nx, Os, Mr, 600.0 / 100.0, s_flag,
                  r_flag); /* shift and rotate */

  for (i = 0; i < nx; i++) {
    s += z[i] * z[i];
    p *= cos(z[i] / sqrt(1.0 + i));
  }
  f[0] = 1.0 + s / 4000.0 - p;
}

void cec2019_ackley_func(double *x, double *f, int nx, double *Os, double *Mr,
                         int s_flag, int r_flag) /* Ackley's  */
{
  int i;
  double sum1, sum2;
  sum1 = 0.0;
  sum2 = 0.0;

  cec2019_sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag); /* shift and rotate */

  for (i = 0; i < nx; i++) {
    sum1 += z[i] * z[i];
    sum2 += cos(2.0 * M_PI * z[i]);
  }
  sum1 = -0.2 * sqrt(sum1 / nx);
  sum2 /= nx;
  f[0] = E - 20.0 * exp(sum1) - exp(sum2) + 20.0;
}

void cec2019_weierstrass_func(double *x, double *f, int nx, double *Os,
                              double *Mr, int s_flag,
                              int r_flag) /* Weierstrass's  */
{
  int i, j, k_max;
  double sum, sum2, a, b;
  a = 0.5;
  b = 3.0;
  k_max = 20;
  f[0] = 0.0;

  cec2019_sr_func(x, z, nx, Os, Mr, 0.5 / 100.0, s_flag,
                  r_flag); /* shift and rotate */

  for (i = 0; i < nx; i++) {
    sum = 0.0;
    sum2 = 0.0;
    for (j = 0; j <= k_max; j++) {
      sum += pow(a, j) * cos(2.0 * M_PI * pow(b, j) * (z[i] + 0.5));
      sum2 += pow(a, j) * cos(2.0 * M_PI * pow(b, j) * 0.5);
    }
    f[0] += sum;
  }
  f[0] -= nx * sum2;
}

void cec2019_rastrigin_func(double *x, double *f, int nx, double *Os, double *Mr,
                    int s_flag, int r_flag) /* Rastrigin's  */
{
  int i;
  f[0] = 0.0;

  cec2019_sr_func(x, z, nx, Os, Mr, 5.12 / 100.0, s_flag,
          r_flag); /* shift and rotate */

  for (i = 0; i < nx; i++) {
    f[0] += (z[i] * z[i] - 10.0 * cos(2.0 * M_PI * z[i]) + 10.0);
  }
}

void cec2019_step_rastrigin_func(double *x, double *f, int nx, double *Os, double *Mr,
                         int s_flag,
                         int r_flag) /* Noncontinuous Rastrigin's  */
{
  int i;
  f[0] = 0.0;
  for (i = 0; i < nx; i++) {
    if (fabs(y[i] - Os[i]) > 0.5)
      y[i] = Os[i] + floor(2 * (y[i] - Os[i]) + 0.5) / 2;
  }

  cec2019_sr_func(x, z, nx, Os, Mr, 5.12 / 100.0, s_flag,
          r_flag); /* shift and rotate */

  for (i = 0; i < nx; i++) {
    f[0] += (z[i] * z[i] - 10.0 * cos(2.0 * M_PI * z[i]) + 10.0);
  }
}

void cec2019_schwefel_func(double *x, double *f, int nx, double *Os, double *Mr,
                   int s_flag, int r_flag) /* Schwefel's  */
{
  int i;
  double tmp;
  f[0] = 0.0;

  cec2019_sr_func(x, z, nx, Os, Mr, 1000.0 / 100.0, s_flag,
          r_flag); /* shift and rotate */

  for (i = 0; i < nx; i++) {
    z[i] += 4.209687462275036e+002;
    if (z[i] > 500) {
      f[0] -=
          (500.0 - fmod(z[i], 500)) * sin(pow(500.0 - fmod(z[i], 500), 0.5));
      tmp = (z[i] - 500.0) / 100;
      f[0] += tmp * tmp / nx;
    } else if (z[i] < -500) {
      f[0] -= (-500.0 + fmod(fabs(z[i]), 500)) *
              sin(pow(500.0 - fmod(fabs(z[i]), 500), 0.5));
      tmp = (z[i] + 500.0) / 100;
      f[0] += tmp * tmp / nx;
    } else
      f[0] -= z[i] * sin(pow(fabs(z[i]), 0.5));
  }
  f[0] += 4.189828872724338e+002 * nx;
}

void cec2019_escaffer6_func(double *x, double *f, int nx, double *Os, double *Mr,
                    int s_flag, int r_flag) /* Expanded Scaffer¡¯s F6  */
{
  int i;
  double temp1, temp2;

  cec2019_sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag); /* shift and rotate */

  f[0] = 0.0;
  for (i = 0; i < nx - 1; i++) {
    temp1 = sin(sqrt(z[i] * z[i] + z[i + 1] * z[i + 1]));
    temp1 = temp1 * temp1;
    temp2 = 1.0 + 0.001 * (z[i] * z[i] + z[i + 1] * z[i + 1]);
    f[0] += 0.5 + (temp1 - 0.5) / (temp2 * temp2);
  }
  temp1 = sin(sqrt(z[nx - 1] * z[nx - 1] + z[0] * z[0]));
  temp1 = temp1 * temp1;
  temp2 = 1.0 + 0.001 * (z[nx - 1] * z[nx - 1] + z[0] * z[0]);
  f[0] += 0.5 + (temp1 - 0.5) / (temp2 * temp2);
}

void cec2019_happycat_func(
    double *x, double *f, int nx, double *Os, double *Mr, int s_flag,
    int r_flag) /* HappyCat, provdided by Hans-Georg Beyer (HGB) */
{
  int i;
  double alpha, r2, sum_z;
  alpha = 1.0 / 8.0;

  cec2019_sr_func(x, z, nx, Os, Mr, 5.0 / 100.0, s_flag, r_flag); /* shift and rotate */

  r2 = 0.0;
  sum_z = 0.0;
  for (i = 0; i < nx; i++) {
    z[i] = z[i] - 1.0; // shift to orgin
    r2 += z[i] * z[i];
    sum_z += z[i];
  }
  f[0] = pow(fabs(r2 - nx), 2 * alpha) + (0.5 * r2 + sum_z) / nx + 0.5;
}

void cec2019_shiftfunc(double *x, double *xshift, int nx, double *Os) {
  int i;
  for (i = 0; i < nx; i++) {
    xshift[i] = x[i] - Os[i];
  }
}

void cec2019_rotatefunc(double *x, double *xrot, int nx, double *Mr) {
  int i, j;
  for (i = 0; i < nx; i++) {
    xrot[i] = 0;
    for (j = 0; j < nx; j++) {
      xrot[i] = xrot[i] + x[j] * Mr[i * nx + j];
    }
  }
}

void cec2019_sr_func(double *x, double *sr_x, int nx, double *Os, double *Mr,
             double sh_rate, int s_flag, int r_flag) /* shift and rotate */
{
  int i;
  if (s_flag == 1) {
    if (r_flag == 1) {
      cec2019_shiftfunc(x, y, nx, Os);
      for (i = 0; i < nx; i++) // shrink to the orginal search range
      {
        y[i] = y[i] * sh_rate;
      }
      cec2019_rotatefunc(y, sr_x, nx, Mr);
    } else {
      cec2019_shiftfunc(x, sr_x, nx, Os);
      for (i = 0; i < nx; i++) // shrink to the orginal search range
      {
        sr_x[i] = sr_x[i] * sh_rate;
      }
    }
  } else {

    if (r_flag == 1) {
      for (i = 0; i < nx; i++) // shrink to the orginal search range
      {
        y[i] = x[i] * sh_rate;
      }
      cec2019_rotatefunc(y, sr_x, nx, Mr);
    } else
      for (i = 0; i < nx; i++) // shrink to the orginal search range
      {
        sr_x[i] = x[i] * sh_rate;
      }
  }
}

void cec2019_asyfunc(double *x, double *xasy, int nx, double beta) {
  int i;
  for (i = 0; i < nx; i++) {
    if (x[i] > 0)
      xasy[i] = pow(x[i], 1.0 + beta * i / (nx - 1) * pow(x[i], 0.5));
  }
}

void cec2019_oszfunc(double *x, double *xosz, int nx) {
  int i, sx;
  double c1, c2, xx;
  for (i = 0; i < nx; i++) {
    if (i == 0 || i == nx - 1) {
      if (x[i] != 0)
        xx = log(fabs(x[i]));
      if (x[i] > 0) {
        c1 = 10;
        c2 = 7.9;
      } else {
        c1 = 5.5;
        c2 = 3.1;
      }
      if (x[i] > 0)
        sx = 1;
      else if (x[i] == 0)
        sx = 0;
      else
        sx = -1;
      xosz[i] = sx * exp(xx + 0.049 * (sin(c1 * xx) + sin(c2 * xx)));
    } else
      xosz[i] = x[i];
  }
}

void cec2019_Lennard_Jones(
    double *x, int D,
    double *f) // find the atomic configuration with minimum energy
{
  /* valid for any dimension, D=3*k, k=2,3,4,...,25.   k is the number of atoms
  in 3-D space constraints: unconstrained type: multi-modal with one global
  minimum; non-separable initial upper bound = 4, initial lower bound = -4
  value-to-reach = minima[k-2]+.0001
  f(x*) = minima[k-2]; see array of minima below; additional minima available at
  the Cambridge cluster database:
  http://www-wales.ch.cam.ac.uk/~jon/structures/LJ/tables.150.html
  */

  f[0] = 0;
  int i, j, k, a, b;
  long double xd, yd, zd, ed, ud, sum = 0;

  static double minima[] = {
      -1.,         -3.,        -6.,        -9.103852,  -12.712062, -16.505384,
      -19.821489,  -24.113360, -28.422532, -32.765970, -37.967600, -44.326801,
      -47.845157,  -52.322627, -56.815742, -61.317995, -66.530949, -72.659782,
      -77.1777043, -81.684571, -86.809782, -02.844472, -97.348815, -102.372663};

  k = D / 3;
  if (k < 2) // default if k<2
  {
    k = 2;
    D = 6;
  }

  for (i = 0; i < k - 1; i++) {
    for (j = i + 1; j < k; j++) {
      a = 3 * i;
      b = 3 * j;
      xd = x[a] - x[b];
      yd = x[a + 1] - x[b + 1];
      zd = x[a + 2] - x[b + 2];
      ed = xd * xd + yd * yd + zd * zd;
      ud = ed * ed * ed;
      if (ud > 1.0e-10)
        sum += (1.0 / ud - 2.0) / ud;
      else
        sum += 1.0e20;
    }
  }

  f[0] += sum;
  f[0] += 12.7120622568;
}

void cec2019_Hilbert(
    double *x, int D,
    double *f) // find the inverse of the (ill-conditioned) Hilbert matrix
{
  /* valid for any dimension, n=k*k, k=2,3,4,...
  constraints: unconstrained
  type: multi-modal with one global minimum; non-separable
  initial upper bound = 2^n, initial lower bound = -2^n
  value-to-reach = f(x*)+1.0e-8
  f(x*) = 0.0; x*={{9,-36,30},{-36,192,-180},{30,-180,180}} (n=9)
  x*={{16,-120,240,-140},{-120,1200,-2700,1680},{240,-2700,6480,4200},{-140,1680,-4200,2800}}
  (n=16)
  */

  f[0] = 0;
  int i, j, k, b;

  long double sum = 0;

  static long double hilbert[10][10],
      y[10][10]; // Increase matrix size if D > 100

  b = (int)sqrt((double)D);

  for (i = 0; i < b; i++) {
    for (j = 0; j < b; j++) {
      hilbert[i][j] =
          1. / (double)(i + j + 1); // Create a static Hilbert matrix
    }
  }

  for (j = 0; j < b; j++) {
    for (k = 0; k < b; k++) {
      y[j][k] = 0;
      for (i = 0; i < b; i++) {
        y[j][k] += hilbert[j][i] * x[k + b * i]; // Compute matrix product H*x
      }
    }
  }

  for (i = 0; i < b; i++) {
    for (j = 0; j < b; j++) {
      if (i == j)
        sum += fabs(y[i][j] - 1); // Sum absolute value of deviations
      else
        sum += fabs(y[i][j]);
    }
  }

  f[0] += sum;
}

void cec2019_Chebyshev(
    double *x, int D,
    double *f) // Storn's Tchebychev - a 2nd ICEO function - generalized version
{
  /* Valid for any D>2
  constraints: unconstrained
  type: multi-modal with one global minimum; non-separable
  initial upper bound = 2^D, initial lower bound = -D^n
  value-to-reach = f(x*)+1.0e-8
  f(x*)=0.0; x*=(128,0,-256,0,160,0,-32,0,1) (n=9)
  x*=(32768,0,-131072,0,212992,0,-180224,0,84480,0,-21504,0,2688,0,-128,0,1)
  (n=17)
  */

  f[0] = 0.0;
  int i, j;
  static int sample;
  long double a = 1., b = 1.2, px, y = -1, sum = 0;
  static long double dx, dy;

  for (j = 0; j < D - 2; j++) {
    dx = 2.4 * b - a;
    a = b;
    b = dx;
  }

  sample = 32 * D;
  dy = 2. / (long double)sample;

  for (i = 0; i <= sample; i++) {
    px = x[0];
    for (j = 1; j < D; j++) {
      px = y * px + x[j];
    }
    if (px < -1 || px > 1)
      sum += (1. - fabs(px)) * (1. - fabs(px));
    y += dy;
  }

  for (i = -1; i <= 1; i += 2) {
    px = x[0];
    for (j = 1; j < D; j++) {
      px = 1.2 * px + x[j];
    }

    if (px < dx)
      sum += px * px;
  }

  f[0] += sum;
}

#endif // CEC2019_FUNCTIONS_H_
