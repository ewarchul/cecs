/*
  CEC13 Test Function Suite for Single Objective Bound Constrained Numerical
  Optimization
*/

#include "cec2013_functions.h"

void cec2013_sphere_func(double *x, double *f, int nx, double *Os, double *Mr,
                         int r_flag) /* Sphere */
{
  int i;
  cec2013_shiftfunc(x, y, nx, Os);
  if (r_flag == 1)
    cec2013_rotatefunc(y, z, nx, Mr);
  else
    for (i = 0; i < nx; i++)
      z[i] = y[i];
  f[0] = 0.0;
  for (i = 0; i < nx; i++) {
    f[0] += z[i] * z[i];
  }
}

void cec2013_ellips_func(double *x, double *f, int nx, double *Os, double *Mr,
                         int r_flag) /* Ellipsoidal */
{
  int i;
  cec2013_shiftfunc(x, y, nx, Os);
  if (r_flag == 1)
    cec2013_rotatefunc(y, z, nx, Mr);
  else
    for (i = 0; i < nx; i++)
      z[i] = y[i];
  cec2013_oszfunc(z, y, nx);
  f[0] = 0.0;
  for (i = 0; i < nx; i++) {
    f[0] += pow(10.0, 6.0 * i / (nx - 1)) * y[i] * y[i];
  }
}

void cec2013_bent_cigar_func(double *x, double *f, int nx, double *Os,
                             double *Mr, int r_flag) /* Bent_Cigar */
{
  int i;
  double beta = 0.5;
  cec2013_shiftfunc(x, y, nx, Os);
  if (r_flag == 1)
    cec2013_rotatefunc(y, z, nx, Mr);
  else
    for (i = 0; i < nx; i++)
      z[i] = y[i];
  cec2013_asyfunc(z, y, nx, beta);
  if (r_flag == 1)
    cec2013_rotatefunc(y, z, nx, &Mr[nx * nx]);
  else
    for (i = 0; i < nx; i++)
      z[i] = y[i];

  f[0] = z[0] * z[0];
  for (i = 1; i < nx; i++) {
    f[0] += pow(10.0, 6.0) * z[i] * z[i];
  }
}

void cec2013_discus_func(double *x, double *f, int nx, double *Os, double *Mr,
                         int r_flag) /* Discus */
{
  int i;
  cec2013_shiftfunc(x, y, nx, Os);
  if (r_flag == 1)
    cec2013_rotatefunc(y, z, nx, Mr);
  else
    for (i = 0; i < nx; i++)
      z[i] = y[i];
  cec2013_oszfunc(z, y, nx);

  f[0] = pow(10.0, 6.0) * y[0] * y[0];
  for (i = 1; i < nx; i++) {
    f[0] += y[i] * y[i];
  }
}

void cec2013_dif_powers_func(double *x, double *f, int nx, double *Os,
                             double *Mr, int r_flag) /* Different Powers */
{
  int i;
  cec2013_shiftfunc(x, y, nx, Os);
  if (r_flag == 1)
    cec2013_rotatefunc(y, z, nx, Mr);
  else
    for (i = 0; i < nx; i++)
      z[i] = y[i];
  f[0] = 0.0;
  for (i = 0; i < nx; i++) {
    f[0] += pow(fabs(z[i]), 2 + 4 * i / (nx - 1));
  }
  f[0] = pow(f[0], 0.5);
}

void cec2013_rosenbrock_func(double *x, double *f, int nx, double *Os,
                             double *Mr, int r_flag) /* Rosenbrock's */
{
  int i;
  double tmp1, tmp2;
  cec2013_shiftfunc(x, y, nx, Os); // shift
  for (i = 0; i < nx; i++)         // shrink to the orginal search range
  {
    y[i] = y[i] * 2.048 / 100;
  }
  if (r_flag == 1)
    cec2013_rotatefunc(y, z, nx, Mr); // rotate
  else
    for (i = 0; i < nx; i++)
      z[i] = y[i];
  for (i = 0; i < nx; i++) // shift to orgin
  {
    z[i] = z[i] + 1;
  }

  f[0] = 0.0;
  for (i = 0; i < nx - 1; i++) {
    tmp1 = z[i] * z[i] - z[i + 1];
    tmp2 = z[i] - 1.0;
    f[0] += 100.0 * tmp1 * tmp1 + tmp2 * tmp2;
  }
}

void cec2013_schaffer_F7_func(double *x, double *f, int nx, double *Os,
                              double *Mr, int r_flag) /* Schwefel's 1.2  */
{
  int i;
  double tmp;
  cec2013_shiftfunc(x, y, nx, Os);
  if (r_flag == 1)
    cec2013_rotatefunc(y, z, nx, Mr);
  else
    for (i = 0; i < nx; i++)
      z[i] = y[i];
  cec2013_asyfunc(z, y, nx, 0.5);
  for (i = 0; i < nx; i++)
    z[i] = y[i] * pow(10.0, 1.0 * i / (nx - 1) / 2.0);
  if (r_flag == 1)
    cec2013_rotatefunc(z, y, nx, &Mr[nx * nx]);
  else
    for (i = 0; i < nx; i++)
      y[i] = z[i];

  for (i = 0; i < nx - 1; i++)
    z[i] = pow(y[i] * y[i] + y[i + 1] * y[i + 1], 0.5);
  f[0] = 0.0;
  for (i = 0; i < nx - 1; i++) {
    tmp = sin(50.0 * pow(z[i], 0.2));
    f[0] += pow(z[i], 0.5) + pow(z[i], 0.5) * tmp * tmp;
  }
  f[0] = f[0] * f[0] / (nx - 1) / (nx - 1);
}

void cec2013_ackley_func(double *x, double *f, int nx, double *Os, double *Mr,
                         int r_flag) /* Ackley's  */
{
  int i;
  double sum1, sum2;

  cec2013_shiftfunc(x, y, nx, Os);
  if (r_flag == 1)
    cec2013_rotatefunc(y, z, nx, Mr);
  else
    for (i = 0; i < nx; i++)
      z[i] = y[i];

  cec2013_asyfunc(z, y, nx, 0.5);
  for (i = 0; i < nx; i++)
    z[i] = y[i] * pow(10.0, 1.0 * i / (nx - 1) / 2.0);
  if (r_flag == 1)
    cec2013_rotatefunc(z, y, nx, &Mr[nx * nx]);
  else
    for (i = 0; i < nx; i++)
      y[i] = z[i];

  sum1 = 0.0;
  sum2 = 0.0;
  for (i = 0; i < nx; i++) {
    sum1 += y[i] * y[i];
    sum2 += cos(2.0 * M_PI * y[i]);
  }
  sum1 = -0.2 * sqrt(sum1 / nx);
  sum2 /= nx;
  f[0] = E - 20.0 * exp(sum1) - exp(sum2) + 20.0;
}

void cec2013_weierstrass_func(double *x, double *f, int nx, double *Os,
                              double *Mr, int r_flag) /* Weierstrass's  */
{
  int i, j, k_max;
  double sum, sum2 = 0.0, a, b;

  cec2013_shiftfunc(x, y, nx, Os);
  for (i = 0; i < nx; i++) // shrink to the orginal search range
  {
    y[i] = y[i] * 0.5 / 100;
  }
  if (r_flag == 1)
    cec2013_rotatefunc(y, z, nx, Mr);
  else
    for (i = 0; i < nx; i++)
      z[i] = y[i];

  cec2013_asyfunc(z, y, nx, 0.5);
  for (i = 0; i < nx; i++)
    z[i] = y[i] * pow(10.0, 1.0 * i / (nx - 1) / 2.0);
  if (r_flag == 1)
    cec2013_rotatefunc(z, y, nx, &Mr[nx * nx]);
  else
    for (i = 0; i < nx; i++)
      y[i] = z[i];

  a = 0.5;
  b = 3.0;
  k_max = 20;
  f[0] = 0.0;
  for (i = 0; i < nx; i++) {
    sum = 0.0;
    sum2 = 0.0;
    for (j = 0; j <= k_max; j++) {
      sum += pow(a, j) * cos(2.0 * M_PI * pow(b, j) * (y[i] + 0.5));
      sum2 += pow(a, j) * cos(2.0 * M_PI * pow(b, j) * 0.5);
    }
    f[0] += sum;
  }
  f[0] -= nx * sum2;
}

void cec2013_griewank_func(double *x, double *f, int nx, double *Os, double *Mr,
                   int r_flag) /* Griewank's  */
{
  int i;
  double s, p;

  cec2013_shiftfunc(x, y, nx, Os);
  for (i = 0; i < nx; i++) // shrink to the orginal search range
  {
    y[i] = y[i] * 600.0 / 100.0;
  }
  if (r_flag == 1)
    cec2013_rotatefunc(y, z, nx, Mr);
  else
    for (i = 0; i < nx; i++)
      z[i] = y[i];

  for (i = 0; i < nx; i++)
    z[i] = z[i] * pow(100.0, 1.0 * i / (nx - 1) / 2.0);

  s = 0.0;
  p = 1.0;
  for (i = 0; i < nx; i++) {
    s += z[i] * z[i];
    p *= cos(z[i] / sqrt(1.0 + i));
  }
  f[0] = 1.0 + s / 4000.0 - p;
}

void cec2013_rastrigin_func(double *x, double *f, int nx, double *Os,
                            double *Mr, int r_flag) /* Rastrigin's  */
{
  int i;
  double alpha = 10.0, beta = 0.2;
  cec2013_shiftfunc(x, y, nx, Os);
  for (i = 0; i < nx; i++) // shrink to the orginal search range
  {
    y[i] = y[i] * 5.12 / 100;
  }

  if (r_flag == 1)
    cec2013_rotatefunc(y, z, nx, Mr);
  else
    for (i = 0; i < nx; i++)
      z[i] = y[i];

  cec2013_oszfunc(z, y, nx);
  cec2013_asyfunc(y, z, nx, beta);

  if (r_flag == 1)
    cec2013_rotatefunc(z, y, nx, &Mr[nx * nx]);
  else
    for (i = 0; i < nx; i++)
      y[i] = z[i];

  for (i = 0; i < nx; i++) {
    y[i] *= pow(alpha, 1.0 * i / (nx - 1) / 2);
  }

  if (r_flag == 1)
    cec2013_rotatefunc(y, z, nx, Mr);
  else
    for (i = 0; i < nx; i++)
      z[i] = y[i];

  f[0] = 0.0;
  for (i = 0; i < nx; i++) {
    f[0] += (z[i] * z[i] - 10.0 * cos(2.0 * M_PI * z[i]) + 10.0);
  }
}

void cec2013_step_rastrigin_func(double *x, double *f, int nx, double *Os,
                                 double *Mr,
                                 int r_flag) /* Noncontinuous Rastrigin's  */
{
  int i;
  double alpha = 10.0, beta = 0.2;
  cec2013_shiftfunc(x, y, nx, Os);
  for (i = 0; i < nx; i++) // shrink to the orginal search range
  {
    y[i] = y[i] * 5.12 / 100;
  }

  if (r_flag == 1)
    cec2013_rotatefunc(y, z, nx, Mr);
  else
    for (i = 0; i < nx; i++)
      z[i] = y[i];

  for (i = 0; i < nx; i++) {
    if (fabs(z[i]) > 0.5)
      z[i] = floor(2 * z[i] + 0.5) / 2;
  }

  cec2013_oszfunc(z, y, nx);
  cec2013_asyfunc(y, z, nx, beta);

  if (r_flag == 1)
    cec2013_rotatefunc(z, y, nx, &Mr[nx * nx]);
  else
    for (i = 0; i < nx; i++)
      y[i] = z[i];

  for (i = 0; i < nx; i++) {
    y[i] *= pow(alpha, 1.0 * i / (nx - 1) / 2);
  }

  if (r_flag == 1)
    cec2013_rotatefunc(y, z, nx, Mr);
  else
    for (i = 0; i < nx; i++)
      z[i] = y[i];

  f[0] = 0.0;
  for (i = 0; i < nx; i++) {
    f[0] += (z[i] * z[i] - 10.0 * cos(2.0 * M_PI * z[i]) + 10.0);
  }
}

void cec2013_schwefel_func(double *x, double *f, int nx, double *Os, double *Mr,
                           int r_flag) /* Schwefel's  */
{
  int i;
  double tmp;
  cec2013_shiftfunc(x, y, nx, Os);
  for (i = 0; i < nx; i++) // shrink to the orginal search range
  {
    y[i] *= 1000 / 100;
  }
  if (r_flag == 1)
    cec2013_rotatefunc(y, z, nx, Mr);
  else
    for (i = 0; i < nx; i++)
      z[i] = y[i];

  for (i = 0; i < nx; i++)
    y[i] = z[i] * pow(10.0, 1.0 * i / (nx - 1) / 2.0);

  for (i = 0; i < nx; i++)
    z[i] = y[i] + 4.209687462275036e+002;

  f[0] = 0;
  for (i = 0; i < nx; i++) {
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
  f[0] = 4.189828872724338e+002 * nx + f[0];
}

void cec2013_katsuura_func(double *x, double *f, int nx, double *Os, double *Mr,
                           int r_flag) /* Katsuura  */
{
  int i, j;
  double temp, tmp1, tmp2, tmp3;
  tmp3 = pow(1.0 * nx, 1.2);
  cec2013_shiftfunc(x, y, nx, Os);
  for (i = 0; i < nx; i++) // shrink to the orginal search range
  {
    y[i] *= 5.0 / 100.0;
  }
  if (r_flag == 1)
    cec2013_rotatefunc(y, z, nx, Mr);
  else
    for (i = 0; i < nx; i++)
      z[i] = y[i];

  for (i = 0; i < nx; i++)
    z[i] *= pow(100.0, 1.0 * i / (nx - 1) / 2.0);

  if (r_flag == 1)
    cec2013_rotatefunc(z, y, nx, &Mr[nx * nx]);
  else
    for (i = 0; i < nx; i++)
      y[i] = z[i];

  f[0] = 1.0;
  for (i = 0; i < nx; i++) {
    temp = 0.0;
    for (j = 1; j <= 32; j++) {
      tmp1 = pow(2.0, j);
      tmp2 = tmp1 * y[i];
      temp += fabs(tmp2 - floor(tmp2 + 0.5)) / tmp1;
    }
    f[0] *= pow(1.0 + (i + 1) * temp, 10.0 / tmp3);
  }
  tmp1 = 10.0 / nx / nx;
  f[0] = f[0] * tmp1 - tmp1;
}

void cec2013_bi_rastrigin_func(double *x, double *f, int nx, double *Os,
                               double *Mr,
                               int r_flag) /* Lunacek Bi_rastrigin Function */
{
  int i;
  double mu0 = 2.5, d = 1.0, s, mu1, tmp, tmp1, tmp2;
  double *tmpx;
  tmpx = (double *)malloc(sizeof(double) * nx);
  s = 1.0 - 1.0 / (2.0 * pow(nx + 20.0, 0.5) - 8.2);
  mu1 = -pow((mu0 * mu0 - d) / s, 0.5);

  cec2013_shiftfunc(x, y, nx, Os);
  for (i = 0; i < nx; i++) // shrink to the orginal search range
  {
    y[i] *= 10.0 / 100.0;
  }

  for (i = 0; i < nx; i++) {
    tmpx[i] = 2 * y[i];
    if (Os[i] < 0.)
      tmpx[i] *= -1.;
  }

  for (i = 0; i < nx; i++) {
    z[i] = tmpx[i];
    tmpx[i] += mu0;
  }
  if (r_flag == 1)
    cec2013_rotatefunc(z, y, nx, Mr);
  else
    for (i = 0; i < nx; i++)
      y[i] = z[i];

  for (i = 0; i < nx; i++)
    y[i] *= pow(100.0, 1.0 * i / (nx - 1) / 2.0);
  if (r_flag == 1)
    cec2013_rotatefunc(y, z, nx, &Mr[nx * nx]);
  else
    for (i = 0; i < nx; i++)
      z[i] = y[i];

  tmp1 = 0.0;
  tmp2 = 0.0;
  for (i = 0; i < nx; i++) {
    tmp = tmpx[i] - mu0;
    tmp1 += tmp * tmp;
    tmp = tmpx[i] - mu1;
    tmp2 += tmp * tmp;
  }
  tmp2 *= s;
  tmp2 += d * nx;
  tmp = 0;
  for (i = 0; i < nx; i++) {
    tmp += cos(2.0 * M_PI * z[i]);
  }

  if (tmp1 < tmp2)
    f[0] = tmp1;
  else
    f[0] = tmp2;
  f[0] += 10.0 * (nx - tmp);
  free(tmpx);
}

void cec2013_grie_rosen_func(double *x, double *f, int nx, double *Os,
                             double *Mr, int r_flag) /* Griewank-Rosenbrock  */
{
  int i;
  double temp, tmp1, tmp2;

  cec2013_shiftfunc(x, y, nx, Os);
  for (i = 0; i < nx; i++) // shrink to the orginal search range
  {
    y[i] = y[i] * 5 / 100;
  }
  if (r_flag == 1)
    cec2013_rotatefunc(y, z, nx, Mr);
  else
    for (i = 0; i < nx; i++)
      z[i] = y[i];

  for (i = 0; i < nx; i++) // shift to orgin
  {
    z[i] = y[i] + 1;
  }

  f[0] = 0.0;
  for (i = 0; i < nx - 1; i++) {
    tmp1 = z[i] * z[i] - z[i + 1];
    tmp2 = z[i] - 1.0;
    temp = 100.0 * tmp1 * tmp1 + tmp2 * tmp2;
    f[0] += (temp * temp) / 4000.0 - cos(temp) + 1.0;
  }
  tmp1 = z[nx - 1] * z[nx - 1] - z[0];
  tmp2 = z[nx - 1] - 1.0;
  temp = 100.0 * tmp1 * tmp1 + tmp2 * tmp2;
  ;
  f[0] += (temp * temp) / 4000.0 - cos(temp) + 1.0;
}

void cec2013_escaffer6_func(double *x, double *f, int nx, double *Os,
                            double *Mr, int r_flag) /* Expanded Scaffer's F6  */
{
  int i;
  double temp1, temp2;
  cec2013_shiftfunc(x, y, nx, Os);
  if (r_flag == 1)
    cec2013_rotatefunc(y, z, nx, Mr);
  else
    for (i = 0; i < nx; i++)
      z[i] = y[i];

  cec2013_asyfunc(z, y, nx, 0.5);
  if (r_flag == 1)
    cec2013_rotatefunc(y, z, nx, &Mr[nx * nx]);
  else
    for (i = 0; i < nx; i++)
      z[i] = y[i];

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

void cec2013_cf01(double *x, double *f, int nx, double *Os, double *Mr,
                  int r_flag) /* Composition Function 1 */
{
  int i, cf_num = 5;
  double fit[5];
  double delta[5] = {10, 20, 30, 40, 50};
  double bias[5] = {0, 100, 200, 300, 400};

  i = 0;
  cec2013_rosenbrock_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx],
                          r_flag);
  fit[i] = 10000 * fit[i] / 1e+4;
  i = 1;
  cec2013_dif_powers_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx],
                          r_flag);
  fit[i] = 10000 * fit[i] / 1e+10;
  i = 2;
  cec2013_bent_cigar_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx],
                          r_flag);
  fit[i] = 10000 * fit[i] / 1e+30;
  i = 3;
  cec2013_discus_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], r_flag);
  fit[i] = 10000 * fit[i] / 1e+10;
  i = 4;
  cec2013_sphere_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 0);
  fit[i] = 10000 * fit[i] / 1e+5;
  cec2013_cf_cal(x, f, nx, Os, delta, bias, fit, cf_num);
}

void cec2013_cf02(double *x, double *f, int nx, double *Os, double *Mr,
                  int r_flag) /* Composition Function 2 */
{
  int i, cf_num = 3;
  double fit[3];
  double delta[3] = {20, 20, 20};
  double bias[3] = {0, 100, 200};
  for (i = 0; i < cf_num; i++) {
    cec2013_schwefel_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx],
                          r_flag);
  }
  cec2013_cf_cal(x, f, nx, Os, delta, bias, fit, cf_num);
}

void cec2013_cf03(double *x, double *f, int nx, double *Os, double *Mr,
                  int r_flag) /* Composition Function 3 */
{
  int i, cf_num = 3;
  double fit[3];
  double delta[3] = {20, 20, 20};
  double bias[3] = {0, 100, 200};
  for (i = 0; i < cf_num; i++) {
    cec2013_schwefel_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx],
                          r_flag);
  }
  cec2013_cf_cal(x, f, nx, Os, delta, bias, fit, cf_num);
}

void cec2013_cf04(double *x, double *f, int nx, double *Os, double *Mr,
                  int r_flag) /* Composition Function 4 */
{
  int i, cf_num = 3;
  double fit[3];
  double delta[3] = {20, 20, 20};
  double bias[3] = {0, 100, 200};
  i = 0;
  cec2013_schwefel_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], r_flag);
  fit[i] = 1000 * fit[i] / 4e+3;
  i = 1;
  cec2013_rastrigin_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], r_flag);
  fit[i] = 1000 * fit[i] / 1e+3;
  i = 2;
  cec2013_weierstrass_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx],
                           r_flag);
  fit[i] = 1000 * fit[i] / 400;
  cec2013_cf_cal(x, f, nx, Os, delta, bias, fit, cf_num);
}

void cec2013_cf05(double *x, double *f, int nx, double *Os, double *Mr,
                  int r_flag) /* Composition Function 4 */
{
  int i, cf_num = 3;
  double fit[3];
  double delta[3] = {10, 30, 50};
  double bias[3] = {0, 100, 200};
  i = 0;
  cec2013_schwefel_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], r_flag);
  fit[i] = 1000 * fit[i] / 4e+3;
  i = 1;
  cec2013_rastrigin_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], r_flag);
  fit[i] = 1000 * fit[i] / 1e+3;
  i = 2;
  cec2013_weierstrass_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx],
                           r_flag);
  fit[i] = 1000 * fit[i] / 400;
  cec2013_cf_cal(x, f, nx, Os, delta, bias, fit, cf_num);
}

void cec2013_cf06(double *x, double *f, int nx, double *Os, double *Mr,
                  int r_flag) /* Composition Function 6 */
{
  int i, cf_num = 5;
  double fit[5];
  double delta[5] = {10, 10, 10, 10, 10};
  double bias[5] = {0, 100, 200, 300, 400};
  i = 0;
  cec2013_schwefel_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], r_flag);
  fit[i] = 1000 * fit[i] / 4e+3;
  i = 1;
  cec2013_rastrigin_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], r_flag);
  fit[i] = 1000 * fit[i] / 1e+3;
  i = 2;
  cec2013_ellips_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], r_flag);
  fit[i] = 1000 * fit[i] / 1e+10;
  i = 3;
  cec2013_weierstrass_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx],
                           r_flag);
  fit[i] = 1000 * fit[i] / 400;
  i = 4;
  cec2013_griewank_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], r_flag);
  fit[i] = 1000 * fit[i] / 100;
  cec2013_cf_cal(x, f, nx, Os, delta, bias, fit, cf_num);
}

void cec2013_cf07(double *x, double *f, int nx, double *Os, double *Mr,
                  int r_flag) /* Composition Function 7 */
{
  int i, cf_num = 5;
  double fit[5];
  double delta[5] = {10, 10, 10, 20, 20};
  double bias[5] = {0, 100, 200, 300, 400};
  i = 0;
  cec2013_griewank_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], r_flag);
  fit[i] = 10000 * fit[i] / 100;
  i = 1;
  cec2013_rastrigin_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], r_flag);
  fit[i] = 10000 * fit[i] / 1e+3;
  i = 2;
  cec2013_schwefel_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], r_flag);
  fit[i] = 10000 * fit[i] / 4e+3;
  i = 3;
  cec2013_weierstrass_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx],
                           r_flag);
  fit[i] = 10000 * fit[i] / 400;
  i = 4;
  cec2013_sphere_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 0);
  fit[i] = 10000 * fit[i] / 1e+5;
  cec2013_cf_cal(x, f, nx, Os, delta, bias, fit, cf_num);
}

void cec2013_cf08(double *x, double *f, int nx, double *Os, double *Mr,
                  int r_flag) /* Composition Function 8 */
{
  int i, cf_num = 5;
  double fit[5];
  double delta[5] = {10, 20, 30, 40, 50};
  double bias[5] = {0, 100, 200, 300, 400};
  i = 0;
  cec2013_grie_rosen_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx],
                          r_flag);
  fit[i] = 10000 * fit[i] / 4e+3;
  i = 1;
  cec2013_schaffer_F7_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx],
                           r_flag);
  fit[i] = 10000 * fit[i] / 4e+6;
  i = 2;
  cec2013_schwefel_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], r_flag);
  fit[i] = 10000 * fit[i] / 4e+3;
  i = 3;
  cec2013_escaffer6_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], r_flag);
  fit[i] = 10000 * fit[i] / 2e+7;
  i = 4;
  cec2013_sphere_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 0);
  fit[i] = 10000 * fit[i] / 1e+5;
  cec2013_cf_cal(x, f, nx, Os, delta, bias, fit, cf_num);
}

void cec2013_shiftfunc(double *x, double *xshift, int nx, double *Os) {
  int i;
  for (i = 0; i < nx; i++) {
    xshift[i] = x[i] - Os[i];
  }
}

void cec2013_rotatefunc(double *x, double *xrot, int nx, double *Mr) {
  int i, j;
  for (i = 0; i < nx; i++) {
    xrot[i] = 0;
    for (j = 0; j < nx; j++) {
      xrot[i] = xrot[i] + x[j] * Mr[i * nx + j];
    }
  }
}

void cec2013_asyfunc(double *x, double *xasy, int nx, double beta) {
  int i;
  for (i = 0; i < nx; i++) {
    if (x[i] > 0)
      xasy[i] = pow(x[i], 1.0 + beta * i / (nx - 1) * pow(x[i], 0.5));
  }
}

void cec2013_oszfunc(double *x, double *xosz, int nx) {
  int i, sx;
  double c1, c2, xx = 0;
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

void cec2013_cf_cal(double *x, double *f, int nx, double *Os, double *delta,
                    double *bias, double *fit, int cf_num) {
  int i, j;
  double *w;
  double w_max = 0, w_sum = 0;
  w = (double *)malloc(cf_num * sizeof(double));
  for (i = 0; i < cf_num; i++) {
    fit[i] += bias[i];
    w[i] = 0;
    for (j = 0; j < nx; j++) {
      w[i] += pow(x[j] - Os[i * nx + j], 2.0);
    }
    if (w[i] != 0)
      w[i] = pow(1.0 / w[i], 0.5) * exp(-w[i] / 2.0 / nx / pow(delta[i], 2.0));
    else
      w[i] = INF;
    if (w[i] > w_max)
      w_max = w[i];
  }

  for (i = 0; i < cf_num; i++) {
    w_sum = w_sum + w[i];
  }
  if (w_max == 0) {
    for (i = 0; i < cf_num; i++)
      w[i] = 1;
    w_sum = cf_num;
  }
  f[0] = 0.0;
  for (i = 0; i < cf_num; i++) {
    f[0] = f[0] + w[i] / w_sum * fit[i];
  }
  free(w);
}
