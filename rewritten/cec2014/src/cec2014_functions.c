#include "../include/cec2014_functions.h"

void cec2014_bent_cigar_func(double *x, double *f, int nx, double *Os,
                             double *Mr, int s_flag,
                             int r_flag) /* Bent_Cigar */
{
  int i;
  double *y = malloc(nx * sizeof(double));
  double *z = malloc(nx * sizeof(double));
  cec2014_sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag,
                  y); /* shift and rotate */

  f[0] = z[0] * z[0];
  for (i = 1; i < nx; i++) {
    f[0] += pow(10.0, 6.0) * z[i] * z[i];
  }
  free(y);
  free(z);
}

void cec2014_shiftfunc(double *x, double *xshift, int nx, double *Os) {
  int i;
  for (i = 0; i < nx; i++) {
    xshift[i] = x[i] - Os[i];
  }
}

void cec2014_rotatefunc(double *x, double *xrot, int nx, double *Mr) {
  int i, j;
  for (i = 0; i < nx; i++) {
    xrot[i] = 0;
    for (j = 0; j < nx; j++) {
      xrot[i] = xrot[i] + x[j] * Mr[i * nx + j];
    }
  }
}

void cec2014_sr_func(double *x, double *sr_x, int nx, double *Os, double *Mr,
                     double sh_rate, int s_flag, int r_flag, double *y) {
  int i;
  if (s_flag == 1) {
    if (r_flag == 1) {
      cec2014_shiftfunc(x, y, nx, Os);
      for (i = 0; i < nx; i++) // shrink to the original search range
      {
        y[i] = y[i] * sh_rate;
      }
      cec2014_rotatefunc(y, sr_x, nx, Mr);
    } else {
      cec2014_shiftfunc(x, sr_x, nx, Os);
      for (i = 0; i < nx; i++) // shrink to the original search range
      {
        sr_x[i] = sr_x[i] * sh_rate;
      }
    }
  } else {

    if (r_flag == 1) {
      for (i = 0; i < nx; i++) // shrink to the original search range
      {
        y[i] = x[i] * sh_rate;
      }
      cec2014_rotatefunc(y, sr_x, nx, Mr);
    } else
      for (i = 0; i < nx; i++) // shrink to the original search range
      {
        sr_x[i] = x[i] * sh_rate;
      }
  }
}

void cec2014_sphere_func(double *x, double *f, int nx, double *Os, double *Mr,
                         int s_flag, int r_flag) {
  int i;
  f[0] = 0.0;
  double *y = malloc(nx * sizeof(double));
  double *z = malloc(nx * sizeof(double));
  cec2014_sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag, y);
  for (i = 0; i < nx; i++) {
    f[0] += z[i] * z[i];
  }
  free(y);
  free(z);
}

void cec2014_ellips_func(double *x, double *f, int nx, double *Os, double *Mr,
                         int s_flag, int r_flag) {
  int i;
  f[0] = 0.0;
  double *y = malloc(nx * sizeof(double));
  double *z = malloc(nx * sizeof(double));
  cec2014_sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag, y);
  for (i = 0; i < nx; i++) {
    f[0] += pow(10.0, 6.0 * i / (nx - 1)) * z[i] * z[i];
  }
  free(y);
  free(z);
}

void cec2014_sum_diff_pow_func(double *x, double *f, int nx, double *Os,
                               double *Mr, int s_flag, int r_flag) {
  int i;
  double *y = malloc(nx * sizeof(double));
  double *z = malloc(nx * sizeof(double));
  cec2014_sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag, y);
  f[0] = 0.0;
  double sum = 0.0;
  for (i = 0; i < nx; i++) {
    double xi = z[i];
    double newv = pow((fabs(xi)), (i + 1));
    sum = sum + newv;
  }
  f[0] = sum;
  free(y);
  free(z);
}

void cec2014_zakharov_func(double *x, double *f, int nx, double *Os, double *Mr,
                           int s_flag, int r_flag) {
  int i;
  double *y = malloc(nx * sizeof(double));
  double *z = malloc(nx * sizeof(double));
  cec2014_sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag, y);
  f[0] = 0.0;
  double sum1 = 0.0;
  double sum2 = 0.0;
  for (i = 0; i < nx; i++) {
    double xi = z[i];
    sum1 = sum1 + pow(xi, 2);
    sum2 = sum2 + 0.5 * (i + 1) * xi;
  }
  f[0] = sum1 + pow(sum2, 2) + pow(sum2, 4);
  free(y);
  free(z);
}

void cec2014_levy_func(double *x, double *f, int nx, double *Os, double *Mr,
                       int s_flag, int r_flag) {
  int i;
  f[0] = 0.0;
  double *y = malloc(nx * sizeof(double));
  double *z = malloc(nx * sizeof(double));
  cec2014_sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag, y);
  double *w;
  w = (double *)malloc(sizeof(double) * nx);
  for (i = 0; i < nx; i++) {
    w[i] = 1.0 + (z[i] - 1.0) / 4.0;
  }
  double term1 = pow((sin(M_PI * w[0])), 2);
  double term3 =
      pow((w[nx - 1] - 1), 2) * (1 + pow((sin(2 * M_PI * w[nx - 1])), 2));
  double sum = 0.0;
  for (i = 0; i < nx - 1; i++) {
    double wi = w[i];
    double newv = pow((wi - 1), 2) * (1 + 10 * pow((sin(M_PI * wi + 1)), 2));
    sum = sum + newv;
  }
  f[0] = term1 + sum + term3;
  free(w);
  free(y);
  free(z);
}

void cec2014_dixon_price_func(double *x, double *f, int nx, double *Os,
                              double *Mr, int s_flag, int r_flag) {
  int i;
  double *y = malloc(nx * sizeof(double));
  double *z = malloc(nx * sizeof(double));
  cec2014_sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag, y);
  f[0] = 0;
  double x1 = z[0];
  ;
  double term1 = pow((x1 - 1), 2);

  double sum = 0;
  for (i = 1; i < nx; i++) {
    double xi = z[i];
    double xold = z[i - 1];
    double newv = i * pow((pow(2 * xi, 2) - xold), 2);
    sum = sum + newv;
  }
  f[0] = term1 + sum;
  free(y);
  free(z);
}

void cec2014_discus_func(double *x, double *f, int nx, double *Os, double *Mr,
                         int s_flag, int r_flag) {
  int i;
  double *y = malloc(nx * sizeof(double));
  double *z = malloc(nx * sizeof(double));
  cec2014_sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag, y);
  f[0] = pow(10.0, 6.0) * z[0] * z[0];
  for (i = 1; i < nx; i++) {
    f[0] += z[i] * z[i];
  }
  free(y);
  free(z);
}

void cec2014_dif_powers_func(double *x, double *f, int nx, double *Os,
                             double *Mr, int s_flag, int r_flag) {
  int i;
  f[0] = 0.0;
  double *y = malloc(nx * sizeof(double));
  double *z = malloc(nx * sizeof(double));
  cec2014_sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag, y);
  for (i = 0; i < nx; i++) {
    f[0] += pow(fabs(z[i]), 2 + 4 * i / (nx - 1));
  }
  f[0] = pow(f[0], 0.5);
  free(y);
  free(z);
}

void cec2014_rosenbrock_func(double *x, double *f, int nx, double *Os,
                             double *Mr, int s_flag, int r_flag) {
  int i;
  double tmp1, tmp2;
  f[0] = 0.0;
  double *y = malloc(nx * sizeof(double));
  double *z = malloc(nx * sizeof(double));
  cec2014_sr_func(x, z, nx, Os, Mr, 2.048 / 100.0, s_flag, r_flag, y);
  z[0] += 1.0; // shift to orgin
  for (i = 0; i < nx - 1; i++) {
    z[i + 1] += 1.0; // shift to orgin
    tmp1 = z[i] * z[i] - z[i + 1];
    tmp2 = z[i] - 1.0;
    f[0] += 100.0 * tmp1 * tmp1 + tmp2 * tmp2;
  }
  free(y);
  free(z);
}

void cec2014_schaffer_F7_func(double *x, double *f, int nx, double *Os,
                              double *Mr, int s_flag, int r_flag) {
  int i;
  double tmp;
  f[0] = 0.0;
  double *y = malloc(nx * sizeof(double));
  double *z = malloc(nx * sizeof(double));
  cec2014_sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag, y); /* shift and rotate
                                                              */
  for (i = 0; i < nx - 1; i++) {
    z[i] = pow(y[i] * y[i] + y[i + 1] * y[i + 1], 0.5);
    tmp = sin(50.0 * pow(z[i], 0.2));
    f[0] += pow(z[i], 0.5) + pow(z[i], 0.5) * tmp * tmp;
  }
  f[0] = f[0] * f[0] / (nx - 1) / (nx - 1);
  free(y);
  free(z);
}

void cec2014_ackley_func(double *x, double *f, int nx, double *Os, double *Mr,
                         int s_flag, int r_flag) {
  int i;
  double sum1 = 0.0;
  double sum2 = 0.0;
  double *y = malloc(nx * sizeof(double));
  double *z = malloc(nx * sizeof(double));

  cec2014_sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag, y); /* shift and rotate
                                                              */

  for (i = 0; i < nx; i++) {
    sum1 += z[i] * z[i];
    sum2 += cos(2.0 * M_PI * z[i]);
  }
  sum1 = -0.2 * sqrt(sum1 / nx);
  sum2 /= nx;
  f[0] = E - 20.0 * exp(sum1) - exp(sum2) + 20.0;
  free(y);
  free(z);
}

void cec2014_weierstrass_func(double *x, double *f, int nx, double *Os,
                              double *Mr, int s_flag, int r_flag) {
  int i, j, k_max;
  double sum = 0.0;
  double sum2 = 0.0;
  double a = 0.5;
  double b = 3.0;
  k_max = 20;
  f[0] = 0.0;
  double *y = malloc(nx * sizeof(double));
  double *z = malloc(nx * sizeof(double));
  cec2014_sr_func(x, z, nx, Os, Mr, 0.5 / 100.0, s_flag, r_flag, y);
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
  free(y);
  free(z);
}

void cec2014_griewank_func(double *x, double *f, int nx, double *Os, double *Mr,
                           int s_flag, int r_flag) {
  int i;
  double s, p;
  s = 0.0;
  p = 1.0;
  double *y = malloc(nx * sizeof(double));
  double *z = malloc(nx * sizeof(double));
  cec2014_sr_func(x, z, nx, Os, Mr, 600.0 / 100.0, s_flag, r_flag, y);
  for (i = 0; i < nx; i++) {
    s += z[i] * z[i];
    p *= cos(z[i] / sqrt(1.0 + i));
  }
  f[0] = 1.0 + s / 4000.0 - p;
  free(y);
  free(z);
}

void cec2014_rastrigin_func(double *x, double *f, int nx, double *Os,
                            double *Mr, int s_flag, int r_flag) {
  int i;
  f[0] = 0.0;
  double *y = malloc(nx * sizeof(double));
  double *z = malloc(nx * sizeof(double));
  cec2014_sr_func(x, z, nx, Os, Mr, 5.12 / 100.0, s_flag, r_flag, y);
  for (i = 0; i < nx; i++) {
    f[0] += (z[i] * z[i] - 10.0 * cos(2.0 * M_PI * z[i]) + 10.0);
  }
  free(y);
  free(z);
}

void cec2014_schwefel_func(double *x, double *f, int nx, double *Os, double *Mr,
                           int s_flag, int r_flag) {
  int i;
  double tmp;
  f[0] = 0.0;
  double *y = malloc(nx * sizeof(double));
  double *z = malloc(nx * sizeof(double));

  cec2014_sr_func(x, z, nx, Os, Mr, 1000.0 / 100.0, s_flag, r_flag, y);

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
  free(y);
  free(z);
}

void cec2014_katsuura_func(double *x, double *f, int nx, double *Os, double *Mr,
                           int s_flag, int r_flag) {
  int i, j;
  double temp, tmp1, tmp2, tmp3;
  f[0] = 1.0;
  tmp3 = pow(1.0 * nx, 1.2);
  double *y = malloc(nx * sizeof(double));
  double *z = malloc(nx * sizeof(double));
  cec2014_sr_func(x, z, nx, Os, Mr, 5.0 / 100.0, s_flag, r_flag, y);

  for (i = 0; i < nx; i++) {
    temp = 0.0;
    for (j = 1; j <= 32; j++) {
      tmp1 = pow(2.0, j);
      tmp2 = tmp1 * z[i];
      temp += fabs(tmp2 - floor(tmp2 + 0.5)) / tmp1;
    }
    f[0] *= pow(1.0 + (i + 1) * temp, 10.0 / tmp3);
  }
  tmp1 = 10.0 / nx / nx;
  f[0] = f[0] * tmp1 - tmp1;
  free(y);
  free(z);
}

void cec2014_grie_rosen_func(double *x, double *f, int nx, double *Os,
                             double *Mr, int s_flag, int r_flag) {
  int i;
  double temp, tmp1, tmp2;
  f[0] = 0.0;
  double *y = malloc(nx * sizeof(double));
  double *z = malloc(nx * sizeof(double));
  cec2014_sr_func(x, z, nx, Os, Mr, 5.0 / 100.0, s_flag, r_flag, y);
  z[0] += 1.0; // shift to orgin
  for (i = 0; i < nx - 1; i++) {
    z[i + 1] += 1.0; // shift to orgin
    tmp1 = z[i] * z[i] - z[i + 1];
    tmp2 = z[i] - 1.0;
    temp = 100.0 * tmp1 * tmp1 + tmp2 * tmp2;
    f[0] += (temp * temp) / 4000.0 - cos(temp) + 1.0;
  }
  tmp1 = z[nx - 1] * z[nx - 1] - z[0];
  tmp2 = z[nx - 1] - 1.0;
  temp = 100.0 * tmp1 * tmp1 + tmp2 * tmp2;
  f[0] += (temp * temp) / 4000.0 - cos(temp) + 1.0;
  free(y);
  free(z);
}

void cec2014_escaffer6_func(double *x, double *f, int nx, double *Os,
                            double *Mr, int s_flag, int r_flag) {
  int i;
  double temp1, temp2;
  double *y = malloc(nx * sizeof(double));
  double *z = malloc(nx * sizeof(double));
  cec2014_sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag, y); /* shift and rotate
                                                              */

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
  free(z);
  free(y);
}

void cec2014_happycat_func(double *x, double *f, int nx, double *Os, double *Mr,
                           int s_flag, int r_flag) {
  int i;
  double alpha, r2, sum_z;
  alpha = 1.0 / 8.0;
  double *y = malloc(nx * sizeof(double));
  double *z = malloc(nx * sizeof(double));
  cec2014_sr_func(x, z, nx, Os, Mr, 5.0 / 100.0, s_flag, r_flag, y);
  r2 = 0.0;
  sum_z = 0.0;
  for (i = 0; i < nx; i++) {
    z[i] = z[i] - 1.0; // shift to orgin
    r2 += z[i] * z[i];
    sum_z += z[i];
  }
  f[0] = pow(fabs(r2 - nx), 2 * alpha) + (0.5 * r2 + sum_z) / nx + 0.5;
  free(y);
  free(z);
}

void cec2014_hgbat_func(double *x, double *f, int nx, double *Os, double *Mr,
                        int s_flag, int r_flag) {
  int i;
  double alpha, r2, sum_z;
  alpha = 1.0 / 4.0;
  double *y = malloc(nx * sizeof(double));
  double *z = malloc(nx * sizeof(double));

  cec2014_sr_func(x, z, nx, Os, Mr, 5.0 / 100.0, s_flag, r_flag, y);
  r2 = 0.0;
  sum_z = 0.0;
  for (i = 0; i < nx; i++) {
    z[i] = z[i] - 1.0; // shift to orgin
    r2 += z[i] * z[i];
    sum_z += z[i];
  }
  f[0] = pow(fabs(pow(r2, 2.0) - pow(sum_z, 2.0)), 2 * alpha) +
         (0.5 * r2 + sum_z) / nx + 0.5;
  free(y);
  free(z);
}

void cec2014_hf01(double *x, double *f, int nx, double *Os, double *Mr, int *S,
          int s_flag, int r_flag) {
  int i, tmp, cf_num = 3;
  double fit[3];
  int G[3], G_nx[3];
  double Gp[3] = {0.3, 0.3, 0.4};
  tmp = 0;
  for (i = 0; i < cf_num - 1; i++) {
    G_nx[i] = ceil(Gp[i] * nx);
    tmp += G_nx[i];
  }
  G_nx[cf_num - 1] = nx - tmp;
  G[0] = 0;
  for (i = 1; i < cf_num; i++) {
    G[i] = G[i - 1] + G_nx[i - 1];
  }
  double *y = malloc(nx * sizeof(double));
  double *z = malloc(nx * sizeof(double));
  cec2014_sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag, y); /* shift and rotate */
  for (i = 0; i < nx; i++) {
    y[i] = z[S[i] - 1];
  }
  i = 0;
  cec2014_schwefel_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 1;
  cec2014_rastrigin_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 2;
  cec2014_ellips_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  f[0] = 0.0;
  for (i = 0; i < cf_num; i++) {
    f[0] += fit[i];
  }
  free(y);
  free(z);
}

void cec2014_hf02(double *x, double *f, int nx, double *Os, double *Mr, int *S,
          int s_flag, int r_flag) /* Hybrid Function 2 */
{
  int i, tmp, cf_num = 3;
  double fit[3];
  int G[3], G_nx[3];
  double Gp[3] = {0.3, 0.3, 0.4};

  tmp = 0;
  for (i = 0; i < cf_num - 1; i++) {
    G_nx[i] = ceil(Gp[i] * nx);
    tmp += G_nx[i];
  }
  G_nx[cf_num - 1] = nx - tmp;

  G[0] = 0;
  for (i = 1; i < cf_num; i++) {
    G[i] = G[i - 1] + G_nx[i - 1];
  }
  double *y = malloc(nx * sizeof(double));
  double *z = malloc(nx * sizeof(double));

  cec2014_sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag, y); /* shift and rotate */

  for (i = 0; i < nx; i++) {
    y[i] = z[S[i] - 1];
  }
  i = 0;
  cec2014_bent_cigar_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 1;
  cec2014_hgbat_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 2;
  cec2014_rastrigin_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);

  f[0] = 0.0;
  for (i = 0; i < cf_num; i++) {
    f[0] += fit[i];
  }
free(y);
  free(z);
}

void cec2014_hf03(double *x, double *f, int nx, double *Os, double *Mr, int *S,
          int s_flag, int r_flag) /* Hybrid Function 3 */
{
  int i, tmp, cf_num = 4;
  double fit[4];
  int G[4], G_nx[4];
  double Gp[4] = {0.2, 0.2, 0.3, 0.3};

  tmp = 0;
  for (i = 0; i < cf_num - 1; i++) {
    G_nx[i] = ceil(Gp[i] * nx);
    tmp += G_nx[i];
  }
  G_nx[cf_num - 1] = nx - tmp;

  G[0] = 0;
  for (i = 1; i < cf_num; i++) {
    G[i] = G[i - 1] + G_nx[i - 1];
  }
  double *y = malloc(nx * sizeof(double));
  double *z = malloc(nx * sizeof(double));

  cec2014_sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag, y); /* shift and rotate */

  for (i = 0; i < nx; i++) {
    y[i] = z[S[i] - 1];
  }
  i = 0;
  cec2014_griewank_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 1;
  cec2014_weierstrass_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 2;
  cec2014_rosenbrock_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 3;
  cec2014_escaffer6_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);

  f[0] = 0.0;
  for (i = 0; i < cf_num; i++) {
    f[0] += fit[i];
  }
free(y);
  free(z);
}

void cec2014_hf04(double *x, double *f, int nx, double *Os, double *Mr, int *S,
          int s_flag, int r_flag) /* Hybrid Function 4 */
{
  int i, tmp, cf_num = 4;
  double fit[4];
  int G[4], G_nx[4];
  double Gp[4] = {0.2, 0.2, 0.3, 0.3};

  tmp = 0;
  for (i = 0; i < cf_num - 1; i++) {
    G_nx[i] = ceil(Gp[i] * nx);
    tmp += G_nx[i];
  }
  G_nx[cf_num - 1] = nx - tmp;

  G[0] = 0;
  for (i = 1; i < cf_num; i++) {
    G[i] = G[i - 1] + G_nx[i - 1];
  }
  double *y = malloc(nx * sizeof(double));
  double *z = malloc(nx * sizeof(double));

  cec2014_sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag, y); /* shift and rotate */

  for (i = 0; i < nx; i++) {
    y[i] = z[S[i] - 1];
  }
  i = 0;
  cec2014_hgbat_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 1;
  cec2014_discus_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 2;
  cec2014_grie_rosen_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 3;
  cec2014_rastrigin_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);

  f[0] = 0.0;
  for (i = 0; i < cf_num; i++) {
    f[0] += fit[i];
  }
free(y);
  free(z);
}
void cec2014_hf05(double *x, double *f, int nx, double *Os, double *Mr, int *S,
          int s_flag, int r_flag) /* Hybrid Function 5 */
{
  int i, tmp, cf_num = 5;
  double fit[5];
  int G[5], G_nx[5];
  double Gp[5] = {0.1, 0.2, 0.2, 0.2, 0.3};

  tmp = 0;
  for (i = 0; i < cf_num - 1; i++) {
    G_nx[i] = ceil(Gp[i] * nx);
    tmp += G_nx[i];
  }
  G_nx[cf_num - 1] = nx - tmp;

  G[0] = 0;
  for (i = 1; i < cf_num; i++) {
    G[i] = G[i - 1] + G_nx[i - 1];
  }
  double *y = malloc(nx * sizeof(double));
  double *z = malloc(nx * sizeof(double));

  cec2014_sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag, y); /* shift and rotate */

  for (i = 0; i < nx; i++) {
    y[i] = z[S[i] - 1];
  }
  i = 0;
  cec2014_escaffer6_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 1;
  cec2014_hgbat_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 2;
  cec2014_rosenbrock_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 3;
  cec2014_schwefel_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 4;
  cec2014_ellips_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);

  f[0] = 0.0;
  for (i = 0; i < cf_num; i++) {
    f[0] += fit[i];
  }
free(y);
  free(z);
}

void cec2014_hf06(double *x, double *f, int nx, double *Os, double *Mr, int *S,
          int s_flag, int r_flag) /* Hybrid Function 6 */
{
  int i, tmp, cf_num = 5;
  double fit[5];
  int G[5], G_nx[5];
  double Gp[5] = {0.1, 0.2, 0.2, 0.2, 0.3};

  tmp = 0;
  for (i = 0; i < cf_num - 1; i++) {
    G_nx[i] = ceil(Gp[i] * nx);
    tmp += G_nx[i];
  }
  G_nx[cf_num - 1] = nx - tmp;

  G[0] = 0;
  for (i = 1; i < cf_num; i++) {
    G[i] = G[i - 1] + G_nx[i - 1];
  }
  double *y = malloc(nx * sizeof(double));
  double *z = malloc(nx * sizeof(double));

  cec2014_sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag, y); /* shift and rotate */

  for (i = 0; i < nx; i++) {
    y[i] = z[S[i] - 1];
  }
  i = 0;
  cec2014_katsuura_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 1;
  cec2014_happycat_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 2;
  cec2014_grie_rosen_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 3;
  cec2014_schwefel_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  i = 4;
  cec2014_ackley_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
  f[0] = 0.0;
  for (i = 0; i < cf_num; i++) {
    f[0] += fit[i];
  }
  free(y);
  free(z);
}

void cec2014_cf01(double *x, double *f, int nx, double *Os, double *Mr,
          int r_flag) /* Composition Function 1 */
{
  int i, cf_num = 5;
  double fit[5];
  double delta[5] = {10, 20, 30, 40, 50};
  double bias[5] = {0, 100, 200, 300, 400};

  i = 0;
  cec2014_rosenbrock_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000 * fit[i] / 1e+4;
  i = 1;
  cec2014_ellips_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000 * fit[i] / 1e+10;
  i = 2;
  cec2014_bent_cigar_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000 * fit[i] / 1e+30;
  i = 3;
  cec2014_discus_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000 * fit[i] / 1e+10;
  i = 4;
  cec2014_ellips_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, 0);
  fit[i] = 10000 * fit[i] / 1e+10;
  cec2014_cf_cal(x, f, nx, Os, delta, bias, fit, cf_num);
}

void cec2014_cf02(double *x, double *f, int nx, double *Os, double *Mr,
          int r_flag) /* Composition Function 2 */
{
  int i, cf_num = 3;
  double fit[3];
  double delta[3] = {20, 20, 20};
  double bias[3] = {0, 100, 200};

  i = 0;
  cec2014_schwefel_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, 0);
  i = 1;
  cec2014_rastrigin_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  i = 2;
  cec2014_hgbat_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  cec2014_cf_cal(x, f, nx, Os, delta, bias, fit, cf_num);
}

void cec2014_cf03(double *x, double *f, int nx, double *Os, double *Mr,
          int r_flag) /* Composition Function 3 */
{
  int i, cf_num = 3;
  double fit[3];
  double delta[3] = {10, 30, 50};
  double bias[3] = {0, 100, 200};
  i = 0;
  cec2014_schwefel_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 1000 * fit[i] / 4e+3;
  i = 1;
  cec2014_rastrigin_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 1000 * fit[i] / 1e+3;
  i = 2;
  cec2014_ellips_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 1000 * fit[i] / 1e+10;
  cec2014_cf_cal(x, f, nx, Os, delta, bias, fit, cf_num);
}

void cec2014_cf04(double *x, double *f, int nx, double *Os, double *Mr,
          int r_flag) /* Composition Function 4 */
{
  int i, cf_num = 5;
  double fit[5];
  double delta[5] = {10, 10, 10, 10, 10};
  double bias[5] = {0, 100, 200, 300, 400};
  i = 0;
  cec2014_schwefel_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 1000 * fit[i] / 4e+3;
  i = 1;
  cec2014_happycat_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 1000 * fit[i] / 1e+3;
  i = 2;
  cec2014_ellips_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 1000 * fit[i] / 1e+10;
  i = 3;
  cec2014_weierstrass_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 1000 * fit[i] / 400;
  i = 4;
  cec2014_griewank_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 1000 * fit[i] / 100;
  cec2014_cf_cal(x, f, nx, Os, delta, bias, fit, cf_num);
}

void cec2014_cf05(double *x, double *f, int nx, double *Os, double *Mr,
          int r_flag) /* Composition Function 4 */
{
  int i, cf_num = 5;
  double fit[5];
  double delta[5] = {10, 10, 10, 20, 20};
  double bias[5] = {0, 100, 200, 300, 400};
  i = 0;
  cec2014_hgbat_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000 * fit[i] / 1000;
  i = 1;
  cec2014_rastrigin_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000 * fit[i] / 1e+3;
  i = 2;
  cec2014_schwefel_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000 * fit[i] / 4e+3;
  i = 3;
  cec2014_weierstrass_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000 * fit[i] / 400;
  i = 4;
  cec2014_ellips_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000 * fit[i] / 1e+10;
  cec2014_cf_cal(x, f, nx, Os, delta, bias, fit, cf_num);
}

void cec2014_cf06(double *x, double *f, int nx, double *Os, double *Mr,
          int r_flag) /* Composition Function 6 */
{
  int i, cf_num = 5;
  double fit[5];
  double delta[5] = {10, 20, 30, 40, 50};
  double bias[5] = {0, 100, 200, 300, 400};
  i = 0;
  cec2014_grie_rosen_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000 * fit[i] / 4e+3;
  i = 1;
  cec2014_happycat_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000 * fit[i] / 1e+3;
  i = 2;
  cec2014_schwefel_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000 * fit[i] / 4e+3;
  i = 3;
  cec2014_escaffer6_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000 * fit[i] / 2e+7;
  i = 4;
  cec2014_ellips_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = 10000 * fit[i] / 1e+10;
  cec2014_cf_cal(x, f, nx, Os, delta, bias, fit, cf_num);
}

void cec2014_cf07(double *x, double *f, int nx, double *Os, double *Mr, int *SS,
          int r_flag) /* Composition Function 7 */
{
  int i, cf_num = 3;
  double fit[3];
  double delta[3] = {10, 30, 50};
  double bias[3] = {0, 100, 200};
  i = 0;
  cec2014_hf01(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], &SS[i * nx], 1, r_flag);
  i = 1;
  cec2014_hf02(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], &SS[i * nx], 1, r_flag);
  i = 2;
  cec2014_hf03(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], &SS[i * nx], 1, r_flag);
  cec2014_cf_cal(x, f, nx, Os, delta, bias, fit, cf_num);
}

void cec2014_cf08(double *x, double *f, int nx, double *Os, double *Mr, int *SS,
          int r_flag) /* Composition Function 8 */
{
  int i, cf_num = 3;
  double fit[3];
  double delta[3] = {10, 30, 50};
  double bias[3] = {0, 100, 200};
  i = 0;
  cec2014_hf04(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], &SS[i * nx], 1, r_flag);
  i = 1;
  cec2014_hf05(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], &SS[i * nx], 1, r_flag);
  i = 2;
  cec2014_hf06(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], &SS[i * nx], 1, r_flag);
  cec2014_cf_cal(x, f, nx, Os, delta, bias, fit, cf_num);
}

void cec2014_asyfunc(double *x, double *xasy, int nx, double beta) {
  int i;
  for (i = 0; i < nx; i++) {
    if (x[i] > 0)
      xasy[i] = pow(x[i], 1.0 + beta * i / (nx - 1) * pow(x[i], 0.5));
  }
}

void cec2014_oszfunc(double *x, double *xosz, int nx) {
  int i, sx;
  double c1, c2;
  double xx = 0.0;
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

void cec2014_cf_cal(double *x, double *f, int nx, double *Os, double *delta,
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
