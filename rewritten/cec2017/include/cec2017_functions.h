#ifndef CEC2017_FUNCTIONS_H
#define CEC2017_FUNCTIONS_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define INF 1.0e99
#define EPS 1.0e-14
#define E 2.7182818284590452353602874713526625
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795029
#endif


void cec2017_sphere_func(double *, double *, int, double *, double *, int,
                         int); /* Sphere */
void cec2017_ellips_func(double *, double *, int, double *, double *, int, int);
void cec2017_bent_cigar_func(double *, double *, int, double *, double *, int,
                             int);
void cec2017_discus_func(double *, double *, int, double *, double *, int, int);
void cec2017_dif_powers_func(double *, double *, int, double *, double *, int,
                             int);
void cec2017_rosenbrock_func(double *, double *, int, double *, double *, int,
                             int);
void cec2017_schaffer_F7_func(double *, double *, int, double *, double *, int,
                              int);
void cec2017_ackley_func(double *, double *, int, double *, double *, int, int);
void cec2017_rastrigin_func(double *, double *, int, double *, double *, int,
                            int);
void cec2017_weierstrass_func(double *, double *, int, double *, double *, int,
                              int);
void cec2017_griewank_func(double *, double *, int, double *, double *, int,
                           int);
void cec2017_schwefel_func(double *, double *, int, double *, double *, int,
                           int);
void cec2017_katsuura_func(double *, double *, int, double *, double *, int,
                           int);
void cec2017_bi_rastrigin_func(double *, double *, int, double *, double *, int,
                               int);
void cec2017_grie_rosen_func(double *, double *, int, double *, double *, int,
                             int);
void cec2017_escaffer6_func(double *, double *, int, double *, double *, int,
                            int);
void cec2017_step_rastrigin_func(double *, double *, int, double *, double *,
                                 int, int);
void cec2017_happycat_func(double *, double *, int, double *, double *, int,
                           int);
void cec2017_hgbat_func(double *, double *, int, double *, double *, int, int);

void cec2017_sum_diff_pow_func(double *, double *, int, double *, double *, int,
                               int);
void cec2017_zakharov_func(double *, double *, int, double *, double *, int,
                           int);
void cec2017_levy_func(double *, double *, int, double *, double *, int, int);
void cec2017_dixon_price_func(double *, double *, int, double *, double *, int,
                              int);

void cec2017_hf01(double *, double *, int, double *, double *, int *, int, int);
void cec2017_hf02(double *, double *, int, double *, double *, int *, int, int);
void cec2017_hf03(double *, double *, int, double *, double *, int *, int, int);
void cec2017_hf04(double *, double *, int, double *, double *, int *, int, int);
void cec2017_hf05(double *, double *, int, double *, double *, int *, int, int);
void cec2017_hf06(double *, double *, int, double *, double *, int *, int, int);
void cec2017_hf07(double *, double *, int, double *, double *, int *, int, int);
void cec2017_hf08(double *, double *, int, double *, double *, int *, int, int);
void cec2017_hf09(double *, double *, int, double *, double *, int *, int, int);
void cec2017_hf10(double *, double *, int, double *, double *, int *, int, int);

void cec2017_cf01(double *, double *, int, double *, double *, int);
void cec2017_cf02(double *, double *, int, double *, double *, int);
void cec2017_cf03(double *, double *, int, double *, double *, int);
void cec2017_cf04(double *, double *, int, double *, double *, int);
void cec2017_cf05(double *, double *, int, double *, double *, int);
void cec2017_cf06(double *, double *, int, double *, double *, int);
void cec2017_cf07(double *, double *, int, double *, double *, int);
void cec2017_cf08(double *, double *, int, double *, double *, int);
void cec2017_cf09(double *, double *, int, double *, double *, int *, int);
void cec2017_cf10(double *, double *, int, double *, double *, int *, int);

void cec2017_shiftfunc(double *, double *, int, double *);
void cec2017_rotatefunc(double *, double *, int, double *);
void cec2017_sr_func(double *, double *, int, double *, double *, double, int,
                     int, double *);
void cec2017_asyfunc(double *, double *x, int, double);
void cec2017_oszfunc(double *, double *, int);
void cec2017_cf_cal(double *, double *, int, double *, double *, double *,
                    double *, int);

#endif // CEC2017_FUNCTIONS_H
