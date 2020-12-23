#ifndef CEC2014_FUNCTIONS_H_
#define CEC2014_FUNCTIONS_H_

#include "consts.h"
#include <malloc.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "globals.h"
void cec2014_sphere_func(double *, double *, int, double *, double *, int,
                         int); /* Sphere */
void cec2014_ellips_func(double *, double *, int, double *, double *, int,
                         int); /* Ellipsoidal */
void cec2014_bent_cigar_func(double *, double *, int, double *, double *, int,
                             int); /* Discus */
void cec2014_discus_func(double *, double *, int, double *, double *, int,
                         int); /* Bent_Cigar */
void cec2014_dif_powers_func(double *, double *, int, double *, double *, int,
                             int); /* Different Powers */
void cec2014_rosenbrock_func(double *, double *, int, double *, double *, int,
                             int); /* Rosenbrock's */
void cec2014_schaffer_F7_func(double *, double *, int, double *, double *, int,
                              int); /* Schwefel's F7 */
void cec2014_ackley_func(double *, double *, int, double *, double *, int,
                         int); /* Ackley's */
void cec2014_rastrigin_func(double *, double *, int, double *, double *, int,
                            int); /* Rastrigin's  */
void cec2014_weierstrass_func(double *, double *, int, double *, double *, int,
                              int); /* Weierstrass's  */
void cec2014_griewank_func(double *, double *, int, double *, double *, int,
                           int); /* Griewank's  */
void cec2014_schwefel_func(double *, double *, int, double *, double *, int,
                           int); /* Schwefel's */
void cec2014_katsuura_func(double *, double *, int, double *, double *, int,
                           int); /* Katsuura */
void cec2014_bi_rastrigin_func(double *, double *, int, double *, double *, int,
                               int); /* Lunacek Bi_rastrigin */
void cec2014_grie_rosen_func(double *, double *, int, double *, double *, int,
                             int); /* Griewank-Rosenbrock  */
void cec2014_escaffer6_func(double *, double *, int, double *, double *, int,
                            int); /* Expanded Scaffer??s F6  */
void cec2014_step_rastrigin_func(double *, double *, int, double *, double *,
                                 int, int); /* Noncontinuous Rastrigin's  */
void cec2014_happycat_func(double *, double *, int, double *, double *, int,
                           int); /* HappyCat */
void cec2014_hgbat_func(double *, double *, int, double *, double *, int,
                        int); /* HGBat  */

void cec2014_hf01(double *, double *, int, double *, double *, int *, int,
                  int); /* Hybrid Function 1 */
void cec2014_hf02(double *, double *, int, double *, double *, int *, int,
                  int); /* Hybrid Function 2 */
void cec2014_hf03(double *, double *, int, double *, double *, int *, int,
                  int); /* Hybrid Function 3 */
void cec2014_hf04(double *, double *, int, double *, double *, int *, int,
                  int); /* Hybrid Function 4 */
void cec2014_hf05(double *, double *, int, double *, double *, int *, int,
                  int); /* Hybrid Function 5 */
void cec2014_hf06(double *, double *, int, double *, double *, int *, int,
                  int); /* Hybrid Function 6 */

void cec2014_cf01(double *, double *, int, double *, double *,
                  int); /* Composition Function 1 */
void cec2014_cf02(double *, double *, int, double *, double *,
                  int); /* Composition Function 2 */
void cec2014_cf03(double *, double *, int, double *, double *,
                  int); /* Composition Function 3 */
void cec2014_cf04(double *, double *, int, double *, double *,
                  int); /* Composition Function 4 */
void cec2014_cf05(double *, double *, int, double *, double *,
                  int); /* Composition Function 5 */
void cec2014_cf06(double *, double *, int, double *, double *,
                  int); /* Composition Function 6 */
void cec2014_cf07(double *, double *, int, double *, double *, int *,
                  int); /* Composition Function 7 */
void cec2014_cf08(double *, double *, int, double *, double *, int *,
                  int); /* Composition Function 8 */

void cec2014_shiftfunc(double *, double *, int, double *);
void cec2014_rotatefunc(double *, double *, int, double *);
void cec2014_sr_func(double *, double *, int, double *, double *, double, int,
                     int); /* shift and rotate */
void cec2014_asyfunc(double *, double *x, int, double);
void cec2014_oszfunc(double *, double *, int);
void cec2014_cf_cal(double *, double *, int, double *, double *, double *,
                    double *, int);

#endif // CEC2014_FUNCTIONS_H_
