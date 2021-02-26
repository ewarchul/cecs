#ifndef CEC2021_FUNCTIONS_H_
#define CEC2021_FUNCTIONS_H_

#include <math.h>
#include <stdio.h>

#include "../consts.h"
#include "../globals.h"


void cec2021_ellips_func(double *, double *, int, double *, double *, int,
                         int); /* Ellipsoidal */
void cec2021_bent_cigar_func(double *, double *, int, double *, double *, int,
                             int); /* Discus */
void cec2021_discus_func(double *, double *, int, double *, double *, int,
                         int); /* Bent_Cigar */
void cec2021_rosenbrock_func(double *, double *, int, double *, double *, int,
                             int); /* Rosenbrock's */
void cec2021_ackley_func(double *, double *, int, double *, double *, int,
                         int); /* Ackley's */
void cec2021_rastrigin_func(double *, double *, int, double *, double *, int,
                            int); /* Rastrigin's  */
void cec2021_griewank_func(double *, double *, int, double *, double *, int,
                           int); /* Griewank's  */
void cec2021_schwefel_func(double *, double *, int, double *, double *, int,
                           int); /* Schwefel's */
void cec2021_bi_rastrigin_func(double *, double *, int, double *, double *, int,
                               int); /* Lunacek Bi_rastrigin */
void cec2021_grie_rosen_func(double *, double *, int, double *, double *, int,
                             int); /* Griewank-Rosenbrock  */
void cec2021_escaffer6_func(double *, double *, int, double *, double *, int,
                            int); /* Expanded Scaffer��s F6  */
void cec2021_happycat_func(double *, double *, int, double *, double *, int,
                           int); /* HappyCat */
void cec2021_hgbat_func(double *, double *, int, double *, double *, int,
                        int); /* HGBat  */
void cec2021_hf01(double *, double *, int, double *, double *, int *, int,
                  int); /* Hybrid Function 1 */
void cec2021_hf05(double *, double *, int, double *, double *, int *, int,
                  int); /* Hybrid Function 5 */
void cec2021_hf06(double *, double *, int, double *, double *, int *, int,
                  int); /* Hybrid Function 6 */
void cec2021_cf02(double *, double *, int, double *, double *,
                  int); /* Composition Function 2 */
void cec2021_cf04(double *, double *, int, double *, double *,
                  int); /* Composition Function 4 */
void cec2021_cf05(double *, double *, int, double *, double *,
                  int); /* Composition Function 5 */
void cec2021_cf02_s(double *, double *, int, double *, double *,
                    int); /* Composition Function 2 for shift case*/
void cec2021_cf04_s(double *, double *, int, double *, double *,
                    int); /* Composition Function 4 for shift case*/
void cec2021_cf05_s(double *, double *, int, double *, double *,
                    int); /* Composition Function 5 for shift case*/
void cec2021_shiftfunc(double *, double *, int, double *);
void cec2021_rotatefunc(double *, double *, int, double *);
void cec2021_sr_func(double *, double *, int, double *, double *, double, int,
                     int); /* shift and rotate */
void cec2021_asyfunc(double *, double *x, int, double);
void cec2021_oszfunc(double *, double *, int);
void cec2021_cf_cal(double *, double *, int, double *, double *, double *,
                    double *, int);

#endif // CEC2021_FUNCTIONS_H_
