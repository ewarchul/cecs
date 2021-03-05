#ifndef CEC2015_FUNCTIONS_H
#define CEC2015_FUNCTIONS_H

void cec2015_sphere_func(double *, double *, int, double *, double *, int,
                 int); /* Sphere */
void cec2015_ellips_func(double *, double *, int, double *, double *, int,
                 int); /* Ellipsoidal */
void cec2015_bent_cigar_func(double *, double *, int, double *, double *, int,
                     int); /* Bent_Cigar */
void cec2015_discus_func(double *, double *, int, double *, double *, int,
                 int); /* Discus */
void cec2015_dif_powers_func(double *, double *, int, double *, double *, int,
                     int); /* Different Powers */
void cec2015_rosenbrock_func(double *, double *, int, double *, double *, int,
                     int); /* Rosenbrock's */
void cec2015_schaffer_F7_func(double *, double *, int, double *, double *, int,
                      int); /* Schwefel's F7 */
void cec2015_ackley_func(double *, double *, int, double *, double *, int,
                 int); /* Ackley's */
void cec2015_rastrigin_func(double *, double *, int, double *, double *, int,
                    int); /* Rastrigin's  */
void cec2015_weierstrass_func(double *, double *, int, double *, double *, int,
                      int); /* Weierstrass's  */
void cec2015_griewank_func(double *, double *, int, double *, double *, int,
                   int); /* Griewank's  */
void cec2015_schwefel_func(double *, double *, int, double *, double *, int,
                   int); /* Schwefel's */
void cec2015_katsuura_func(double *, double *, int, double *, double *, int,
                   int); /* Katsuura */
void cec2015_bi_rastrigin_func(double *, double *, int, double *, double *, int,
                       int); /* Lunacek Bi_rastrigin */
void cec2015_grie_rosen_func(double *, double *, int, double *, double *, int,
                     int); /* Griewank-Rosenbrock  */
void cec2015_escaffer6_func(double *, double *, int, double *, double *, int,
                    int); /* Expanded Scaffer¡¯s F6  */
void cec2015_step_rastrigin_func(double *, double *, int, double *, double *, int,
                         int); /* Noncontinuous Rastrigin's  */
void cec2015_happycat_func(double *, double *, int, double *, double *, int,
                   int); /* HappyCat */
void cec2015_hgbat_func(double *, double *, int, double *, double *, int,
                int); /* HGBat  */
void cec2015_hf01(double *, double *, int, double *, double *, int *, int,
          int); /* Hybrid Function 1 */
void cec2015_hf02(double *, double *, int, double *, double *, int *, int,
          int); /* Hybrid Function 2 */
void cec2015_hf03(double *, double *, int, double *, double *, int *, int,
          int); /* Hybrid Function 3 */
void cec2015_cf01(double *, double *, int, double *, double *, double *,
          int); /* Composition Function 1 */
void cec2015_cf02(double *, double *, int, double *, double *, int *, double *,
          int); /* Composition Function 2 */
void cec2015_cf03(double *, double *, int, double *, double *, double *,
          int); /* Composition Function 3 */
void cec2015_cf04(double *, double *, int, double *, double *, double *,
          int); /* Composition Function 4 */
void cec2015_cf05(double *, double *, int, double *, double *, int *, double *,
          int); /* Composition Function 5 */
void cec2015_cf06(double *, double *, int, double *, double *, double *,
          int); /* Composition Function 6 */
void cec2015_cf07(double *, double *, int, double *, double *, double *,
          int); /* Composition Function 7 */
void cec2015_shiftfunc(double *, double *, int, double *);
void cec2015_rotatefunc(double *, double *, int, double *);
void cec2015_sr_func(double *, double *, int, double *, double *, double, int,
             int); /* shift and rotate */
void cec2015_cf_cal(double *, double *, int, double *, double *, double *, double *,
            int);
#endif // CEC2015_FUNCTION_H
