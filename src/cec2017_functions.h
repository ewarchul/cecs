#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <malloc.h>

void cec2017_sphere_func (double *, double *, int , double *,double *, int, int); /* Sphere */
void cec2017_ellips_func(double *, double *, int , double *,double *, int, int); /* Ellipsoidal */
void cec2017_bent_cigar_func(double *, double *, int , double *,double *, int, int); /* Discus */
void cec2017_discus_func(double *, double *, int , double *,double *, int, int);  /* Bent_Cigar */
void cec2017_dif_powers_func(double *, double *, int , double *,double *, int, int);  /* Different Powers */
void cec2017_rosenbrock_func (double *, double *, int , double *,double *, int, int); /* Rosenbrock's */
void cec2017_schaffer_F7_func (double *, double *, int , double *,double *, int, int); /* Schwefel's F7 */
void cec2017_ackley_func (double *, double *, int , double *,double *, int, int); /* Ackley's */
void cec2017_rastrigin_func (double *, double *, int , double *,double *, int, int); /* Rastrigin's  */
void cec2017_weierstrass_func (double *, double *, int , double *,double *, int, int); /* Weierstrass's  */
void cec2017_griewank_func (double *, double *, int , double *,double *, int, int); /* Griewank's  */
void cec2017_schwefel_func (double *, double *, int , double *,double *, int, int); /* Schwefel's */
void cec2017_katsuura_func (double *, double *, int , double *,double *, int, int); /* Katsuura */
void cec2017_bi_rastrigin_func (double *, double *, int , double *,double *, int, int); /* Lunacek Bi_rastrigin */
void cec2017_grie_rosen_func (double *, double *, int , double *,double *, int, int); /* Griewank-Rosenbrock  */
void cec2017_escaffer6_func (double *, double *, int , double *,double *, int, int); /* Expanded Scaffer¡¯s F6  */
void cec2017_step_rastrigin_func (double *, double *, int , double *,double *, int, int); /* Noncontinuous Rastrigin's  */
void cec2017_happycat_func (double *, double *, int , double *,double *, int, int); /* HappyCat */
void cec2017_hgbat_func (double *, double *, int , double *,double *, int, int); /* HGBat  */

/* New functions Noor Changes */
void cec2017_sum_diff_pow_func(double *, double *, int , double *,double *, int, int); /* Sum of different power */
void cec2017_zakharov_func(double *, double *, int , double *,double *, int, int); /* ZAKHAROV */
void cec2017_levy_func(double *, double *, int , double *,double *, int, int); /* Levy */
void cec2017_dixon_price_func(double *, double *, int , double *,double *, int, int); /* Dixon and Price */

void cec2017_hf01 (double *, double *, int, double *,double *, int *,int, int); /* Hybrid Function 1 */
void cec2017_hf02 (double *, double *, int, double *,double *, int *,int, int); /* Hybrid Function 2 */
void cec2017_hf03 (double *, double *, int, double *,double *, int *,int, int); /* Hybrid Function 3 */
void cec2017_hf04 (double *, double *, int, double *,double *, int *,int, int); /* Hybrid Function 4 */
void cec2017_hf05 (double *, double *, int, double *,double *, int *,int, int); /* Hybrid Function 5 */
void cec2017_hf06 (double *, double *, int, double *,double *, int *,int, int); /* Hybrid Function 6 */
void cec2017_hf07 (double *, double *, int, double *,double *, int *,int, int); /* Hybrid Function 7 */
void cec2017_hf08 (double *, double *, int, double *,double *, int *,int, int); /* Hybrid Function 8 */
void cec2017_hf09 (double *, double *, int, double *,double *, int *,int, int); /* Hybrid Function 9 */
void cec2017_hf10 (double *, double *, int, double *,double *, int *,int, int); /* Hybrid Function 10 */

void cec2017_cf01 (double *, double *, int , double *,double *, int); /* Composition Function 1 */
void cec2017_cf02 (double *, double *, int , double *,double *, int); /* Composition Function 2 */
void cec2017_cf03 (double *, double *, int , double *,double *, int); /* Composition Function 3 */
void cec2017_cf04 (double *, double *, int , double *,double *, int); /* Composition Function 4 */
void cec2017_cf05 (double *, double *, int , double *,double *, int); /* Composition Function 5 */
void cec2017_cf06 (double *, double *, int , double *,double *, int); /* Composition Function 6 */
void cec2017_cf07 (double *, double *, int , double *,double *, int); /* Composition Function 7 */
void cec2017_cf08 (double *, double *, int , double *,double *, int); /* Composition Function 8 */
void cec2017_cf09 (double *, double *, int , double *,double *, int *, int); /* Composition Function 9 */
void cec2017_cf10 (double *, double *, int , double *,double *, int *, int); /* Composition Function 10 */

void cec2017_shiftfunc (double*,double*,int,double*);
void cec2017_rotatefunc (double*,double*,int, double*);
void cec2017_sr_func (double *, double *, int, double*, double*, double, int, int); /* shift and rotate */
void cec2017_asyfunc (double *, double *x, int, double);
void cec2017_oszfunc (double *, double *, int);
void cec2017_cf_cal(double *, double *, int, double *,double *,double *,double *,int);