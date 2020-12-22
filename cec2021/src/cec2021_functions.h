#ifndef CEC2021_FUNCTIONS_H_
#define CEC2021_FUNCTIONS_H_

#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include "consts.h"
#include "interface.h"

void ellips_func(double *, double *, int , double *,double *, int, int); /* Ellipsoidal */
void bent_cigar_func(double *, double *, int , double *,double *, int, int); /* Discus */
void discus_func(double *, double *, int , double *,double *, int, int);  /* Bent_Cigar */
void rosenbrock_func (double *, double *, int , double *,double *, int, int); /* Rosenbrock's */
void ackley_func (double *, double *, int , double *,double *, int, int); /* Ackley's */
void rastrigin_func (double *, double *, int , double *,double *, int, int); /* Rastrigin's  */
void griewank_func (double *, double *, int , double *,double *, int, int); /* Griewank's  */
void schwefel_func (double *, double *, int , double *,double *, int, int); /* Schwefel's */
void bi_rastrigin_func (double *, double *, int , double *,double *, int, int); /* Lunacek Bi_rastrigin */
void grie_rosen_func (double *, double *, int , double *,double *, int, int); /* Griewank-Rosenbrock  */
void escaffer6_func (double *, double *, int , double *,double *, int, int); /* Expanded Scaffer��s F6  */
void happycat_func (double *, double *, int , double *,double *, int, int); /* HappyCat */
void hgbat_func (double *, double *, int , double *,double *, int, int); /* HGBat  */
void hf01 (double *, double *, int, double *,double *, int *,int, int); /* Hybrid Function 1 */
void hf05 (double *, double *, int, double *,double *, int *,int, int); /* Hybrid Function 5 */
void hf06 (double *, double *, int, double *,double *, int *,int, int); /* Hybrid Function 6 */
void cf02 (double *, double *, int , double *,double *, int); /* Composition Function 2 */
void cf04 (double *, double *, int , double *,double *, int); /* Composition Function 4 */
void cf05 (double *, double *, int , double *,double *, int); /* Composition Function 5 */
void cf02_s (double *, double *, int , double *,double *, int); /* Composition Function 2 for shift case*/
void cf04_s (double *, double *, int , double *,double *, int); /* Composition Function 4 for shift case*/
void cf05_s (double *, double *, int , double *,double *, int); /* Composition Function 5 for shift case*/
void shiftfunc (double*,double*,int,double*);
void rotatefunc (double*,double*,int, double*);
void sr_func (double *, double *, int, double*, double*, double, int, int); /* shift and rotate */
void asyfunc (double *, double *x, int, double);
void oszfunc (double *, double *, int);
void cf_cal(double *, double *, int, double *,double *,double *,double *,int);

#endif // CEC2021_FUNCTIONS_H_