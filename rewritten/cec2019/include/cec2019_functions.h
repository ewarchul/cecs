#ifndef CEC2019_FUNCTIONS_H
#define CEC2019_FUNCTIONS_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define INF 1.0e99
#define EPS 1.0e-14
#define E 2.7182818284590452353602874713526625
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795029
#endif

void Lennard_Jones(double *, int, double *); /* Lennard Jones */
void Hilbert(double *, int, double *);       /* Hilbert */
void Chebyshev(double *, int, double *);     /* Chebyshev */
void cec2019_schaffer_F7_func(double *, double *, int, double *, double *, int,
                      int); /* Schwefel's F7 */
void cec2019_ackley_func(double *, double *, int, double *, double *, int,
                 int); /* Ackley's */
void cec2019_rastrigin_func(double *, double *, int, double *, double *, int,
                    int); /* Rastrigin's  */
void cec2019_weierstrass_func(double *, double *, int, double *, double *, int,
                      int); /* Weierstrass's  */
void cec2019_schwefel_func(double *, double *, int, double *, double *, int,
                   int); /* Schwefel's */
void cec2019_escaffer6_func(double *, double *, int, double *, double *, int,
                    int); /* Expanded Scaffer¡¯s F6  */
void cec2019_happycat_func(double *, double *, int, double *, double *, int,
                   int); /* HappyCat */
void cec2019_griewank_func(double *, double *, int, double *, double *, int,
                   int); /* Griewank's  */

void cec2019_shiftfunc(double *, double *, int, double *);
void cec2019_rotatefunc(double *, double *, int, double *);
void cec2019_sr_func(double *, double *, int, double *, double *, double, int,
             int, double *); /* shift and rotate */
void cec2019_asyfunc(double *, double *x, int, double);
void cec2019_oszfunc(double *, double *, int);

#endif // CEC2019_FUNCTIONS_H
