
/*
  CEC15 Test Function Suite for Single Objective Optimization
  Jane Jing Liang (email: liangjing@zzu.edu.cn; liangjing@pmail.ntu.edu.cn) 
  Nov. 12th 2014
*/

#include <WINDOWS.H>    
#include <stdio.h>
#include <math.h>
#include <malloc.h>


void cec15_nich_func(double *, double *,int,int,int);

double *OShift,*M,*y,*z,*x_bound;
int ini_flag=0,n_flag,func_flag;

void main()
{
	int i,j,k,n,m,func_num;
	double *f,*x;
	FILE *fpt;
	char FileName[30];
	m=2;
	n=2;
	x=(double *)malloc(m*n*sizeof(double));
	f=(double *)malloc(sizeof(double)  *  m);
	for (i = 0; i < 15; i++)
	{
		func_num=i+1;
		sprintf(FileName, "input_data/shift_data_%d.txt", func_num);
		fpt = fopen(FileName,"r");
		if (fpt==NULL)
		{
			printf("\n Error: Cannot open input file for reading \n");
		}
		
		if (x==NULL)
			printf("\nError: there is insufficient memory available!\n");

		for(k=0;k<n;k++)
		{
				fscanf(fpt,"%Lf",&x[k]);
				/*printf("%Lf\n",x[k]);*/
		}

		fclose(fpt);

			for (j = 0; j < n; j++)
			{
				x[1*n+j]=0.0;
				/*printf("%Lf\n",x[1*n+j]);*/
			}
		
		
		for (k = 0; k < 1; k++)
		{
			cec15_nich_func(x, f, n,m,func_num);
			for (j = 0; j < 2; j++)
			{
				printf(" f%d(x[%d]) = %Lf,",func_num,j+1,f[j]);
			}
			printf("\n");
		}
	
	}
	free(x);
	free(f);
	free(y);
	free(z);
	free(M);
	free(OShift);
	free(x_bound);
}




/*
  CEC15 Niching Test Function Suite
  Jane Jing Liang (email: liangjing@zzu.edu.cn)
  Nov. 21th 2014
  Modified on Dec.26th 2014 (Load shift_data)

  Reference:
  B. Y. Qu, J. J. Liang, P. N. Suganthan, Q. Chen, "Problem Definitions and
  Evaluation Criteria for the CEC 2015 Special Session and Competition on
  Niching Numerical Optimization",Technical Report201411B,Computational
  Intelligence Laboratory, Zhengzhou University, Zhengzhou China and Technical
  Report, Nanyang Technological University, Singapore, November 2014
*/

#include <WINDOWS.H>
#include <malloc.h>
#include <math.h>
#include <stdio.h>

#include <WINDOWS.H>
#include <malloc.h>
#include <math.h>
#include <stdio.h>

#define INF 1.0e99
#define EPS 1.0e-14
#define E 2.7182818284590452353602874713526625
#define PI 3.1415926535897932384626433832795029

void twopeaks_func(double *, double *, int, double *, double *, int,
                   int); /* Expanded Two-Peak Trap */
void fiveuneven_func(double *, double *, int, double *, double *, int,
                     int); /* Expanded Five-Uneven-Peak Trap */
void equalmin_func(double *, double *, int, double *, double *, int,
                   int); /* Expanded Equal Minima */
void decreasemin_func(double *, double *, int, double *, double *, int,
                      int); /* Expanded Decreasing Minima  */
void unevenmin_func(double *, double *, int, double *, double *, int,
                    int); /* Expanded Uneven Minima  */
void himmelblau_func(double *, double *, int, double *, double *, int,
                     int); /* Expanded Himmelblau¡¯s Function */
void camelback_func(double *, double *, int, double *, double *, int,
                    int); /* Expanded Six-Hump Camel Back */
void vincent_func(double *, double *, int, double *, double *, int,
                  int); /* Modified Vincent Function */

void sphere_func(double *, double *, int, double *, double *, int,
                 int); /* Sphere */
void ellips_func(double *, double *, int, double *, double *, int,
                 int); /* Ellipsoidal */
void bent_cigar_func(double *, double *, int, double *, double *, int,
                     int); /* Bent_Cigar */
void discus_func(double *, double *, int, double *, double *, int,
                 int); /* Discus */
void dif_powers_func(double *, double *, int, double *, double *, int,
                     int); /* Different Powers */
void rosenbrock_func(double *, double *, int, double *, double *, int,
                     int); /* Rosenbrock's */
void schaffer_F7_func(double *, double *, int, double *, double *, int,
                      int); /* Schwefel's F7 */
void ackley_func(double *, double *, int, double *, double *, int,
                 int); /* Ackley's */
void rastrigin_func(double *, double *, int, double *, double *, int,
                    int); /* Rastrigin's  */
void weierstrass_func(double *, double *, int, double *, double *, int,
                      int); /* Weierstrass's  */
void griewank_func(double *, double *, int, double *, double *, int,
                   int); /* Griewank's  */
void schwefel_func(double *, double *, int, double *, double *, int,
                   int); /* Schwefel's */
void katsuura_func(double *, double *, int, double *, double *, int,
                   int); /* Katsuura */
void bi_rastrigin_func(double *, double *, int, double *, double *, int,
                       int); /* Lunacek Bi_rastrigin */
void grie_rosen_func(double *, double *, int, double *, double *, int,
                     int); /* Griewank-Rosenbrock  */
void escaffer6_func(double *, double *, int, double *, double *, int,
                    int); /* Expanded Scaffer¡¯s F6  */
void step_rastrigin_func(double *, double *, int, double *, double *, int,
                         int); /* Noncontinuous Rastrigin's  */
void happycat_func(double *, double *, int, double *, double *, int,
                   int); /* HappyCat */
void hgbat_func(double *, double *, int, double *, double *, int,
                int); /* HGBat  */

void cf01(double *, double *, int, double *, double *,
          int); /* Composition Function 1 */
void cf02(double *, double *, int, double *, double *,
          int); /* Composition Function 2 */
void cf03(double *, double *, int, double *, double *,
          int); /* Composition Function 3 */
void cf04(double *, double *, int, double *, double *,
          int); /* Composition Function 4 */
void cf05(double *, double *, int, double *, double *,
          int); /* Composition Function 5 */
void cf06(double *, double *, int, double *, double *,
          int); /* Composition Function 6 */
void cf07(double *, double *, int, double *, double *,
          int); /* Composition Function 7 */

void shiftfunc(double *, double *, int, double *);
void rotatefunc(double *, double *, int, double *);
void sr_func(double *, double *, int, double *, double *, double, int,
             int); /* shift and rotate */
void cf_cal(double *, double *, int, double *, double *, double *, double *,
            int);
void cec15_nich_func(double *, double *, int, int, int);

extern double *OShift, *M, *y, *z, *x_bound;
extern int ini_flag, n_flag, func_flag;

int cf_nums[] = {0, 1, 1, 1, 1, 1, 1, 1, 1, 10, 10, 10, 10, 10, 10, 10};

void cec15_nich_func(double *x, double *f, int nx, int mx, int func_num) {
  int cf_num = 10, i, j;
  if (ini_flag == 1) {
    if ((n_flag != nx) || (func_flag != func_num)) {
      ini_flag = 0;
    }
  }

  if (ini_flag == 0) {
    FILE *fpt;
    char FileName[256];
    free(M);
    free(OShift);
    free(y);
    free(z);
    free(x_bound);
    y = (double *)malloc(sizeof(double) * nx);
    z = (double *)malloc(sizeof(double) * nx);
    x_bound = (double *)malloc(sizeof(double) * nx);
    for (i = 0; i < nx; i++)
      x_bound[i] = 100.0;

    /* Load Matrix M*/
    sprintf(FileName, "input_data/M_%d_D%d.txt", func_num, nx);
    fpt = fopen(FileName, "r");
    if (fpt == NULL) {
      mexErrMsgTxt("\n Error: Cannot open input file for reading \n");
    }
    M = (double *)malloc(cf_nums[func_num] * nx * nx * sizeof(double));
    if (M == NULL)
      mexErrMsgTxt("\nError: there is insufficient memory available!\n");
    for (i = 0; i < cf_nums[func_num] * nx * nx; i++)
      fscanf(fpt, "%lf", &M[i]);
    fclose(fpt);

    /* Load shift_data */
    sprintf(FileName, "input_data/shift_data_%d.txt", func_num);
    fpt = fopen(FileName, "r");
    if (fpt == NULL) {
      mexErrMsgTxt("\n Error: Cannot open input file for reading \n");
    }
    OShift = (double *)malloc(nx * cf_nums[func_num] * sizeof(double));
    if (OShift == NULL)
      mexErrMsgTxt("\nError: there is insufficient memory available!\n");
    if (func_num < 9) {
      for (i = 0; i < nx * cf_nums[func_num]; i++) {
        fscanf(fpt, "%Lf", &OShift[i]);
      }
    } else {

      for (i = 0; i < cf_nums[func_num] - 1; i++) {
        for (j = 0; j < nx; j++) {
          fscanf(fpt, "%Lf", &OShift[i * nx + j]);
        }
        fscanf(fpt, "%*[^\n]%*c");
      }
      for (j = 0; j < nx; j++) {
        fscanf(fpt, "%Lf", &OShift[(cf_nums[func_num] - 1) * nx + j]);
      }
    }
    fclose(fpt);

    n_flag = nx;
    func_flag = func_num;
    ini_flag = 1;
    // printf("Function has been initialized!\n");
  }

  for (i = 0; i < mx; i++) {
    switch (func_num) {
    case 1:
      twopeaks_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      f[i] += 100.0;
      break;
    case 2:
      fiveuneven_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      f[i] += 200.0;
      break;
    case 3:
      equalmin_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      f[i] += 300.0;
      break;
    case 4:
      decreasemin_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      f[i] += 400.0;
      break;
    case 5:
      unevenmin_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      f[i] += 500.0;
      break;
    case 6:
      himmelblau_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      f[i] += 600.0;
      break;
    case 7:
      camelback_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      f[i] += 700.0;
      break;
    case 8:
      vincent_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      f[i] += 800.0;
      break;
    case 9:
      cf01(&x[i * nx], &f[i], nx, OShift, M, 1);
      f[i] += 900.0;
      break;
    case 10:
      cf02(&x[i * nx], &f[i], nx, OShift, M, 1);
      f[i] += 1000.0;
      break;
    case 11:
      cf03(&x[i * nx], &f[i], nx, OShift, M, 1);
      f[i] += 1100.0;
      break;
    case 12:
      cf04(&x[i * nx], &f[i], nx, OShift, M, 1);
      f[i] += 1200.0;
      break;
    case 13:
      cf05(&x[i * nx], &f[i], nx, OShift, M, 1);
      f[i] += 1300.0;
      break;
    case 14:
      cf06(&x[i * nx], &f[i], nx, OShift, M, 1);
      f[i] += 1400.0;
      break;
    case 15:
      cf07(&x[i * nx], &f[i], nx, OShift, M, 1);
      f[i] += 1500.0;
      break;

    default:
      printf("\nError: There are only 15 test functions in this test suite!\n");
      f[i] = 0.0;
      break;
    }
  }
}

void twopeaks_func(double *x, double *f, int nx, double *Os, double *Mr,
                   int s_flag, int r_flag) /* Expanded Two-Peak Trap */
{
  int i;
  f[0] = 0.0;
  sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag); /* shift and rotate */
  for (i = 0; i < nx; i++) {
    z[i] += 20.0; // shift to orgin
    if ((z[i] < 15.0) & (z[i] >= 0.0)) {
      f[0] += -(160.0 / 15.0) * (15.0 - z[i]);
    } else if ((z[i] <= 20.0) & (z[i] >= 15.0)) {
      f[0] += -40.0 * (z[i] - 15.0);
    } else if (z[i] < 0.0) {
      f[0] += -160.0 + pow(z[i], 2.0);
    } else {
      f[0] += -200.0 + pow(z[i] - 20.0, 2.0);
    }
  }
  f[0] += 200.0 * nx;
}

void fiveuneven_func(double *x, double *f, int nx, double *Os, double *Mr,
                     int s_flag,
                     int r_flag) /* Expanded Five-Uneven-Peak Trap */
{
  int i;
  f[0] = 0.0;
  sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag); /* shift and rotate */
  for (i = 0; i < nx; i++) {
    if (z[i] < 0) {
      f[0] += -200.0 + pow(z[i], 2.0);
    } else if (z[i] < 2.5) {
      f[0] += -80.0 * (2.5 - z[i]);
    } else if (z[i] < 5.0) {
      f[0] += -64.0 * (z[i] - 2.5);
    } else if (z[i] < 7.5) {
      f[0] += -160.0 + pow(z[i], 2.0);
    } else if (z[i] < 12.5) {
      f[0] += -28.0 * (z[i] - 7.5);
    } else if (z[i] < 17.5) {
      f[0] += -28.0 * (17.5 - z[i]);
    } else if (z[i] < 22.5) {
      f[0] += -32.0 * (z[i] - 17.5);
    } else if (z[i] < 27.5) {
      f[0] += -32.0 * (27.5 - z[i]);
    } else if (z[i] <= 30.0) {
      f[0] += -80.0 * (z[i] - 27.5);
    } else {
      f[0] += -200.0 + pow(z[i] - 30.0, 2.0);
    }
  }
  f[0] += 200.0 * nx;
}

void equalmin_func(double *x, double *f, int nx, double *Os, double *Mr,
                   int s_flag, int r_flag) /* Expanded Equal Minima */
{
  int i;
  f[0] = 0.0;
  sr_func(x, z, nx, Os, Mr, 1.0 / 20.0, s_flag, r_flag); /* shift and rotate */
  for (i = 0; i < nx; i++) {
    z[i] += 0.1; // shift to orgin
    if ((z[i] <= 1.0) & (z[i] >= 0.0)) {
      f[0] += -pow((sin(5 * PI * z[i])), 6.0);
    } else {
      f[0] += pow(z[i], 2.0);
    }
  }
  f[0] += 1.0 * nx;
}

void decreasemin_func(double *x, double *f, int nx, double *Os, double *Mr,
                      int s_flag, int r_flag) /* Expanded Decreasing Minima */
{
  int i;
  f[0] = 0.0;
  sr_func(x, z, nx, Os, Mr, 1.0 / 20.0, s_flag, r_flag); /* shift and rotate */
  for (i = 0; i < nx; i++) {
    z[i] += 0.1; // shift to orgin
    if ((z[i] <= 1.0) & (z[i] >= 0.0)) {
      f[0] += -exp(-2.0 * log(2.0) * pow((z[i] - 0.1) / 0.8, 2.0)) *
              pow(sin(5.0 * PI * z[i]), 6.0);
    } else {
      f[0] += pow(z[i], 2.0);
    }
  }
  f[0] += 1.0 * nx;
}

void unevenmin_func(double *x, double *f, int nx, double *Os, double *Mr,
                    int s_flag, int r_flag) /* Expanded Uneven Minima */
{
  int i;
  f[0] = 0.0;
  sr_func(x, z, nx, Os, Mr, 1.0 / 20.0, s_flag, r_flag); /* shift and rotate */
  for (i = 0; i < nx; i++) {
    z[i] += 0.079699392688696; // pow(0.15,4.0/3.0);//shift to orgin
    if ((z[i] <= 1.0) & (z[i] >= 0.0)) {
      f[0] -= pow(sin(5.0 * PI * (pow(z[i], 0.75) - 0.05)), 6.0);
    } else {
      f[0] += pow(z[i], 2.0);
    }
  }
  f[0] += 1.0 * nx;
}

void himmelblau_func(double *x, double *f, int nx, double *Os, double *Mr,
                     int s_flag,
                     int r_flag) /* Expanded Himmelblau¡¯s Function */
{
  int i;
  f[0] = 0.0;
  sr_func(x, z, nx, Os, Mr, 1.0 / 5.0, s_flag, r_flag); /* shift and rotate */
  for (i = 0; i < nx - 1; i = i + 2) {
    z[i] += 3.0;
    z[i + 1] += 2.0; // shift to orgin
    f[0] += pow((z[i] * z[i] + z[i + 1] - 11.0), 2.0) +
            pow((z[i] + z[i + 1] * z[i + 1] - 7.0), 2.0);
  }
}

void camelback_func(double *x, double *f, int nx, double *Os, double *Mr,
                    int s_flag, int r_flag) /* Expanded Six-Hump Camel Back */
{
  int i;
  f[0] = 0.0;
  sr_func(x, z, nx, Os, Mr, 1.0 / 20.0, s_flag, r_flag); /* shift and rotate */
  for (i = 0; i < nx - 1; i = i + 2) {
    z[i] += 0.089842;
    z[i + 1] += -0.712656; // shift to orgin
    f[0] +=
        ((4.0 - 2.1 * pow(z[i], 2.0) + pow(z[i], 4.0) / 3.0) * pow(z[i], 2.0) +
         z[i] * z[i + 1] +
         ((-4.0 + 4.0 * pow(z[i + 1], 2.0)) * pow(z[i + 1], 2.0))) *
        4.0;
  }
  f[0] += 4.126514 * nx / 2.0;
}

void vincent_func(double *x, double *f, int nx, double *Os, double *Mr,
                  int s_flag, int r_flag) /* Modified Vincent Function */
{
  // orginal bound [0.25, 10], optima=[0.333;
  // 0.6242; 1.1701; 2.1933; 4.1112; 7.7063]
  int i;
  f[0] = 0.0;
  sr_func(x, z, nx, Os, Mr, 1.0 / 5.0, s_flag, r_flag); /* shift and rotate */

  for (i = 0; i < nx; i++) {
    z[i] += 4.1112; // shift to orgin
    if ((z[i] >= 0.25) & (z[i] <= 10.0)) {
      f[0] += -sin(10.0 * logf(z[i]));
    } else if (z[i] < 0.25) {
      f[0] += pow(0.25 - z[i], 2.0) - sin(10.0 * logf(2.5));
    } else {
      f[0] += pow(z[i] - 10, 2.0) - sin(10.0 * logf(10.0));
    }
  }
  f[0] = f[0] / nx;
  f[0] += 1.0;
}

void sphere_func(double *x, double *f, int nx, double *Os, double *Mr,
                 int s_flag, int r_flag) /* Sphere */
{
  int i;
  f[0] = 0.0;
  sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag); /* shift and rotate */
  for (i = 0; i < nx; i++) {
    f[0] += z[i] * z[i];
  }
}

void ellips_func(double *x, double *f, int nx, double *Os, double *Mr,
                 int s_flag, int r_flag) /* Ellipsoidal */
{
  int i;
  f[0] = 0.0;
  sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag); /* shift and rotate */
  for (i = 0; i < nx; i++) {
    f[0] += pow(10.0, 6.0 * i / (nx - 1)) * z[i] * z[i];
  }
}

void bent_cigar_func(double *x, double *f, int nx, double *Os, double *Mr,
                     int s_flag, int r_flag) /* Bent_Cigar */
{
  int i;
  sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag); /* shift and rotate */
  f[0] = z[0] * z[0];
  for (i = 1; i < nx; i++) {
    f[0] += pow(10.0, 6.0) * z[i] * z[i];
  }
}

void discus_func(double *x, double *f, int nx, double *Os, double *Mr,
                 int s_flag, int r_flag) /* Discus */
{
  int i;
  sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag); /* shift and rotate */
  f[0] = pow(10.0, 6.0) * z[0] * z[0];
  for (i = 1; i < nx; i++) {
    f[0] += z[i] * z[i];
  }
}

void dif_powers_func(double *x, double *f, int nx, double *Os, double *Mr,
                     int s_flag, int r_flag) /* Different Powers */
{
  int i;
  f[0] = 0.0;
  sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag); /* shift and rotate */

  for (i = 0; i < nx; i++) {
    f[0] += pow(fabs(z[i]), 2.0 + 4.0 * i / (nx - 1.0));
  }
  f[0] = pow(f[0], 0.5);
}

void rosenbrock_func(double *x, double *f, int nx, double *Os, double *Mr,
                     int s_flag, int r_flag) /* Rosenbrock's */
{
  int i;
  double tmp1, tmp2;
  f[0] = 0.0;
  sr_func(x, z, nx, Os, Mr, 2.048 / 100.0, s_flag,
          r_flag); /* shift and rotate */
  z[0] += 1.0;     // shift to orgin
  for (i = 0; i < nx - 1; i++) {
    z[i + 1] += 1.0; // shift to orgin
    tmp1 = z[i] * z[i] - z[i + 1];
    tmp2 = z[i] - 1.0;
    f[0] += 100.0 * tmp1 * tmp1 + tmp2 * tmp2;
  }
}

void schaffer_F7_func(double *x, double *f, int nx, double *Os, double *Mr,
                      int s_flag, int r_flag) /* Schwefel's 1.2  */
{
  int i;
  double tmp;
  f[0] = 0.0;
  sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag); /* shift and rotate */
  for (i = 0; i < nx - 1; i++) {
    z[i] = pow(y[i] * y[i] + y[i + 1] * y[i + 1], 0.5);
    tmp = sin(50.0 * pow(z[i], 0.2));
    f[0] += pow(z[i], 0.5) + pow(z[i], 0.5) * tmp * tmp;
  }
  f[0] = f[0] * f[0] / (nx - 1.0) / (nx - 1.0);
}

void ackley_func(double *x, double *f, int nx, double *Os, double *Mr,
                 int s_flag, int r_flag) /* Ackley's  */
{
  int i;
  double sum1, sum2;
  sum1 = 0.0;
  sum2 = 0.0;

  sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag); /* shift and rotate */

  for (i = 0; i < nx; i++) {
    sum1 += z[i] * z[i];
    sum2 += cos(2.0 * PI * z[i]);
  }
  sum1 = -0.2 * sqrt(sum1 / nx);
  sum2 /= nx;
  f[0] = E - 20.0 * exp(sum1) - exp(sum2) + 20.0;
}

void weierstrass_func(double *x, double *f, int nx, double *Os, double *Mr,
                      int s_flag, int r_flag) /* Weierstrass's  */
{
  int i, j, k_max;
  double sum, sum2, a, b;
  a = 0.5;
  b = 3.0;
  k_max = 20;
  f[0] = 0.0;

  sr_func(x, z, nx, Os, Mr, 0.5 / 100.0, s_flag, r_flag); /* shift and rotate */

  for (i = 0; i < nx; i++) {
    sum = 0.0;
    sum2 = 0.0;
    for (j = 0; j <= k_max; j++) {
      sum += pow(a, j) * cos(2.0 * PI * pow(b, j) * (z[i] + 0.5));
      sum2 += pow(a, j) * cos(2.0 * PI * pow(b, j) * 0.5);
    }
    f[0] += sum;
  }
  f[0] -= nx * sum2;
}

void griewank_func(double *x, double *f, int nx, double *Os, double *Mr,
                   int s_flag, int r_flag) /* Griewank's  */
{
  int i;
  double s, p;
  s = 0.0;
  p = 1.0;

  sr_func(x, z, nx, Os, Mr, 600.0 / 100.0, s_flag,
          r_flag); /* shift and rotate */

  for (i = 0; i < nx; i++) {
    s += z[i] * z[i];
    p *= cos(z[i] / sqrt(1.0 + i));
  }
  f[0] = 1.0 + s / 4000.0 - p;
}

void rastrigin_func(double *x, double *f, int nx, double *Os, double *Mr,
                    int s_flag, int r_flag) /* Rastrigin's  */
{
  int i;
  f[0] = 0.0;

  sr_func(x, z, nx, Os, Mr, 5.12 / 100.0, s_flag,
          r_flag); /* shift and rotate */

  for (i = 0; i < nx; i++) {
    f[0] += (z[i] * z[i] - 10.0 * cos(2.0 * PI * z[i]) + 10.0);
  }
}

void step_rastrigin_func(double *x, double *f, int nx, double *Os, double *Mr,
                         int s_flag,
                         int r_flag) /* Noncontinuous Rastrigin's  */
{
  int i;
  f[0] = 0.0;
  for (i = 0; i < nx; i++) {
    if (fabs(y[i] - Os[i]) > 0.5)
      y[i] = Os[i] + floor(2.0 * (y[i] - Os[i]) + 0.5) / 2.0;
  }

  sr_func(x, z, nx, Os, Mr, 5.12 / 100.0, s_flag,
          r_flag); /* shift and rotate */

  for (i = 0; i < nx; i++) {
    f[0] += (z[i] * z[i] - 10.0 * cos(2.0 * PI * z[i]) + 10.0);
  }
}

void schwefel_func(double *x, double *f, int nx, double *Os, double *Mr,
                   int s_flag, int r_flag) /* Schwefel's  */
{
  int i;
  double tmp;
  f[0] = 0.0;

  sr_func(x, z, nx, Os, Mr, 1000.0 / 100.0, s_flag,
          r_flag); /* shift and rotate */

  for (i = 0; i < nx; i++) {
    z[i] += 4.209687462275036e+002;
    if (z[i] > 500) {
      f[0] -= (500.0 - fmod(z[i], 500.0)) *
              sin(pow(500.0 - fmod(z[i], 500.0), 0.5));
      tmp = (z[i] - 500.0) / 100.0;
      f[0] += tmp * tmp / nx;
    } else if (z[i] < -500) {
      f[0] -= (-500.0 + fmod(fabs(z[i]), 500.0)) *
              sin(pow(500.0 - fmod(fabs(z[i]), 500.0), 0.5));
      tmp = (z[i] + 500.0) / 100.0;
      f[0] += tmp * tmp / nx;
    } else
      f[0] -= z[i] * sin(pow(fabs(z[i]), 0.5));
  }
  f[0] += 4.189828872724338e+002 * nx;
}

void katsuura_func(double *x, double *f, int nx, double *Os, double *Mr,
                   int s_flag, int r_flag) /* Katsuura  */
{
  int i, j;
  double temp, tmp1, tmp2, tmp3;
  f[0] = 1.0;
  tmp3 = pow(1.0 * nx, 1.2);

  sr_func(x, z, nx, Os, Mr, 5.0 / 100.0, s_flag, r_flag); /* shift and rotate */

  for (i = 0; i < nx; i++) {
    temp = 0.0;
    for (j = 1; j <= 32; j++) {
      tmp1 = pow(2.0, j);
      tmp2 = tmp1 * z[i];
      temp += fabs(tmp2 - floor(tmp2 + 0.5)) / tmp1;
    }
    f[0] *= pow(1.0 + (i + 1.0) * temp, 10.0 / tmp3);
  }
  tmp1 = 10.0 / nx / nx;
  f[0] = f[0] * tmp1 - tmp1;
}

void bi_rastrigin_func(double *x, double *f, int nx, double *Os, double *Mr,
                       int s_flag,
                       int r_flag) /* Lunacek Bi_rastrigin Function */
{
  int i;
  double mu0 = 2.5, d = 1.0, s, mu1, tmp, tmp1, tmp2;
  double *tmpx;
  tmpx = (double *)malloc(sizeof(double) * nx);
  s = 1.0 - 1.0 / (2.0 * pow(nx + 20.0, 0.5) - 8.2);
  mu1 = -pow((mu0 * mu0 - d) / s, 0.5);

  if (s_flag == 1)
    shiftfunc(x, y, nx, Os);
  else {
    for (i = 0; i < nx; i++) // shrink to the orginal search range
    {
      y[i] = x[i];
    }
  }
  for (i = 0; i < nx; i++) // shrink to the orginal search range
  {
    y[i] *= 10.0 / 100.0;
  }

  for (i = 0; i < nx; i++) {
    tmpx[i] = 2 * y[i];
    if (Os[i] < 0.0)
      tmpx[i] *= -1.0;
  }
  for (i = 0; i < nx; i++) {
    z[i] = tmpx[i];
    tmpx[i] += mu0;
  }
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
  tmp = 0.0;

  if (r_flag == 1) {
    rotatefunc(z, y, nx, Mr);
    for (i = 0; i < nx; i++) {
      tmp += cos(2.0 * PI * y[i]);
    }
    if (tmp1 < tmp2)
      f[0] = tmp1;
    else
      f[0] = tmp2;
    f[0] += 10.0 * (nx - tmp);
  } else {
    for (i = 0; i < nx; i++) {
      tmp += cos(2.0 * PI * z[i]);
    }
    if (tmp1 < tmp2)
      f[0] = tmp1;
    else
      f[0] = tmp2;
    f[0] += 10.0 * (nx - tmp);
  }

  free(tmpx);
}

void grie_rosen_func(double *x, double *f, int nx, double *Os, double *Mr,
                     int s_flag, int r_flag) /* Griewank-Rosenbrock  */
{
  int i;
  double temp, tmp1, tmp2;
  f[0] = 0.0;

  sr_func(x, z, nx, Os, Mr, 5.0 / 100.0, s_flag, r_flag); /* shift and rotate */

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
}

void escaffer6_func(double *x, double *f, int nx, double *Os, double *Mr,
                    int s_flag, int r_flag) /* Expanded Scaffer¡¯s F6  */
{
  int i;
  double temp1, temp2;

  sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag); /* shift and rotate */

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

void happycat_func(
    double *x, double *f, int nx, double *Os, double *Mr, int s_flag,
    int r_flag) /* HappyCat, provdided by Hans-Georg Beyer (HGB) */
/* original global optimum: [-1,-1,...,-1] */
{
  int i;
  double alpha, r2, sum_z;
  alpha = 1.0 / 8.0;

  sr_func(x, z, nx, Os, Mr, 5.0 / 100.0, s_flag, r_flag); /* shift and rotate */

  r2 = 0.0;
  sum_z = 0.0;
  for (i = 0; i < nx; i++) {
    z[i] = z[i] - 1.0; // shift to orgin
    r2 += z[i] * z[i];
    sum_z += z[i];
  }
  f[0] = pow(fabs(r2 - nx), 2.0 * alpha) + (0.5 * r2 + sum_z) / nx + 0.5;
}

void hgbat_func(double *x, double *f, int nx, double *Os, double *Mr,
                int s_flag,
                int r_flag) /* HGBat, provdided by Hans-Georg Beyer (HGB)*/
/* original global optimum: [-1,-1,...,-1] */
{
  int i;
  double alpha, r2, sum_z;
  alpha = 1.0 / 4.0;

  sr_func(x, z, nx, Os, Mr, 5.0 / 100.0, s_flag, r_flag); /* shift and rotate */

  r2 = 0.0;
  sum_z = 0.0;
  for (i = 0; i < nx; i++) {
    z[i] = z[i] - 1.0; // shift to orgin
    r2 += z[i] * z[i];
    sum_z += z[i];
  }
  f[0] = pow(fabs(pow(r2, 2.0) - pow(sum_z, 2.0)), 2 * alpha) +
         (0.5 * r2 + sum_z) / nx + 0.5;
}

void cf01(double *x, double *f, int nx, double *Os, double *Mr,
          int r_flag) /* Composition Function 1 */
{
  int i, cf_num = 10;
  double fit[10];
  double delta[10] = {10, 20, 10, 20, 10, 20, 10, 20, 10, 20};
  double bias[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  i = 0;
  sphere_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  i = 1;
  sphere_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  i = 2;
  ellips_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = fit[i] / 1e+6;
  i = 3;
  ellips_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = fit[i] / 1e+6;
  i = 4;
  bent_cigar_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = fit[i] / 1e+6;
  i = 5;
  bent_cigar_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = fit[i] / 1e+6;
  i = 6;
  discus_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = fit[i] / 1e+4;
  i = 7;
  discus_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = fit[i] / 1e+4;
  i = 8;
  dif_powers_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = fit[i] / 1e+5;
  i = 9;
  dif_powers_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = fit[i] / 1e+5;
  cf_cal(x, f, nx, Os, delta, bias, fit, cf_num);
}

void cf02(double *x, double *f, int nx, double *Os, double *Mr,
          int r_flag) /* Composition Function 3 */
{
  int i, cf_num = 10;
  double fit[10];
  double delta[10] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
  double bias[10] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90};
  i = 0;
  ellips_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = fit[i] / 1e+6;
  i = 1;
  ellips_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = fit[i] / 1e+6;
  i = 2;
  dif_powers_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = fit[i] / 1e+5;
  i = 3;
  dif_powers_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = fit[i] / 1e+5;
  i = 4;
  bent_cigar_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = fit[i] / 1e+6;
  i = 5;
  bent_cigar_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = fit[i] / 1e+6;
  i = 6;
  discus_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = fit[i] / 1e+4;
  i = 7;
  discus_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = fit[i] / 1e+4;
  i = 8;
  sphere_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  i = 9;
  sphere_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  cf_cal(x, f, nx, Os, delta, bias, fit, cf_num);
}

void cf03(double *x, double *f, int nx, double *Os, double *Mr,
          int r_flag) /* Composition Function 4 */
{
  int i, cf_num = 10;
  double fit[10];
  double delta[10] = {10, 10, 10, 10, 10, 10, 10, 10, 10, 10};
  double bias[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  ;
  i = 0;
  rosenbrock_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = fit[i] / 10;
  i = 1;
  rosenbrock_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = fit[i] / 10;
  i = 2;
  rastrigin_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = fit[i] * 10;
  i = 3;
  rastrigin_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = fit[i] * 10;
  i = 4;
  happycat_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = fit[i] * 10;
  i = 5;
  happycat_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = fit[i] * 10;
  i = 6;
  escaffer6_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = fit[i] * 100;
  i = 7;
  escaffer6_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = fit[i] * 100;
  i = 8;
  schwefel_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  i = 9;
  schwefel_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  cf_cal(x, f, nx, Os, delta, bias, fit, cf_num);
}

void cf04(double *x, double *f, int nx, double *Os, double *Mr,
          int r_flag) /* Composition Function 5 */
{
  int i, cf_num = 10;
  double fit[10];
  double delta[10] = {10, 10, 20, 20, 30, 30, 40, 40, 50, 50};
  double bias[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  i = 0;
  rosenbrock_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = fit[i] / 10;
  i = 1;
  rosenbrock_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = fit[i] / 10;
  i = 2;
  rastrigin_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = fit[i] * 10;
  i = 3;
  rastrigin_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = fit[i] * 10;
  i = 4;
  happycat_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = fit[i] * 10;
  i = 5;
  happycat_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = fit[i] * 10;
  i = 6;
  escaffer6_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = fit[i] * 100;
  i = 7;
  escaffer6_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = fit[i] * 100;
  i = 8;
  schwefel_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  i = 9;
  schwefel_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  cf_cal(x, f, nx, Os, delta, bias, fit, cf_num);
}

void cf05(double *x, double *f, int nx, double *Os, double *Mr,
          int r_flag) /* Composition Function 6 */
{
  int i, cf_num = 10;
  double fit[10];
  double delta[10] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
  double bias[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  i = 0;
  rosenbrock_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = fit[i] / 10.0;
  i = 1;
  hgbat_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = fit[i] * 10.0;
  i = 2;
  rastrigin_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = fit[i] * 10.0;
  i = 3;
  ackley_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = fit[i] / 10.0;
  i = 4;
  weierstrass_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = fit[i] * 2.5;
  i = 5;
  katsuura_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = fit[i] / 1e+3;
  i = 6;
  escaffer6_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = fit[i] * 100;
  i = 7;
  grie_rosen_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = fit[i] * 2.5;
  i = 8;
  happycat_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = fit[i] * 10.0;
  i = 9;
  schwefel_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  cf_cal(x, f, nx, Os, delta, bias, fit, cf_num);
}

void cf06(double *x, double *f, int nx, double *Os, double *Mr,
          int r_flag) /* Composition Function 7 */
{
  int i, cf_num = 10;
  double fit[10];
  double delta[10] = {10, 10, 20, 20, 30, 30, 40, 40, 50, 50};
  double bias[10] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180};
  i = 0;
  rastrigin_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = fit[i] * 10.0;
  i = 1;
  schwefel_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  i = 2;
  rastrigin_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = fit[i] * 10.0;
  i = 3;
  schwefel_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  i = 4;
  rastrigin_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = fit[i] * 10.0;
  i = 5;
  schwefel_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  i = 6;
  rastrigin_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = fit[i] * 10.0;
  i = 7;
  schwefel_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  i = 8;
  rastrigin_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = fit[i] * 10.0;
  i = 9;
  schwefel_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  cf_cal(x, f, nx, Os, delta, bias, fit, cf_num);
}

void cf07(double *x, double *f, int nx, double *Os, double *Mr,
          int r_flag) /* Composition Function 8 */
{
  int i, cf_num = 10;
  double fit[10];
  double delta[10] = {10, 10, 20, 20, 30, 30, 40, 40, 50, 50};
  double bias[10] = {0, 20, 40, 60, 80, 100, 120, 140, 160, 180};
  i = 0;
  rosenbrock_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = fit[i] / 10.0;
  i = 1;
  hgbat_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = fit[i] * 10.0;
  i = 2;
  rastrigin_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = fit[i] * 10.0;
  i = 3;
  ackley_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = fit[i] / 10.0;
  i = 4;
  weierstrass_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = fit[i] * 2.5;
  i = 5;
  katsuura_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = fit[i] / 1e+3;
  i = 6;
  escaffer6_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = fit[i] * 100;
  i = 7;
  grie_rosen_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = fit[i] * 2.5;
  i = 8;
  happycat_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  fit[i] = fit[i] * 10.0;
  i = 9;
  schwefel_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
  cf_cal(x, f, nx, Os, delta, bias, fit, cf_num);
}

void shiftfunc(double *x, double *xshift, int nx, double *Os) {
  int i;
  for (i = 0; i < nx; i++) {
    xshift[i] = x[i] - Os[i];
  }
}

void rotatefunc(double *x, double *xrot, int nx, double *Mr) {
  int i, j;
  for (i = 0; i < nx; i++) {
    xrot[i] = 0;
    for (j = 0; j < nx; j++) {
      xrot[i] = xrot[i] + x[j] * Mr[i * nx + j];
    }
  }
}

void sr_func(double *x, double *sr_x, int nx, double *Os, double *Mr,
             double sh_rate, int s_flag, int r_flag) /* shift and rotate */
{
  int i;
  if (s_flag == 1) {
    if (r_flag == 1) {
      shiftfunc(x, y, nx, Os);
      for (i = 0; i < nx; i++) // shrink to the orginal search range
      {
        y[i] = y[i] * sh_rate;
      }
      rotatefunc(y, sr_x, nx, Mr);
    } else {
      shiftfunc(x, sr_x, nx, Os);
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
      rotatefunc(y, sr_x, nx, Mr);
    } else
      for (i = 0; i < nx; i++) {
        for (i = 0; i < nx; i++) // shrink to the orginal search range
        {
          sr_x[i] = x[i] * sh_rate;
        }
      }
  }
}

void cf_cal(double *x, double *f, int nx, double *Os, double *delta,
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
