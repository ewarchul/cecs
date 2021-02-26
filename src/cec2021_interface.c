#include "cec2021/cec2021_interface.h"
#include <string.h>

void cec2021_func(double *x, double *f, int nx, int mx, int func_num,
                  char **suite) {

  char *suitex = *suite;

  if (!strcmp(suitex, "basic")) {

    cec2021_basic_func(x, f, nx, mx, func_num);
  } else if (!strcmp(suitex, "shift")) {

    cec2021_shift_func(x, f, nx, mx, func_num);
  } else if (!strcmp(suitex, "rot")) {

    cec2021_rot_func(x, f, nx, mx, func_num);
  } else if (!strcmp(suitex, "bias")) {

    cec2021_bias_func(x, f, nx, mx, func_num);
  } else if (!strcmp(suitex, "shift_rot")) {

    cec2021_shift_rot_func(x, f, nx, mx, func_num);
  } else if (!strcmp(suitex, "bias_rot")) {

    cec2021_bias_rot_func(x, f, nx, mx, func_num);
  } else if (!strcmp(suitex, "bias_shift")) {

    cec2021_bias_shift_func(x, f, nx, mx, func_num);
  } else if (!strcmp(suitex, "bias_shift_rot")) {

    cec2021_bias_shift_rot_func(x, f, nx, mx, func_num);
  } else {
    return;
  }
}

void cec2021_basic_func(double *x, double *f, int nx, int mx, int func_num) {

  int cf_num = 10, i, j;

  if (func_num < 1 || func_num > 10) {
    perror("\nError: Test function is not defined.\n");
  }
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

    if (!(nx == 2 || nx == 10 || nx == 20)) {
      perror("\nError: Test functions are only defined for D=2,10,20.\n");
    }
    if (nx == 2 && (func_num == 5 || func_num == 6 || func_num == 7)) {
      perror("\nError:  NOT defined for D=2.\n");
    }

    /* Load Matrix M*/
    sprintf(FileName, "%s/M_%d_D%d_nr.txt", extdata, func_num, nx);
    fpt = fopen(FileName, "r");
    if (fpt == NULL) {
      perror("\n Error: Cannot open M_d_Dd_nr.txt's for reading \n");
    }
    if (func_num < 7) {
      M = (double *)malloc(nx * nx * sizeof(double));
      if (M == NULL)
        perror("\nError: there is insufficient memory available!\n");
      for (i = 0; i < nx * nx; i++) {
        if (!fscanf(fpt, "%lf", &M[i])) {
          perror("\nError\n");
        }
      }
    } else {
      M = (double *)malloc(cf_num * nx * nx * sizeof(double));
      if (M == NULL)
        perror("\nError: there is insufficient memory available!\n");
      for (i = 0; i < cf_num * nx * nx; i++) {
        if (!fscanf(fpt, "%lf", &M[i])) {
          perror("\nError\n");
        }
      }
    }
    fclose(fpt);

    /* Load shift_data */
    sprintf(FileName, "%s/shift_data_%d_ns.txt", extdata, func_num);
    fpt = fopen(FileName, "r");
    if (fpt == NULL) {
      perror("\n Error: Cannot open shift_data_d_ns.txt for reading \n");
    }

    if (func_num < 7) {
      OShift = (double *)malloc(nx * sizeof(double));
      if (OShift == NULL)
        perror("\nError: there is insufficient memory available!\n");
      for (i = 0; i < nx; i++) {
        if (!fscanf(fpt, "%lf", &OShift[i])) {
          perror("\nError\n");
        }
      }
    } else {
      // OShift=(double *)malloc(nx*sizeof(double));
      OShift = (double *)malloc(nx * cf_num * sizeof(double));
      if (OShift == NULL)
        perror("\nError: there is insufficient memory available!\n");
      for (i = 0; i < cf_num - 1; i++) {
        for (j = 0; j < nx; j++) {
          if (!fscanf(fpt, "%lf", &OShift[i * nx + j])) {

            perror("\nError\n");
          }
        }
      //  if (fscanf(fpt, "%*[^\n]%*c") != 1) {
        //  perror("\nError\n");
     //   }
      }
      for (j = 0; j < nx; j++) {
        if (!fscanf(fpt, "%lf", &OShift[nx * (cf_num - 1) + j])) {
          perror("\nError\n");
        }
      }
    }
    fclose(fpt);

    /* Load Shuffle_data */

    if (func_num >= 5 && func_num <= 7) {
      sprintf(FileName, "%s/shuffle_data_%d_D%d.txt", extdata, func_num, nx);
      fpt = fopen(FileName, "r");
      if (fpt == NULL) {
        perror("\n Error: Cannot open shuffle_data_d_Dd.txt's for reading \n");
      }
      SS = (int *)malloc(nx * sizeof(int));
      if (SS == NULL)
        perror("\nError: there is insufficient memory available!\n");
      for (i = 0; i < nx; i++) {
        if (!fscanf(fpt, "%d", &SS[i])) {
          perror("\nError\n");
        }
      }
      fclose(fpt);
    }

    n_flag = nx;
    func_flag = func_num;
    ini_flag = 1;
    // perror("Function has been initialized!\n");
  }

  for (i = 0; i < mx; i++) {
    switch (func_num) {
    case 1:
      cec2021_bent_cigar_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      break;
    case 2:
      cec2021_schwefel_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      break;
    case 3:
      cec2021_bi_rastrigin_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      break;
    case 4:
      cec2021_grie_rosen_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      break;
    case 5:
      cec2021_hf01(&x[i * nx], &f[i], nx, OShift, M, SS, 1, 1);
      break;
    case 6:
      cec2021_hf06(&x[i * nx], &f[i], nx, OShift, M, SS, 1, 1);
      break;
    case 7:
      cec2021_hf05(&x[i * nx], &f[i], nx, OShift, M, SS, 1, 1);
      break;
    case 8:
      cec2021_cf02(&x[i * nx], &f[i], nx, OShift, M, 1);
      break;
    case 9:
      cec2021_cf04(&x[i * nx], &f[i], nx, OShift, M, 1);
      break;
    case 10:
      cec2021_cf05(&x[i * nx], &f[i], nx, OShift, M, 1);
      break;
    default:
      perror("\nError: There are only 10 test functions in this test suite!\n");
      break;
    }
  }
}

void cec2021_bias_rot_func(double *x, double *f, int nx, int mx, int func_num) {

  int cf_num = 10, i, j;
  if (func_num < 1 || func_num > 10) {
    perror("\nError: Test function is not defined.\n");
  }
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

    if (!(nx == 2 || nx == 10 || nx == 20)) {
      perror("\nError: Test functions are only defined for D=2,10,20.\n");
    }
    if (nx == 2 && (func_num == 5 || func_num == 6 || func_num == 7)) {
      perror("\nError:  NOT defined for D=2.\n");
    }

    /* Load Matrix M*/
    sprintf(FileName, "%s/M_%d_D%d.txt", extdata, func_num, nx);
    fpt = fopen(FileName, "r");
    if (fpt == NULL) {
      perror("\n Error: Cannot open M_d_Dd.txt for reading \n");
    }
    if (func_num < 7) {
      M = (double *)malloc(nx * nx * sizeof(double));
      if (M == NULL)
        perror("\nError: there is insufficient memory available!\n");
      for (i = 0; i < nx * nx; i++) {
        if (!fscanf(fpt, "%lf", &M[i])) {

          perror("\nError\n");
        }
      }
    } else {
      M = (double *)malloc(cf_num * nx * nx * sizeof(double));
      if (M == NULL)
        perror("\nError: there is insufficient memory available!\n");
      for (i = 0; i < cf_num * nx * nx; i++) {
        if (!fscanf(fpt, "%lf", &M[i])) {
          perror("\nError\n");
        }
      }
    }
    fclose(fpt);

    /* Load shift_data */
    sprintf(FileName, "%s/shift_data_%d_ns.txt", extdata, func_num);
    fpt = fopen(FileName, "r");
    if (fpt == NULL) {
      perror("\n Error: Cannot open shift_data_d_ns.txt for reading \n");
    }

    if (func_num < 7) {
      OShift = (double *)malloc(nx * sizeof(double));
      if (OShift == NULL)
        perror("\nError: there is insufficient memory available!\n");
      for (i = 0; i < nx; i++) {
        if (!fscanf(fpt, "%lf", &OShift[i])) {

          perror("\nError\n");
        }
      }
    } else {
      // OShift=(double *)malloc(nx*sizeof(double));
      OShift = (double *)malloc(nx * cf_num * sizeof(double));
      if (OShift == NULL)
        perror("\nError: there is insufficient memory available!\n");
      for (i = 0; i < cf_num - 1; i++) {
        for (j = 0; j < nx; j++) {
          if (!fscanf(fpt, "%lf", &OShift[i * nx + j])) {

            perror("\nError\n");
          }
        }
       // if (fscanf(fpt, "%*[^\n]%*c") != 1) {

         // perror("\nError\n");
        //}
      }
      for (j = 0; j < nx; j++) {
        if (!fscanf(fpt, "%lf", &OShift[nx * (cf_num - 1) + j])) {

          perror("\nError\n");
        }
      }
    }
    fclose(fpt);

    /* Load Shuffle_data */

    if (func_num >= 5 && func_num <= 7) {
      sprintf(FileName, "%s/shuffle_data_%d_D%d.txt", extdata, func_num, nx);
      fpt = fopen(FileName, "r");
      if (fpt == NULL) {
        perror("\n Error: Cannot open shuffle_data_d_Dd.txt's for reading \n");
      }
      SS = (int *)malloc(nx * sizeof(int));
      if (SS == NULL)
        perror("\nError: there is insufficient memory available!\n");
      for (i = 0; i < nx; i++) {
        if (!fscanf(fpt, "%d", &SS[i])) {

          perror("\nError\n");
        }
      }
      fclose(fpt);
    }

    n_flag = nx;
    func_flag = func_num;
    ini_flag = 1;
    // perror("Function has been initialized!\n");
  }

  for (i = 0; i < mx; i++) {
    switch (func_num) {
    case 1:
      cec2021_bent_cigar_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      f[i] += 100.0;
      break;
    case 2:
      cec2021_schwefel_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      f[i] += 1100.0;
      break;
    case 3:
      cec2021_bi_rastrigin_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      f[i] += 700.0;
      break;
    case 4:
      cec2021_grie_rosen_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      f[i] += 1900.0;
      break;
    case 5:
      cec2021_hf01(&x[i * nx], &f[i], nx, OShift, M, SS, 1, 1);
      f[i] = f[i] + 1700.0;
      break;
    case 6:
      cec2021_hf06(&x[i * nx], &f[i], nx, OShift, M, SS, 1, 1);
      f[i] += 1600.0;
      break;
    case 7:
      cec2021_hf05(&x[i * nx], &f[i], nx, OShift, M, SS, 1, 1);
      f[i] += 2100.0;
      break;
    case 8:
      cec2021_cf02(&x[i * nx], &f[i], nx, OShift, M, 1);
      f[i] += 2200.0;
      break;
    case 9:
      cec2021_cf04(&x[i * nx], &f[i], nx, OShift, M, 1);
      f[i] += 2400.0;
      break;
    case 10:
      cec2021_cf05(&x[i * nx], &f[i], nx, OShift, M, 1);
      f[i] += 2500.0;
      break;
    default:
      perror("\nError: There are only 10 test functions in this test suite!\n");
      f[i] = 0.0;
      break;
    }
  }
}

void cec2021_bias_shift_rot_func(double *x, double *f, int nx, int mx,
                                 int func_num) {

  int cf_num = 10, i, j;
  if (func_num < 1 || func_num > 10) {
    perror("\nError: Test function is not defined.\n");
  }
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

    if (!(nx == 2 || nx == 10 || nx == 20)) {
      perror("\nError: Test functions are only defined for D=2,10,20.\n");
    }
    if (nx == 2 && (func_num == 5 || func_num == 6 || func_num == 7)) {
      perror("\nError:  NOT defined for D=2.\n");
    }

    /* Load Matrix M*/
    sprintf(FileName, "%s/M_%d_D%d.txt", extdata, func_num, nx);
    fpt = fopen(FileName, "r");
    if (fpt == NULL) {
      perror("\n Error: Cannot open M_d_Dd.txt's for reading \n");
    }
    if (func_num < 7) {
      M = (double *)malloc(nx * nx * sizeof(double));
      if (M == NULL)
        perror("\nError: there is insufficient memory available!\n");
      for (i = 0; i < nx * nx; i++) {
        if (!fscanf(fpt, "%lf", &M[i])) {

          perror("\nError\n");
        }
      }
    } else {
      M = (double *)malloc(cf_num * nx * nx * sizeof(double));
      if (M == NULL)
        perror("\nError: there is insufficient memory available!\n");
      for (i = 0; i < cf_num * nx * nx; i++) {
        if (!fscanf(fpt, "%lf", &M[i])) {

          perror("\nError\n");
        }
      }
    }
    fclose(fpt);

    /* Load shift_data */
    sprintf(FileName, "%s/shift_data_%d.txt", extdata, func_num);
    fpt = fopen(FileName, "r");
    if (fpt == NULL) {
      perror("\n Error: Cannot open shift_data_d.txt for reading \n");
    }

    if (func_num < 7) {
      OShift = (double *)malloc(nx * sizeof(double));
      if (OShift == NULL)
        perror("\nError: there is insufficient memory available!\n");
      for (i = 0; i < nx; i++) {
        if (!fscanf(fpt, "%lf", &OShift[i])) {

          perror("\nError\n");
        }
      }
    } else {
      // OShift=(double *)malloc(nx*sizeof(double));
      OShift = (double *)malloc(nx * cf_num * sizeof(double));
      if (OShift == NULL)
        perror("\nError: there is insufficient memory available!\n");
      for (i = 0; i < cf_num - 1; i++) {
        for (j = 0; j < nx; j++) {
          if (!fscanf(fpt, "%lf", &OShift[i * nx + j])) {

            perror("\nError\n");
          }
        }
//        if (fscanf(fpt, "%*[^\n]%*c") != 1) {
//
//        perror("\nError\n");
//        }
      }
      for (j = 0; j < nx; j++) {
        if (!fscanf(fpt, "%lf", &OShift[nx * (cf_num - 1) + j])) {

          perror("\nError\n");
        }
      }
    }
    fclose(fpt);

    /* Load Shuffle_data */

    if (func_num >= 5 && func_num <= 7) {
      sprintf(FileName, "%s/shuffle_data_%d_D%d.txt", extdata, func_num, nx);
      fpt = fopen(FileName, "r");
      if (fpt == NULL) {
        perror("\n Error: Cannot open shuffle_data_d_Dd.txt for reading \n");
      }
      SS = (int *)malloc(nx * sizeof(int));
      if (SS == NULL)
        perror("\nError: there is insufficient memory available!\n");
      for (i = 0; i < nx; i++) {
        if (!fscanf(fpt, "%d", &SS[i])) {

          perror("\nError\n");
        }
      }
      fclose(fpt);
    }

    n_flag = nx;
    func_flag = func_num;
    ini_flag = 1;
    // perror("Function has been initialized!\n");
  }

  for (i = 0; i < mx; i++) {
    switch (func_num) {
    case 1:
      cec2021_bent_cigar_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      f[i] += 100.0;
      break;
    case 2:
      cec2021_schwefel_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      f[i] += 1100.0;
      break;
    case 3:
      cec2021_bi_rastrigin_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      f[i] += 700.0;
      break;
    case 4:
      cec2021_grie_rosen_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      f[i] += 1900.0;
      break;
    case 5:
      cec2021_hf01(&x[i * nx], &f[i], nx, OShift, M, SS, 1, 1);
      f[i] = f[i] + 1700.0;
      break;
    case 6:
      cec2021_hf06(&x[i * nx], &f[i], nx, OShift, M, SS, 1, 1);
      f[i] += 1600.0;
      break;
    case 7:
      cec2021_hf05(&x[i * nx], &f[i], nx, OShift, M, SS, 1, 1);
      f[i] += 2100.0;
      break;
    case 8:
      cec2021_cf02_s(&x[i * nx], &f[i], nx, OShift, M, 1);
      f[i] += 2200.0;
      break;
    case 9:
      cec2021_cf04_s(&x[i * nx], &f[i], nx, OShift, M, 1);
      f[i] += 2400.0;
      break;
    case 10:
      cec2021_cf05_s(&x[i * nx], &f[i], nx, OShift, M, 1);
      f[i] += 2500.0;
      break;
    default:
      perror("\nError: There are only 10 test functions in this test suite!\n");
      f[i] = 0.0;
      break;
    }
  }
}

void cec2021_bias_shift_func(double *x, double *f, int nx, int mx,
                             int func_num) {

  int cf_num = 10, i, j;
  if (func_num < 1 || func_num > 10) {
    perror("\nError: Test function is not defined.\n");
  }
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

    if (!(nx == 2 || nx == 10 || nx == 20)) {
      perror("\nError: Test functions are only defined for D=2,10,20.\n");
    }
    if (nx == 2 && (func_num == 5 || func_num == 6 || func_num == 7)) {
      perror("\nError:  NOT defined for D=2.\n");
    }

    /* Load Matrix M*/
    sprintf(FileName, "%s/M_%d_D%d_nr.txt", extdata, func_num, nx);
    fpt = fopen(FileName, "r");
    if (fpt == NULL) {
      perror("\n Error: Cannot open M_d_Dd_nr.txt for reading \n");
    }
    if (func_num < 7) {
      M = (double *)malloc(nx * nx * sizeof(double));
      if (M == NULL)
        perror("\nError: there is insufficient memory available!\n");
      for (i = 0; i < nx * nx; i++) {
        if (!fscanf(fpt, "%lf", &M[i])) {

          perror("\nError\n");
        }
      }
    } else {
      M = (double *)malloc(cf_num * nx * nx * sizeof(double));
      if (M == NULL)
        perror("\nError: there is insufficient memory available!\n");
      for (i = 0; i < cf_num * nx * nx; i++) {
        if (!fscanf(fpt, "%lf", &M[i])) {

          perror("\nError\n");
        }
      }
    }
    fclose(fpt);

    /* Load shift_data */
    sprintf(FileName, "%s/shift_data_%d.txt", extdata, func_num);
    fpt = fopen(FileName, "r");
    if (fpt == NULL) {
      perror("\n Error: Cannot open shift_data_d.txt for reading \n");
    }

    if (func_num < 7) {
      OShift = (double *)malloc(nx * sizeof(double));
      if (OShift == NULL)
        perror("\nError: there is insufficient memory available!\n");
      for (i = 0; i < nx; i++) {
        if (!fscanf(fpt, "%lf", &OShift[i])) {

          perror("\nError\n");
        }
      }
    } else {
      // OShift=(double *)malloc(nx*sizeof(double));
      OShift = (double *)malloc(nx * cf_num * sizeof(double));
      if (OShift == NULL)
        perror("\nError: there is insufficient memory available!\n");
      for (i = 0; i < cf_num - 1; i++) {
        for (j = 0; j < nx; j++) {
          if (!fscanf(fpt, "%lf", &OShift[i * nx + j])) {

            perror("\nError\n");
          }
        }
//        if (fscanf(fpt, "%*[^\n]%*c") != 1) {
//
//          perror("\nError\n");
//        }
      }
      for (j = 0; j < nx; j++) {
        if (!fscanf(fpt, "%lf", &OShift[nx * (cf_num - 1) + j])) {

          perror("\nError\n");
        }
      }
    }
    fclose(fpt);

    /* Load Shuffle_data */

    if (func_num >= 5 && func_num <= 7) {
      sprintf(FileName, "%s/shuffle_data_%d_D%d.txt", extdata, func_num, nx);
      fpt = fopen(FileName, "r");
      if (fpt == NULL) {
        perror("\n Error: Cannot open shuffle_data_d_Dd.txt for reading \n");
      }
      SS = (int *)malloc(nx * sizeof(int));
      if (SS == NULL)
        perror("\nError: there is insufficient memory available!\n");
      for (i = 0; i < nx; i++) {
        if (!fscanf(fpt, "%d", &SS[i])) {
          perror("\nError\n");
        }
      }
      fclose(fpt);
    }

    n_flag = nx;
    func_flag = func_num;
    ini_flag = 1;
    // perror("Function has been initialized!\n");
  }

  for (i = 0; i < mx; i++) {
    switch (func_num) {
    case 1:
      cec2021_bent_cigar_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      f[i] += 100.0;
      break;
    case 2:
      cec2021_schwefel_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      f[i] += 1100.0;
      break;
    case 3:
      cec2021_bi_rastrigin_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      f[i] += 700.0;
      break;
    case 4:
      cec2021_grie_rosen_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      f[i] += 1900.0;
      break;
    case 5:
      cec2021_hf01(&x[i * nx], &f[i], nx, OShift, M, SS, 1, 1);
      f[i] = f[i] + 1700.0;
      break;
    case 6:
      cec2021_hf06(&x[i * nx], &f[i], nx, OShift, M, SS, 1, 1);
      f[i] += 1600.0;
      break;
    case 7:
      cec2021_hf05(&x[i * nx], &f[i], nx, OShift, M, SS, 1, 1);
      f[i] += 2100.0;
      break;
    case 8:
      cec2021_cf02_s(&x[i * nx], &f[i], nx, OShift, M, 1);
      f[i] += 2200.0;
      break;
    case 9:
      cec2021_cf04_s(&x[i * nx], &f[i], nx, OShift, M, 1);
      f[i] += 2400.0;
      break;
    case 10:
      cec2021_cf05_s(&x[i * nx], &f[i], nx, OShift, M, 1);
      f[i] += 2500.0;
      break;
    default:
      perror("\nError: There are only 10 test functions in this test suite!\n");
      f[i] = 0.0;
      break;
    }
  }
}

void cec2021_bias_func(double *x, double *f, int nx, int mx, int func_num) {

  int cf_num = 10, i, j;

  if (func_num < 1 || func_num > 10) {
    perror("\nError: Test function is not defined.\n");
  }

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

    if (!(nx == 2 || nx == 10 || nx == 20)) {
      perror("\nError: Test functions are only defined for D=2,10,20.\n");
    }
    if (nx == 2 && (func_num == 5 || func_num == 6 || func_num == 7)) {
      perror("\nError:  NOT defined for D=2.\n");
    }

    /* Load Matrix M*/
    sprintf(FileName, "%s/M_%d_D%d_nr.txt", extdata, func_num, nx);
    fpt = fopen(FileName, "r");
    if (fpt == NULL) {
      perror("\n Error: Cannot open M_d_Dd_nr.txt for reading \n");
    }
    if (func_num < 7) {
      M = (double *)malloc(nx * nx * sizeof(double));
      if (M == NULL)
        perror("\nError: there is insufficient memory available!\n");
      for (i = 0; i < nx * nx; i++) {
        if (!fscanf(fpt, "%lf", &M[i])) {

          perror("\nError\n");
        }
      }
    } else {
      M = (double *)malloc(cf_num * nx * nx * sizeof(double));
      if (M == NULL)
        perror("\nError: there is insufficient memory available!\n");
      for (i = 0; i < cf_num * nx * nx; i++) {
        if (!fscanf(fpt, "%lf", &M[i])) {

          perror("\nError\n");
        }
      }
    }
    fclose(fpt);

    /* Load shift_data */
    sprintf(FileName, "%s/shift_data_%d_ns.txt", extdata, func_num);
    fpt = fopen(FileName, "r");
    if (fpt == NULL) {
      perror("\n Error: Cannot open shift_data_d_ns.txt for reading \n");
    }

    if (func_num < 7) {
      OShift = (double *)malloc(nx * sizeof(double));
      if (OShift == NULL)
        perror("\nError: there is insufficient memory available!\n");
      for (i = 0; i < nx; i++) {
        if (!fscanf(fpt, "%lf", &OShift[i])) {

          perror("\nError\n");
        }
      }
    } else {
      // OShift=(double *)malloc(nx*sizeof(double));
      OShift = (double *)malloc(nx * cf_num * sizeof(double));
      if (OShift == NULL)
        perror("\nError: there is insufficient memory available!\n");
      for (i = 0; i < cf_num - 1; i++) {
        for (j = 0; j < nx; j++) {
          if (!fscanf(fpt, "%lf", &OShift[i * nx + j])) {

            perror("\nError\n");
          }
        }
//        if (fscanf(fpt, "%*[^\n]%*c") != 1) {
//
//          perror("\nError\n");
//       }
      }
      for (j = 0; j < nx; j++) {
        if (!fscanf(fpt, "%lf", &OShift[nx * (cf_num - 1) + j])) {

          perror("\nError\n");
        }
      }
    }
    fclose(fpt);

    /* Load Shuffle_data */

    if (func_num >= 5 && func_num <= 7) {
      sprintf(FileName, "%s/shuffle_data_%d_D%d.txt", extdata, func_num, nx);
      fpt = fopen(FileName, "r");
      if (fpt == NULL) {
        perror("\n Error: Cannot open shuffle_data_D.txt for reading \n");
      }
      SS = (int *)malloc(nx * sizeof(int));
      if (SS == NULL)
        perror("\nError: there is insufficient memory available!\n");
      for (i = 0; i < nx; i++) {
        if (!fscanf(fpt, "%d", &SS[i])) {

          perror("\nError\n");
        }
      }
      fclose(fpt);
    }

    n_flag = nx;
    func_flag = func_num;
    ini_flag = 1;
    // perror("Function has been initialized!\n");
  }

  for (i = 0; i < mx; i++) {
    switch (func_num) {
    case 1:
      cec2021_bent_cigar_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      f[i] += 100.0;
      break;
    case 2:
      cec2021_schwefel_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      f[i] += 1100.0;
      break;
    case 3:
      cec2021_bi_rastrigin_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      f[i] += 700.0;
      break;
    case 4:
      cec2021_grie_rosen_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      f[i] += 1900.0;
      break;
    case 5:
      cec2021_hf01(&x[i * nx], &f[i], nx, OShift, M, SS, 1, 1);
      f[i] = f[i] + 1700.0;
      break;
    case 6:
      cec2021_hf06(&x[i * nx], &f[i], nx, OShift, M, SS, 1, 1);
      f[i] += 1600.0;
      break;
    case 7:
      cec2021_hf05(&x[i * nx], &f[i], nx, OShift, M, SS, 1, 1);
      f[i] += 2100.0;
      break;
    case 8:
      cec2021_cf02(&x[i * nx], &f[i], nx, OShift, M, 1);
      f[i] += 2200.0;
      break;
    case 9:
      cec2021_cf04(&x[i * nx], &f[i], nx, OShift, M, 1);
      f[i] += 2400.0;
      break;
    case 10:
      cec2021_cf05(&x[i * nx], &f[i], nx, OShift, M, 1);
      f[i] += 2500.0;
      break;
    default:
      perror("\nError: There are only 10 test functions in this test suite!\n");
      f[i] = 0.0;
      break;
    }
  }
}

void cec2021_rot_func(double *x, double *f, int nx, int mx, int func_num) {

  int cf_num = 10, i, j;
  if (func_num < 1 || func_num > 10) {
    perror("\nError: Test function is not defined.\n");
  }
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

    if (!(nx == 2 || nx == 10 || nx == 20)) {
      perror("\nError: Test functions are only defined for D=2,10,20.\n");
    }
    if (nx == 2 && (func_num == 5 || func_num == 6 || func_num == 7)) {
      perror("\nError:  NOT defined for D=2.\n");
    }

    /* Load Matrix M*/
    sprintf(FileName, "%s/M_%d_D%d.txt", extdata, func_num, nx);
    fpt = fopen(FileName, "r");
    if (fpt == NULL) {
      perror("\n Error: Cannot open M_d_Dd.txt for reading \n");
    }
    if (func_num < 7) {
      M = (double *)malloc(nx * nx * sizeof(double));
      if (M == NULL)
        perror("\nError: there is insufficient memory available!\n");
      for (i = 0; i < nx * nx; i++) {
        if (!fscanf(fpt, "%lf", &M[i])) {

          perror("\nError\n");
        }
      }
    } else {
      M = (double *)malloc(cf_num * nx * nx * sizeof(double));
      if (M == NULL)
        perror("\nError: there is insufficient memory available!\n");
      for (i = 0; i < cf_num * nx * nx; i++) {
        if (!fscanf(fpt, "%lf", &M[i])) {

          perror("\nError\n");
        }
      }
    }
    fclose(fpt);

    /* Load shift_data */
    sprintf(FileName, "%s/shift_data_%d_ns.txt", extdata, func_num);
    fpt = fopen(FileName, "r");
    if (fpt == NULL) {
      perror("\n Error: Cannot open shift_data_d_ns.txt for reading \n");
    }

    if (func_num < 7) {
      OShift = (double *)malloc(nx * sizeof(double));
      if (OShift == NULL)
        perror("\nError: there is insufficient memory available!\n");
      for (i = 0; i < nx; i++) {
        if (!fscanf(fpt, "%lf", &OShift[i])) {

          perror("\nError\n");
        }
      }
    } else {
      // OShift=(double *)malloc(nx*sizeof(double));
      OShift = (double *)malloc(nx * cf_num * sizeof(double));
      if (OShift == NULL)
        perror("\nError: there is insufficient memory available!\n");
      for (i = 0; i < cf_num - 1; i++) {
        for (j = 0; j < nx; j++) {
          if (!fscanf(fpt, "%lf", &OShift[i * nx + j])) {

            perror("\nError\n");
          }
        }
//        if (fscanf(fpt, "%*[^\n]%*c") != 1) {
//
//          perror("\nError\n");
//        }
      }
      for (j = 0; j < nx; j++) {
        if (!fscanf(fpt, "%lf", &OShift[nx * (cf_num - 1) + j])) {

          perror("\nError\n");
        }
      }
    }
    fclose(fpt);

    /* Load Shuffle_data */

    if (func_num >= 5 && func_num <= 7) {
      sprintf(FileName, "%s/shuffle_data_%d_D%d.txt", extdata, func_num, nx);
      fpt = fopen(FileName, "r");
      if (fpt == NULL) {
        perror("\n Error: Cannot open shuffle_data_d_Dd.txt for reading \n");
      }
      SS = (int *)malloc(nx * sizeof(int));
      if (SS == NULL)
        perror("\nError: there is insufficient memory available!\n");
      for (i = 0; i < nx; i++) {
        if (!fscanf(fpt, "%d", &SS[i])) {

          perror("\nError\n");
        }
      }
      fclose(fpt);
    }

    n_flag = nx;
    func_flag = func_num;
    ini_flag = 1;
    // perror("Function has been initialized!\n");
  }

  for (i = 0; i < mx; i++) {
    switch (func_num) {
    case 1:
      cec2021_bent_cigar_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      break;
    case 2:
      cec2021_schwefel_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      break;
    case 3:
      cec2021_bi_rastrigin_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      break;
    case 4:
      cec2021_grie_rosen_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      break;
    case 5:
      cec2021_hf01(&x[i * nx], &f[i], nx, OShift, M, SS, 1, 1);
      break;
    case 6:
      cec2021_hf06(&x[i * nx], &f[i], nx, OShift, M, SS, 1, 1);
      break;
    case 7:
      cec2021_hf05(&x[i * nx], &f[i], nx, OShift, M, SS, 1, 1);
      break;
    case 8:
      cec2021_cf02(&x[i * nx], &f[i], nx, OShift, M, 1);
      break;
    case 9:
      cec2021_cf04(&x[i * nx], &f[i], nx, OShift, M, 1);
      break;
    case 10:
      cec2021_cf05(&x[i * nx], &f[i], nx, OShift, M, 1);
      break;
    default:
      perror("\nError: There are only 10 test functions in this test suite!\n");
      break;
    }
  }
}
void cec2021_shift_rot_func(double *x, double *f, int nx, int mx,
                            int func_num) {

  int cf_num = 10, i, j;
  if (func_num < 1 || func_num > 10) {
    perror("\nError: Test function is not defined.\n");
  }
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

    if (!(nx == 2 || nx == 10 || nx == 20)) {
      perror("\nError: Test functions are only defined for D=2,10,20.\n");
    }
    if (nx == 2 && (func_num == 5 || func_num == 6 || func_num == 7)) {
      perror("\nError:  NOT defined for D=2.\n");
    }

    /* Load Matrix M*/
    sprintf(FileName, "%s/M_%d_D%d.txt", extdata, func_num, nx);
    fpt = fopen(FileName, "r");
    if (fpt == NULL) {
      perror("\n Error: Cannot open M_d_Dd.txt for reading \n");
    }
    if (func_num < 7) {
      M = (double *)malloc(nx * nx * sizeof(double));
      if (M == NULL)
        perror("\nError: there is insufficient memory available!\n");
      for (i = 0; i < nx * nx; i++) {
        if (!fscanf(fpt, "%lf", &M[i])) {

          perror("\nError\n");
        }
      }
    } else {
      M = (double *)malloc(cf_num * nx * nx * sizeof(double));
      if (M == NULL)
        perror("\nError: there is insufficient memory available!\n");
      for (i = 0; i < cf_num * nx * nx; i++) {
        if (!fscanf(fpt, "%lf", &M[i])) {

          perror("\nError\n");
        }
      }
    }
    fclose(fpt);

    /* Load shift_data */
    sprintf(FileName, "%s/shift_data_%d.txt", extdata, func_num);
    fpt = fopen(FileName, "r");
    if (fpt == NULL) {
      perror("\n Error: Cannot open shift_data_d.txt for reading \n");
    }

    if (func_num < 7) {
      OShift = (double *)malloc(nx * sizeof(double));
      if (OShift == NULL)
        perror("\nError: there is insufficient memory available!\n");
      for (i = 0; i < nx; i++) {
        if (!fscanf(fpt, "%lf", &OShift[i])) {

          perror("\nError\n");
        }
      }
    } else {
      // OShift=(double *)malloc(nx*sizeof(double));
      OShift = (double *)malloc(nx * cf_num * sizeof(double));
      if (OShift == NULL)
        perror("\nError: there is insufficient memory available!\n");
      for (i = 0; i < cf_num - 1; i++) {
        for (j = 0; j < nx; j++) {
          if (!fscanf(fpt, "%lf", &OShift[i * nx + j])) {

            perror("\nError\n");
          }
        }
//        if (fscanf(fpt, "%*[^\n]%*c") != 1) {
//
//          perror("\nError\n");
//        }
      }
      for (j = 0; j < nx; j++) {
        if (!fscanf(fpt, "%lf", &OShift[nx * (cf_num - 1) + j])) {

          perror("\nError\n");
        }
      }
    }
    fclose(fpt);

    /* Load Shuffle_data */

    if (func_num >= 5 && func_num <= 7) {
      sprintf(FileName, "%s/shuffle_data_%d_D%d.txt", extdata, func_num, nx);
      fpt = fopen(FileName, "r");
      if (fpt == NULL) {
        perror("\n Error: Cannot open shuffle_data_d_Dd.txt for reading \n");
      }
      SS = (int *)malloc(nx * sizeof(int));
      if (SS == NULL)
        perror("\nError: there is insufficient memory available!\n");
      for (i = 0; i < nx; i++) {
        if (!fscanf(fpt, "%d", &SS[i])) {

          perror("\nError\n");
        }
      }
      fclose(fpt);
    }

    n_flag = nx;
    func_flag = func_num;
    ini_flag = 1;
    // perror("Function has been initialized!\n");
  }

  for (i = 0; i < mx; i++) {
    switch (func_num) {
    case 1:
      cec2021_bent_cigar_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      break;
    case 2:
      cec2021_schwefel_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      break;
    case 3:
      cec2021_bi_rastrigin_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      break;
    case 4:
      cec2021_grie_rosen_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      break;
    case 5:
      cec2021_hf01(&x[i * nx], &f[i], nx, OShift, M, SS, 1, 1);
      break;
    case 6:
      cec2021_hf06(&x[i * nx], &f[i], nx, OShift, M, SS, 1, 1);
      break;
    case 7:
      cec2021_hf05(&x[i * nx], &f[i], nx, OShift, M, SS, 1, 1);
      break;
    case 8:
      cec2021_cf02_s(&x[i * nx], &f[i], nx, OShift, M, 1);
      break;
    case 9:
      cec2021_cf04_s(&x[i * nx], &f[i], nx, OShift, M, 1);
      break;
    case 10:
      cec2021_cf05_s(&x[i * nx], &f[i], nx, OShift, M, 1);
      break;
    default:
      perror("\nError: There are only 10 test functions in this test suite!\n");
      break;
    }
  }
}

void cec2021_shift_func(double *x, double *f, int nx, int mx, int func_num) {

  int cf_num = 10, i, j;
  if (func_num < 1 || func_num > 10) {
    perror("\nError: Test function is not defined.\n");
  }
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

    if (!(nx == 2 || nx == 10 || nx == 20)) {
      perror("\nError: Test functions are only defined for D=2,10,20.\n");
    }
    if (nx == 2 && (func_num == 5 || func_num == 6 || func_num == 7)) {
      perror("\nError:  NOT defined for D=2.\n");
    }

    /* Load Matrix M*/
    sprintf(FileName, "%s/M_%d_D%d_nr.txt", extdata, func_num, nx);
    fpt = fopen(FileName, "r");
    if (fpt == NULL) {
      perror("\n Error: Cannot open M_%d_D%d_nr.txt for reading \n");
    }
    if (func_num < 7) {
      M = (double *)malloc(nx * nx * sizeof(double));
      if (M == NULL)
        perror("\nError: there is insufficient memory available!\n");
      for (i = 0; i < nx * nx; i++) {
        if (!fscanf(fpt, "%lf", &M[i])) {

          perror("\nError\n");
        }
      }
    } else {
      M = (double *)malloc(cf_num * nx * nx * sizeof(double));
      if (M == NULL)
        perror("\nError: there is insufficient memory available!\n");
      for (i = 0; i < cf_num * nx * nx; i++) {
        if (!fscanf(fpt, "%lf", &M[i])) {

          perror("\nError\n");
        }
      }
    }
    fclose(fpt);

    /* Load shift_data */
    sprintf(FileName, "%s/shift_data_%d.txt", extdata, func_num);
    fpt = fopen(FileName, "r");
    if (fpt == NULL) {
      perror("\n Error: Cannot open shift_data_d.txt for reading \n");
    }

    if (func_num < 7) {
      OShift = (double *)malloc(nx * sizeof(double));
      if (OShift == NULL)
        perror("\nError: there is insufficient memory available!\n");
      for (i = 0; i < nx; i++) {
        if (!fscanf(fpt, "%lf", &OShift[i])) {

          perror("\nError\n");
        }
      }
    } else {
      // OShift=(double *)malloc(nx*sizeof(double));
      OShift = (double *)malloc(nx * cf_num * sizeof(double));
      if (OShift == NULL)
        perror("\nError: there is insufficient memory available!\n");
      for (i = 0; i < cf_num - 1; i++) {
        for (j = 0; j < nx; j++) {
          if (!fscanf(fpt, "%lf", &OShift[i * nx + j])) {

            perror("\nError\n");
          }
        }
//        if (fscanf(fpt, "%*[^\n]%*c") != 1) {
//
//         perror("\nError\n");
//        }
      }
      for (j = 0; j < nx; j++) {
        if (!fscanf(fpt, "%lf", &OShift[nx * (cf_num - 1) + j])) {

          perror("\nError\n");
        }
      }
    }
    fclose(fpt);

    /* Load Shuffle_data */

    if (func_num >= 5 && func_num <= 7) {
      sprintf(FileName, "%s/shuffle_data_%d_D%d.txt", extdata, func_num, nx);
      fpt = fopen(FileName, "r");
      if (fpt == NULL) {
        perror("\n Error: Cannot open shuffle_data_d_Dd.txt for reading \n");
      }
      SS = (int *)malloc(nx * sizeof(int));
      if (SS == NULL)
        perror("\nError: there is insufficient memory available!\n");
      for (i = 0; i < nx; i++) {
        if (!fscanf(fpt, "%d", &SS[i])) {

          perror("\nError\n");
        }
      }
      fclose(fpt);
    }

    n_flag = nx;
    func_flag = func_num;
    ini_flag = 1;
    // perror("Function has been initialized!\n");
  }

  for (i = 0; i < mx; i++) {
    switch (func_num) {
    case 1:
      cec2021_bent_cigar_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      break;
    case 2:
      cec2021_schwefel_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      break;
    case 3:
      cec2021_bi_rastrigin_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      break;
    case 4:
      cec2021_grie_rosen_func(&x[i * nx], &f[i], nx, OShift, M, 1, 1);
      break;
    case 5:
      cec2021_hf01(&x[i * nx], &f[i], nx, OShift, M, SS, 1, 1);
      break;
    case 6:
      cec2021_hf06(&x[i * nx], &f[i], nx, OShift, M, SS, 1, 1);
      break;
    case 7:
      cec2021_hf05(&x[i * nx], &f[i], nx, OShift, M, SS, 1, 1);
      break;
    case 8:
      cec2021_cf02_s(&x[i * nx], &f[i], nx, OShift, M, 1);
      break;
    case 9:
      cec2021_cf04_s(&x[i * nx], &f[i], nx, OShift, M, 1);
      break;
    case 10:
      cec2021_cf05_s(&x[i * nx], &f[i], nx, OShift, M, 1);
      break;
    default:
      perror("\nError: There are only 10 test functions in this test suite!\n");
      break;
    }
  }
}
