#include "../include/cec2021.h"
#include <stddef.h>

void print_vec(size_t n, double vec[n]) {
  printf("{");
  for (size_t i = 0; i < n; ++i) {
    printf("%f,", vec[i]);
  }
  printf("}\n");
}

int main() {
  char *dataPath = "/home/ewarchul/cecs/inst/extdata/cec2021";
  for (int fn = 1; fn < 11; ++fn) {
   // printf("Function: %d\n", fn);
    char **ptrDataPath = &dataPath;
    int func = fn;
    int *ptrFunc = &func;
    int row = 1;
    int *ptrRow = &row;
    int N = 10;
    int *ptrN = &N;
    double *f = malloc(N * sizeof(double));
    for (int i = 0; i < N; ++i) {
      f[i] = 0;
    }
    double X[10];
    for (int i = 0; i < N; ++i) {
      X[i] = 0.1;
    }
    cec2021(ptrDataPath, ptrFunc, X, ptrRow, ptrN, f, "bias_shift_rot");
  //  print_vec(N, f);
    free(f);
  }
  return 0;
}
