#include "../include/cec2017.h"
#include <stddef.h>

void print_vec(size_t n, double vec[n]) {
  printf("{");
  for (size_t i = 0; i < n; ++i) {
    printf("%f,", vec[i]); 
  }
  printf("}\n");
}


int main() {
  char *dataPath = "/home/ewarchul/cecs/inst/extdata/cec2017";
  char **ptrDataPath = &dataPath;
  int func = 25;
  int *ptrFunc = &func;
  int row = 1;
  int *ptrRow = &row;
  int N = 100;
  int *ptrN = &N;
  double *f = malloc(N * sizeof(double));
  for (int i = 0; i < N; ++i) {
    f[i] = 0;
  }
  double X[100]; 
  for (int i = 0; i < 100; ++i) {
    X[i] = 0.1;
  }
  cec2017(ptrDataPath, ptrFunc, X, ptrRow, ptrN, f);   
  free(f);
  return 0;
}
