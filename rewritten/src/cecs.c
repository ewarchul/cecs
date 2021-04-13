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
  int func = 1;
  int *ptrFunc = &func;
  int row = 1;
  int *ptrRow = &row;
  int N = 10;
  int *ptrN = &N;
  double *f = malloc(N * sizeof(double));
  for (int i = 0; i < N; ++i) {
    f[i] = 0;
  }
  double X[10] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

  cec2017(ptrDataPath, ptrFunc, X, ptrRow, ptrN, f);   

  free(f);
  return 0;
}
