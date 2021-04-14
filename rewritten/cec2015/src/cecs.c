#include "../include/cec2015.h"
#include <stddef.h>

void print_vec(size_t n, double vec[n]) {
  printf("{");
  for (size_t i = 0; i < n; ++i) {
    printf("%f,", vec[i]); 
  }
  printf("}\n");
}


int main() {
  char *dataPath = "/home/ewarchul/cecs/inst/extdata/cec2015";
  char **ptrDataPath = &dataPath;
  int func = 12;
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
  for (int i = 0; i < N; ++i) {
    X[i] = 0.23;
  }
  cec2015(ptrDataPath, ptrFunc, X, ptrRow, ptrN, f);   
  print_vec(N, f);
  free(f);
  return 0;
}
