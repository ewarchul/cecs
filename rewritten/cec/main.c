#include "cec.h"

int main() {
  char *dataPath = "../data/cec2017";
  int dims[4] = {10, 30, 50, 100};
  for (int d = 0; d < 4; ++d) {
    for (int fn = 1; fn < 31; ++fn) {
      char **ptrDataPath = &dataPath;
      int func = fn;
      int *ptrFunc = &func;
      int row = 1;
      int cec_ = 2017;
      int *ptrCec = &cec_;
      int *ptrRow = &row;
      int N = dims[d];
      int *ptrN = &N;
      double *f = malloc(N * sizeof(double));
      double *X = malloc(dims[d] * sizeof(double));

      for (int i = 0; i < N; ++i) {
        f[i] = 0;
      }

      for (int i = 0; i < N; ++i) {
        X[i] = 5;
      }
      cec(ptrDataPath, ptrCec, ptrFunc, X, ptrRow, ptrN, f, "basic");
      printf("(Dimension, Function, Value): (%d, %d, %.2f)\n", dims[d], fn,
             f[0]);
      free(f);
      free(X);
    }
  }
  return 0;
}
