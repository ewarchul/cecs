#define UNIT_TESTING 1

#include "cec.h"

#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#include <cmocka.h>

typedef struct CecInfo CecInfo;
struct CecInfo { 
  int m_version;
  int m_fnNum;
  int m_size;
  int m_dims[];
};

static CecInfo* mk_CecInfo(size_t t_n, int t_version, int t_fnNum) {
  CecInfo *ci = malloc(sizeof(CecInfo) + t_n * sizeof(int));
  if (!ci) {
    perror("Error: mk_CecInfo constructor");
    exit(-1);
  }
  for (size_t i = 0; i < t_n; ++i) {
    ci->m_dims[i] = 0;
  }
  ci->m_size = t_n;
  ci->m_version = t_version;
  ci->m_fnNum = t_fnNum;
  return ci; 
}

static void destroy_CecInfo(CecInfo* ci) {
  free(ci);
}

static void null_test_success(void **state) { (void)state; }

static void first_cmocka_test(void **state) {
  (void)state;
  assert_int_equal(5, 5);
}


/*static void test_if_cecs2014_are_equal(void **state) {*/
  /*(void)state;*/
/*}*/

int main() {
  CecInfo* bInfo = mk_CecInfo(4, 2014, 30);
  printf("siema %d\n", bInfo->m_version); 
  bInfo->m_dims[0] = 1;
  bInfo->m_dims[1] = 1;
  bInfo->m_dims[2] = 1;
  bInfo->m_dims[3] = 1;
  destroy_CecInfo(bInfo);
  /*char *dataPath = "../data/cec2017";*/
  /*char **ptrDataPath = &dataPath;*/
  /*int func = 7;*/
  /*int *ptrFunc = &func;*/
  /*int row = 1;*/
  /*int cec_ = 2017;*/
  /*int *ptrCec = &cec_;*/
  /*int *ptrRow = &row;*/
  /*int N = 10;*/
  /*int *ptrN = &N;*/
  /*double *f = malloc(N * sizeof(double));*/
  /*double *X = malloc(N * sizeof(double));*/

  /*for (int i = 0; i < N; ++i) {*/
    /*f[i] = 0;*/
  /*}*/

  /*for (int i = 0; i < N; ++i) {*/
    /*X[i] = 5;*/
  /*}*/
  /*cec(ptrDataPath, ptrCec, ptrFunc, X, ptrRow, ptrN, f, "basic");*/

  const struct CMUnitTest tests[] = {
      cmocka_unit_test(null_test_success),
      cmocka_unit_test(first_cmocka_test),
  };
//  free(X);
//  free(f);

  return cmocka_run_group_tests(tests, NULL, NULL);
}
