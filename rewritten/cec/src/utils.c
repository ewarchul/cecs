#include "utils.h"

double getFunctionBias(const int biasFlag, const int fnNumber) {
  double bias = 0.0;
  double fnBiasDict[10] = {100.0,  1100.0, 700.0,  1900.0, 1700.0,
                           1600.0, 2100.0, 2200.0, 2400.0, 2500.0};
  if (biasFlag) {
    bias = fnBiasDict[fnNumber - 1];
  } else {
    bias = 0.0;
  }
  return bias;
}
