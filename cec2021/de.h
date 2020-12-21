/*
  L-SHADE implemented by C++ for Special Session & Competition on Real-Parameter Single Objective Optimization at CEC-2014

  Version: 1.0   Date: 16/Apr/2014
  Written by Ryoji Tanabe (rt.ryoji.tanabe [at] gmail.com)
*/

#ifndef _HEADER_H_
#define _HEADER_H_

#include <stdlib.h>
#include<iostream>
#include <iomanip>
#include <string.h>
#include <vector>
#include<math.h>

using namespace std;

typedef  double variable;
typedef variable *Individual;
typedef  double Fitness;

extern int g_function_number;
extern int g_problem_size;
extern unsigned int g_max_num_evaluations;
extern int function_name;

extern int g_pop_size;
extern int g_memory_size;
extern double g_p_best_rate;
extern double g_arc_rate;

void cec21_bias_shift_rot_func(double *, double *,int,int,int);
void cec21_bias_shift_func(double *, double *,int,int,int);
void cec21_bias_rot_func(double *, double *,int,int,int);
void cec21_shift_rot_func(double *, double *,int,int,int);
void cec21_rot_func(double *, double *,int,int,int);
void cec21_shift_func(double *, double *,int,int,int);
void cec21_bias_func(double *, double *,int,int,int);
void cec21_basic_func(double *, double *,int,int,int);

class searchAlgorithm {
public:
  virtual Fitness run() = 0;
protected:
  void evaluatePopulation(const vector<Individual> &pop, vector<Fitness> &fitness);
  void initializeFitnessFunctionParameters();

  void initializeParameters();
  Individual makeNewIndividual();
  void modifySolutionWithParentMedium(Individual child, Individual parent);
  void setBestSolution(const vector<Individual> &pop, const vector<Fitness> &fitness, Individual &bsf_solution, Fitness &bsf_fitness);

  //Return random value with uniform distribution [0, 1)
  inline double randDouble() {
    return (double)rand() / (double) RAND_MAX;
  }

  /*
    Return random value from Cauchy distribution with mean "mu" and variance "gamma"
    http://www.sat.t.u-tokyo.ac.jp/~omi/random_variables_generation.html#Cauchy
  */
  inline double cauchy_g(double mu, double gamma) {
    return mu + gamma * tan(M_PI*(randDouble() - 0.5));
  }

  /*
    Return random value from normal distribution with mean "mu" and variance "gamma"
    http://www.sat.t.u-tokyo.ac.jp/~omi/random_variables_generation.html#Gauss
  */
  inline double gauss(double mu, double sigma){
    return mu + sigma * sqrt(-2.0 * log(randDouble())) * sin(2.0 * M_PI * randDouble());
  }

  //Recursive quick sort with index array
  template<class VarType>
    void sortIndexWithQuickSort(VarType array[], int first, int last, int index[]) {
    VarType x = array[(first + last) / 2];
    int i = first;
    int j = last;
    VarType temp_var = 0;
    int temp_num = 0;

    while (true) {
      while (array[i] < x) i++;
      while (x < array[j]) j--;
      if (i >= j) break;

      temp_var = array[i];
      array[i] = array[j];
      array[j] = temp_var;

      temp_num = index[i];
      index[i] = index[j];
      index[j] = temp_num;

      i++;
      j--;
    }

    if (first < (i -1)) sortIndexWithQuickSort(array, first, i - 1, index);
    if ((j + 1) < last) sortIndexWithQuickSort(array, j + 1, last, index);
  }

  int function_number;
  int problem_size;
  variable max_region;
  variable min_region;
  Fitness optimum;
  // acceptable error value
  Fitness epsilon;
  unsigned int max_num_evaluations;
  int pop_size;
};

class LSHADE: public searchAlgorithm {
public:
  virtual Fitness run();
  void setSHADEParameters();
  void reducePopulationWithSort(vector<Individual> &pop, vector<Fitness> &fitness);
  void operateCurrentToPBest1BinWithArchive(const vector<Individual> &pop, Individual child, int &target, int &p_best_individual, variable &scaling_factor, variable &cross_rate, const vector<Individual> &archive, int &arc_ind_count);

  int arc_size;
  double arc_rate;
  variable p_best_rate;
  int memory_size;
  int reduction_ind_num;
};

#endif


