#pragma once


#include <iostream>
#include <fstream>
#include <cmath>
//#include <Eigen/Dense>

using namespace std;

double Kernel(double r);
void Neighbours(double* nb_x, double* nb_y, double x, double y);
/**class ComputeEpsilon
{
private:
  int Nodes;
  //Eigen
  Eigen::MatrixXd A;
  Eigen::VectorXd epsilon_sol;
public:
  ComputeEpsilon(int N);
  ~ComputeEpsilon();
  double GetEpsilonK(int k){return epsilon_sol(k);};
  double * GetEpsilon(){return epsilon_sol.data();};
  void BuildMatrix(double * dotsx, double * dotsy, double ds);
  void SolveA();
  void ShowEpsilon(){cout << epsilon_sol;};
  void ShowA(){cout << A;};
};**/
