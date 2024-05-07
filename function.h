#ifndef FUNCTION_H
#define FUNCTION_H

#include "mesh_processor.h"

#include <string>
#include <vector>
#include <tuple>
#include <Eigen/Sparse>

using namespace Eigen;
using namespace std;

double Function(double& X, double& Y);
double FunctionN1(double& X, double& Y, const int& tag, vector<double>& n);
double FunctionN2(double& X, double& Y, const int& tag, vector<double>& n);
double fi1(double& x);
double fi2(double& x);

double dotprod(vector<double>& u, vector<double>& v);

vector<double> FunctionV(double& X, double& Y);



#endif // FUNCTION_H