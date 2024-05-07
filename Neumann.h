#ifndef NEUMANN_H
#define NEUMANN_H

#include <string>
#include <vector>
#include <tuple>
#include <Eigen/Sparse>

using namespace Eigen;
using namespace std;


vector<double> NeumannBoundaryCondition(const string& meshFile, const int& tag);

#endif // NEUMANN_H