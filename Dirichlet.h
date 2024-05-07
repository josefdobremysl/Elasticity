#ifndef DIRICHLET_H
#define DIRICHLET_H

#include "mesh_processor.h"

#include <string>
#include <vector>
#include <tuple>
#include <Eigen/Sparse>

using namespace Eigen;
using namespace std;

SparseMatrixVectorResult DirichletBoundaryCondition(const string& meshFile, const int& tag);

#endif // DIRICHLET_H