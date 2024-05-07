#ifndef STIFFNESSMATRIXQ_H
#define STIFFNESSMATRIXQ_H

#include "mesh_processor.h"

#include <string>
#include <vector>
#include <tuple>
#include <Eigen/Sparse>

using namespace Eigen;
using namespace std;

SparseMatrixVectorResult StiffnessMatrix(const string& meshFile);

#endif // STIFFNESSMATRIXQ_H