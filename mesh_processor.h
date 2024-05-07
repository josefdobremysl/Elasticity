#ifndef MESH_PROCESSOR_H
#define MESH_PROCESSOR_H

#include <string>
#include <vector>
#include <tuple>
#include <Eigen/Sparse>

using namespace Eigen;
using namespace std;

struct MatrixVectorResult {
    vector<vector<int>> matrix1;
    vector<vector<int>> matrix2;
    vector<vector<double>> matrix3;
    vector<int> vector1;
    vector<int> vector2;
    int integerResult;
    int integerResult2;
};
struct SparseMatrixVectorResult {
    SparseMatrix<double> sparsematrix;
    vector<double> right_side_vector;
};
struct twovectors {
    vector<double> bN;
    vector<double> bz;
};

MatrixVectorResult GetTriangles(const string& meshFile);

#endif // MESH_PROCESSOR_H