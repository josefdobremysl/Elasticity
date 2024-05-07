#include "mesh_processor.h"
#include "Dirichlet.h"
#include "Neumann.h"
#include "function.h"
#include "Assemblematrix.h"

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <armadillo>
#include <omp.h>

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <tuple>
#include <cmath>
#include <numeric>

using namespace Eigen;
using namespace std;
using namespace arma;


SparseMatrixVectorResult DirichletBoundaryCondition(const string& meshFile, const int& tag) {

	int num_th = 4;
	omp_set_num_threads(num_th); // Nastavení počtu vláken

	SparseMatrixVectorResult resultS;
	MatrixVectorResult result1 = GetTriangles(meshFile);
	vector<int> vectorOfTags = result1.vector1;
	int nTri = result1.integerResult;		// poèet trojúhelníkùD
	int nNodes = result1.matrix3.size();			// poèet bodù
	int nNodes2 = 2 * nNodes;
	vector<int> I(nNodes2,0);
	vector<int> J(nNodes2, 0);
	vector<double> VAL(nNodes2, 0);

	vector<double> bD = NeumannBoundaryCondition(meshFile, tag);


	vector<double> x(nNodes2, 0);
	double xx = 0;

	double cD = 10e9;
	// Vynásobení vektoru koeficientem
	for (double& element : bD) {
		element *= cD;
	}

	//vector<vector<vector<double>>> basisFE_ref = {
	//	{{(1 - 1 / sqrt(3)) / 2, 0}, {0, 20 + (1 - 1 / sqrt(3)) / 2}, {30 + (1 + 1 / sqrt(3)) / 2, 0}, {0, 40 + (1 + 1 / sqrt(3)) / 2}},
	//	{{(1 + 1 / sqrt(3)) / 2, 0}, {0, 60 + (1 + 1 / sqrt(3)) / 2}, {70 + (1 - 1 / sqrt(3)) / 2, 0}, {0, 80 + (1 - 1 / sqrt(3)) / 2}}
	//};
	vector<vector<vector<double>>> basisFE_ref = {
	{{(1 - 1 / sqrt(3)) / 2, 0}, { (1 + 1 / sqrt(3)) / 2,0}, {0,(1 - 1 / sqrt(3)) / 2}, {0, (1 + 1 / sqrt(3)) / 2}},
	{{(1 + 1 / sqrt(3)) / 2,0}, {  (1 - 1 / sqrt(3)) / 2,0}, {0,(1 + 1 / sqrt(3)) / 2}, {0, (1 - 1 / sqrt(3)) / 2}}
	};
	double nQuadrature = 2;
	double w = 0.5;
	int index_IJ = 0;

	#pragma omp for 
	for (int k = 0; k < vectorOfTags.size(); k++) {

		if (vectorOfTags[k] == tag) {

			int E = result1.matrix2[k][0] - 1;					// Oznaèení vrcholù hran
			int D = result1.matrix2[k][1] - 1;
			vector<int> edgeVertex = { D, E };

			vector<double> cooD = result1.matrix3[D];		// Souøadnice vrcholù hran
			vector<double> cooE = result1.matrix3[E];

			//vector<vector<double>> matA = {
			//	{ cooE[0] - cooD[0] },
			//	{ cooE[1] - cooD[1] },
			//};
			double edgeLenght = sqrt((cooD[0] - cooE[0]) * (cooD[0] - cooE[0]) + (cooD[1] - cooE[1]) * (cooD[1] - cooE[1]));
			vector<int> index = { D, E, D + nNodes, E + nNodes };

			for (int i = 0; i < 4; i++) {
				for (int j = 0; j < 4; j++) {
					double k_val = 0;
					for (int l = 0; l < nQuadrature; l++) {
						vector<double> v = basisFE_ref[l][j];
						//cout <<" i "<< i << " j " << j << " v " << v[0] <<" "<< v[1] << endl;
						//cout << "v  " << " i " << i << " j " << j << " l " << l << v[0] << "  " << v[1] << endl;
						vector<double> u = basisFE_ref[l][i];
						//cout << "	 u " << u[0] <<" "<< u[1] << endl;
						k_val = k_val + edgeLenght * w * dotprod(u, v);
						//cout << "dotprod(u, v)		" << dotprod(u, v) << endl;

					};
					k_val = cD * k_val;
					I[index_IJ] = index[i];
					J[index_IJ] = index[j];
					VAL[index_IJ] = k_val;
					index_IJ++;
				};
			};
		};

	};
	SparseMatrixVectorResult assembledMatrix = AssembledMatrix(I, J, VAL, bD);



	return assembledMatrix;
}