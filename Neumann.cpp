#include "mesh_processor.h"
#include "Dirichlet.h"
#include "Neumann.h"
#include "function.h"

#include <Eigen/Sparse>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <tuple>
#include <cmath>
#include <armadillo>
#include <numeric>
#include <omp.h>

using namespace Eigen;
using namespace std;
using namespace arma;

vector<double> NeumannBoundaryCondition(const string& meshFile, const int& tag) {

	int num_th = 2;
	omp_set_num_threads(num_th); // Nastavení počtu vláken

	MatrixVectorResult loadedMesh = GetTriangles(meshFile);
	vector<int> vectorOfTags = loadedMesh.vector1;
	vector<vector<int>> lineVertex = loadedMesh.matrix2;
	int size = loadedMesh.integerResult2;
	vector<vector<int>> markedEdges(size, vector<int>(2, 0));

	vector<vector<double>> nM(size, vector<double>(2, 0));
	int nNodes = 2 * loadedMesh.matrix3.size();
	int N = loadedMesh.matrix3.size();
	vector<double> bNinPoints(nNodes, 0);

	double w = 1;

	#pragma omp for 
	for (int i = 0; i < size; i++) {
		if (vectorOfTags[i] == tag) {
			int D = lineVertex[i][0] - 1;
			int E = lineVertex[i][1] - 1;
			//cout << "D:" << D << "E:" << E << endl;
			vector<double> cooD = loadedMesh.matrix3[D];		// Souøadnice vrcholù trojúhelníku k
			vector<double> cooE = loadedMesh.matrix3[E];

			vector<double> n = { cooE[1] - cooD[1], cooD[0] - cooE[0] };

			//// nahradit toto matici pro x a y s body kvadratury v globalni siti
			//double cooF1X = cooD[0] + (cooE[0] - cooD[0]) * (0.5 - 1 / (2 * sqrt(3)));
			//double cooF1Y = cooD[1] + (cooE[1] - cooD[1]) * (0.5 - 1 / (2 * sqrt(3)));
			//double cooF2X = cooD[0] + (cooE[0] - cooD[0]) * (0.5 + 1 / (2 * sqrt(3)));
			//double cooF2Y = cooD[1] + (cooE[1] - cooD[1]) * (0.5 + 1 / (2 * sqrt(3)));
			//double edgeLength = sqrt(pow(cooD[0] - cooE[0], 2) + pow(cooD[1] - cooE[1], 2));

			//// hodnoty bazovych funcki na ref. hrane

			//vector<double> fN1 = { FunctionN1(cooF1X, cooF1Y, tag), FunctionN1(cooF2X, cooF2Y, tag) };
			//vector<double> fN2 = { FunctionN2(cooF1X, cooF1Y, tag), FunctionN2(cooF2X, cooF2Y, tag) };

			//// bázové funkce fi

			//vector<double> fi_E1_F1 = { 0.5 - 1 / (2 * sqrt(3)) , 0 };
			//vector<double> fi_E2_F1 = { 0 , 0.5 - 1 / (2 * sqrt(3)) };
			//vector<double> fi_D1_F1 = { 0.5 + 1 / (2 * sqrt(3)) , 0 };
			//vector<double> fi_D2_F1 = { 0 , 0.5 + 1 / (2 * sqrt(3)) };

			//vector<double> fi_E1_F2 = { 0.5 + 1 / (2 * sqrt(3)) , 0 };
			//vector<double> fi_E2_F2 = { 0 , 0.5 + 1 / (2 * sqrt(3)) };
			//vector<double> fi_D1_F2 = { 0.5 - 1 / (2 * sqrt(3)) , 0 };
			//vector<double> fi_D2_F2 = { 0 , 0.5 - 1 / (2 * sqrt(3)) };

			//double bN_E = inner_product(fN1.begin(), fN1.end(), fi_E1_F1.begin(), 0.0) + inner_product(fN1.begin(), fN1.end(), fi_E1_F2.begin(), 0.0);
			//double bN_EpN = inner_product(fN2.begin(), fN2.end(), fi_E2_F1.begin(), 0.0) + inner_product(fN2.begin(), fN2.end(), fi_E2_F2.begin(), 0.0);
			//double bN_D = inner_product(fN1.begin(), fN1.end(), fi_D1_F1.begin(), 0.0) + inner_product(fN1.begin(), fN1.end(), fi_D1_F2.begin(), 0.0);
			//double bN_DpN = inner_product(fN2.begin(), fN2.end(), fi_D2_F1.begin(), 0.0) + inner_product(fN2.begin(), fN2.end(), fi_D2_F2.begin(), 0.0);
			//
			//bNinPoints[E] += edgeLength * w * (bN_E + bN_D);
			//bNinPoints[E + N] += edgeLength * w * (bN_EpN + bN_DpN);
			//bNinPoints[D] += edgeLength * w * (bN_E + bN_D);
			//bNinPoints[D + N] += edgeLength * w * (bN_EpN + bN_DpN);

			//////////////////////////////////////////////////////////////////////////////

			double Q1 = 0.5 - 1 / (2 * sqrt(3));
			double Q2 = 0.5 + 1 / (2 * sqrt(3));

			vector<double> cooF1 = { cooD[0] + (cooE[0] - cooD[0]) * Q1, cooD[1] + (cooE[1] - cooD[1]) * Q1 };
			vector<double> cooF2 = { cooD[0] + (cooE[0] - cooD[0]) * Q2, cooD[1] + (cooE[1] - cooD[1]) * Q2 };
			double edgeLength = sqrt(pow(cooD[0] - cooE[0], 2) + pow(cooD[1] - cooE[1], 2));

			vector<double> fN1 = { FunctionN1(cooF1[0], cooF1[1], tag, n), FunctionN1(cooF2[0], cooF2[1], tag, n) };
			vector<double> fN2 = { FunctionN2(cooF1[0], cooF1[1], tag, n), FunctionN2(cooF2[0], cooF2[1], tag, n) };

			vector<double> fi_D_F1 = { fi1(Q1),fi1(Q2) };
			vector<double> fi_D_F2 = { fi1(Q1),fi1(Q2) };
			vector<double> fi_E_F1 = { fi2(Q1),fi2(Q2) };
			vector<double> fi_E_F2 = { fi2(Q1),fi2(Q2) };

			bNinPoints[D] += edgeLength * w * dotprod(fN1, fi_D_F1);
			bNinPoints[E] += edgeLength * w * dotprod(fN1, fi_E_F1);
			bNinPoints[D + N] += edgeLength * w * dotprod(fN2, fi_D_F2);
			bNinPoints[E + N] += edgeLength * w * dotprod(fN2, fi_E_F2);

			//////////////////////////////////////////////////////////////////////////////

			//double Q = 0.5;

			//vector<double> cooF = { cooD[0] + (cooE[0] - cooD[0]) * Q, cooD[1] + (cooE[1] - cooD[1]) * Q };

			//double edgeLength = sqrt(pow(cooD[0] - cooE[0], 2) + pow(cooD[1] - cooE[1], 2));

			//double fN1 = FunctionN1(cooF[0], cooF[1], tag, n);
			//double fN2 = FunctionN2(cooF[0], cooF[1], tag, n);

			//bNinPoints[D] += edgeLength * w * fN1*0.5;
			//bNinPoints[E] += edgeLength * w * fN1 * 0.5;
			//bNinPoints[D + N] += edgeLength * w * fN2 * 0.5;
			//bNinPoints[E + N] += edgeLength * w * fN2 * 0.5;

		}
	}

	//for (const double& element : bNinPoints) {
	//	std::cout << element << " ";
	//}

	return bNinPoints;

};