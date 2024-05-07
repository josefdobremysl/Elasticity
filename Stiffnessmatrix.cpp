#include "mesh_processor.h"
#include "Dirichlet.h"
#include "Neumann.h"
#include "function.h"
#include "Assemblematrix.h"

#include <Eigen/Sparse>
#include <armadillo>
#include <omp.h>

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <tuple>
#include <numeric>
#include <chrono>	
using namespace Eigen;
using namespace std;
using namespace arma;
using Time = std::chrono::steady_clock;
using ms = std::chrono::milliseconds;


SparseMatrixVectorResult StiffnessMatrix(const string& meshFile) {

	int num_th = 1;
	omp_set_num_threads(num_th); // Nastavení počtu vláken

	double E = 100000;					// Materiálové vlastnosti
	double sigma = 0.3;
	double lambda = E * sigma / ((1 + sigma) * (1 - 2 * sigma));
	double mu = E / (2 * (1 + sigma));

	const auto start = Time::now();

	MatrixVectorResult result = GetTriangles(meshFile);

	const auto end = Time::now();
	const auto diff = std::chrono::duration_cast<ms>(end - start).count();
	std::cout << " GetTriangles time = " << diff << " ms" << "\n";

	int nTri = result.integerResult;		// poèet trojúhelníkù
	int nNodes = result.matrix3.size();			// poèet bodù
	int nNodes2 = 2 * nNodes;
	int nNodesIJ = nNodes2 * 36;
	vector<int> I(nNodesIJ, 0);
	vector<int> J(nNodesIJ, 0);
	vector<double> VAL(nNodesIJ, 0);
	vector<double> right_side_vector(nNodes2, 0);		// vektor pravé strany
	cout << nNodes2 << endl;


	// Vytvoøení matice basisFE_ref
	mat basisFE_ref = { // bazove f-ce v x-smeru a pak baz. fce v y-smeru
		{0.5, 0.0, 0.5,   0,   0,   0},
		{0.5, 0.5, 0.0,   0,   0,   0},
		{0.0, 0.5, 0.5,   0,   0,   0},
		{  0,   0,   0, 0.5, 0.0, 0.5},
		{  0,   0,   0, 0.5, 0.5, 0.0},
		{  0,   0,   0, 0.0, 0.5, 0.5}
	};
	// PD podle X baz. fci v kvadr. uzlech
	mat gradFEx_ref = {
		{-1, -1, -1,  0,  0,  0},
		{ 1,  1,  1,  0,  0,  0},
		{ 0,  0,  0,  0,  0,  0},
		{ 0,  0,  0, -1, -1, -1},
		{ 0,  0,  0,  1,  1,  1},
		{ 0,  0,  0,  0,  0,  0}
	};
	// PD podle Y baz. fci v kvadr. uzlech
	mat gradFEy_ref = {
		{-1, -1, -1,  0,  0,  0},
		{ 0,  0,  0,  0,  0,  0},
		{ 1,  1,  1,  0,  0,  0},
		{ 0,  0,  0, -1, -1, -1},
		{ 0,  0,  0,  0,  0,  0},
		{ 0,  0,  0,  1,  1,  1}
	};

	double nQuadrature = 3;
	double w = 1 / nQuadrature;
	int index_IJ = 0;

	#pragma omp for 
	for (int k = 0; k < nTri; k++) {
		int A = result.matrix1[k][0] - 1;					// Oznaèení vrcholù trojúhelníku k
		int B = result.matrix1[k][1] - 1;
		int C = result.matrix1[k][2] - 1;
		vector<int> verTri = { A, B, C };

		vector<double> cooA = result.matrix3[A];		// Souøadnice vrcholù trojúhelníku k
		vector<double> cooB = result.matrix3[B];
		vector<double> cooC = result.matrix3[C];

		mat CoQ = {
			{ (cooA[0] + cooB[0]) / 2, (cooB[0] + cooC[0]) / 2, (cooC[0] + cooA[0]) / 2},
			{ (cooA[1] + cooB[1]) / 2, (cooB[1] + cooC[1]) / 2, (cooC[1] + cooA[1]) / 2},
		};

		mat matA = {
			{ cooB[0] - cooA[0], cooC[0] - cooA[0] },
			{ cooB[1] - cooA[1], cooC[1] - cooA[1] },
		};

		double detA = matA(0, 0) * matA(1, 1) - matA(0, 1) * matA(1, 0);

		mat invA = {
			{  matA(1,1) / detA, -matA(0,1) / detA },
			{ -matA(1,0) / detA, matA(0,0) / detA  },
		};

		mat gradFEx = invA(0, 0) * gradFEx_ref + invA(1, 0) * gradFEy_ref;
		mat gradFEy = invA(0, 1) * gradFEx_ref + invA(1, 1) * gradFEy_ref;

		vector<int> index = { A, B, C, nNodes + A, nNodes + B, nNodes + C };

		for (int i = 0; i < 6; ++i) {
			for (int j = 0; j < 6; ++j) {
				double k_val = 0;
				for (int l = 0; l < nQuadrature; ++l) {
					vec ee_i = { gradFEx(i,l), gradFEy(i, nQuadrature + l), gradFEx(i, nQuadrature + l) + gradFEy(i,l) };
					double div_i = gradFEx(i, l) + gradFEy(i, nQuadrature + l);
					vec sigma_i = vec{ 1.0, 1.0, 0.0 } *lambda * div_i + 2 * mu * ee_i;
					vector<double> ee_j = { gradFEx(j,l), gradFEy(j, nQuadrature + l), gradFEx(j, nQuadrature + l) + gradFEy(j,l) };

					k_val = k_val + detA / 2 * w * inner_product(sigma_i.begin(), sigma_i.end(), ee_j.begin(), 0.0);
					//cout << "kval" << k_val << endl;
					//cout << sigma_i << endl;


				};
				// zapis vysledneho soucinu vektoru-tenzoru na spravne misto do tripletu 
				I[index_IJ]=index[i];
				J[index_IJ] = index[j];
				VAL[index_IJ] = k_val;
				index_IJ++;
			}
			vector<double> fV(6, 0);
			for (int m = 0; m < nQuadrature; ++m) {
				double X = CoQ(0, m);
				double Y = CoQ(1, m);
				//double fx = 3 * X;	// Function(X, Y);
				//double fy = 2 * Y;
				fV[m] = FunctionV(X, Y)[0];
				fV[m + 3] = FunctionV(X, Y)[1];
			}

			arma::vec basisFE_ref_i = basisFE_ref.col(i);
			right_side_vector[index[i]] = 0.5 * detA * w * inner_product(fV.begin(), fV.end(), basisFE_ref_i.begin(), 0.0);


		};
	};




	SparseMatrixVectorResult assembledMatrix = AssembledMatrix(I, J, VAL, right_side_vector);

	return assembledMatrix;
}