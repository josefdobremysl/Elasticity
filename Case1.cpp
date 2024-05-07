#include "mesh_processor.h"
#include "Dirichlet.h"
#include "Neumann.h"
#include "Stiffnessmatrix.h"

#include <Eigen/Sparse>
#include <Eigen/Cholesky>
#include <Eigen/IterativeLinearSolvers>


#include <iostream>
#include <vector>
#include <tuple>
#include <Eigen/Sparse>
#include "vtk.h"

void Case1() {

	// Definice typu pro řídkou matici s číselnými hodnotami typu double
	typedef Eigen::SparseMatrix<double> SparseMatrixd;

	std::string meshFile = "ctverec.txt";

	SparseMatrixVectorResult stiffness_matrix = StiffnessMatrix(meshFile);
	SparseMatrix<double> K = stiffness_matrix.sparsematrix;
	vector<double> right_side_vector0 = stiffness_matrix.right_side_vector;
	int tagN1 = 100;
	vector<double> bN1 = NeumannBoundaryCondition(meshFile, tagN1);
	int tagN2 = 200;
	vector<double> bN2 = NeumannBoundaryCondition(meshFile, tagN2);
	int tagN3 = 300;
	vector<double> bN3 = NeumannBoundaryCondition(meshFile, tagN3);

	int tagD1 = 400;
	SparseMatrixVectorResult Dirichlet = DirichletBoundaryCondition(meshFile, tagD1);

	SparseMatrix<double> sparsematrixDirichlet = Dirichlet.sparsematrix;
	vector<double> bD = Dirichlet.right_side_vector;
	// Velikost matice
	int numRows = right_side_vector0.size();
	int numCols = right_side_vector0.size();

	vector<double> right_side_vector(numRows, 0);

	for (int i = 0; i < right_side_vector0.size(); i++) {

		right_side_vector[i] = right_side_vector0[i] + bD[i] + bN1[i] + bN2[i] + bN3[i];
		//cout << right_side_vector[i] << "   " << right_side_vector0[i] << "   " << right_side_vectorN1[i] << endl;
	}

	SparseMatrix<double> A = K + sparsematrixDirichlet;

	// Pravá strana soustavy b
	Eigen::VectorXd b = Eigen::Map<Eigen::VectorXd>(right_side_vector.data(), right_side_vector.size());
	//cout << "b" << b << endl;

	//// Inicializace øešení
	Eigen::VectorXd x(numRows);

	ConjugateGradient<SparseMatrix<double>, Eigen::Upper> solver;
	x = solver.compute(A).solve(b);
	string meshFile2 = meshFile;

	// Odvození jmena vystupniho souboru
	size_t dotPos = meshFile2.find_last_of(".");

	if (dotPos != std::string::npos && dotPos != meshFile2.size() - 1) {
		// Nahrazení části řetězce začínající po tečce novým řetězcem
		meshFile2 = meshFile2.substr(0, dotPos) + ".vtk";
	}
	else {
		meshFile2 += ".vtk";
	}
	VTKfile(meshFile2, meshFile, x);

	MatrixVectorResult result = GetTriangles(meshFile);


	//int nTri = result.integerResult;		// poèet trojúhelníkù
	//vector<int> I;
	//vector<int> J;
	//vector<double> VAL;
	int nNodes = result.matrix3.size();			// poèet bodù
	int nNodes2 = 2 * nNodes;
	vector<vector<double>>NodesM = result.matrix3;

	double Xi;
	double Yi;
	vector<double>UU(nNodes2, 0);

	for (int i = 0; i < nNodes; i++) {
		Xi = NodesM[i][0];
		Yi = NodesM[i][1];
		UU[i] = Xi;
		UU[i + nNodes] = 0;
	}

	Eigen::VectorXd UUU = Eigen::Map<Eigen::VectorXd>(UU.data(), UU.size());
	Eigen::VectorXd d = x - UUU;
	string meshFilechyba = "ctverec_chyba.vtk";

	VTKfile(meshFilechyba, meshFile, d);

}