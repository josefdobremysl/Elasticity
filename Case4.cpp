#include "mesh_processor.h"
#include "Dirichlet.h"
#include "Neumann.h"
#include "Stiffnessmatrix.h"
#include "MyMatOper.h"

#include <Eigen/Sparse>
#include <Eigen/Cholesky>
#include <Eigen/IterativeLinearSolvers>


#include <iostream>
#include <vector>
#include <tuple>
#include <Eigen/Sparse>
#include "vtk.h"

void Case4() {

	// Definice typu pro řídkou matici s číselnými hodnotami typu double
	typedef Eigen::SparseMatrix<double> SparseMatrixd;

	std::string meshFile = "ctverec2.txt";

	SparseMatrixVectorResult stiffness_matrix = StiffnessMatrix(meshFile);
	SparseMatrix<double> K = stiffness_matrix.sparsematrix;
	vector<double> right_side_vector0 = stiffness_matrix.right_side_vector;
	int tagD1 = 100;
	SparseMatrixVectorResult Dirichlet1 = DirichletBoundaryCondition(meshFile, tagD1);

	SparseMatrix<double> sparsematrixDirichlet1 = Dirichlet1.sparsematrix;
	vector<double> bD1 = Dirichlet1.right_side_vector;
	int tagD2 = 200;
	SparseMatrixVectorResult Dirichlet2 = DirichletBoundaryCondition(meshFile, tagD2);

	SparseMatrix<double> sparsematrixDirichlet2 = Dirichlet2.sparsematrix;
	vector<double> bD2 = Dirichlet2.right_side_vector;


	int tagD3 = 300;
	SparseMatrixVectorResult Dirichlet3 = DirichletBoundaryCondition(meshFile, tagD3);
	SparseMatrix<double> sparsematrixDirichlet3 = Dirichlet3.sparsematrix;
	vector<double> bD3 = Dirichlet3.right_side_vector;
	
	int tagD4 = 400;
	SparseMatrixVectorResult Dirichlet4 = DirichletBoundaryCondition(meshFile, tagD4);

	SparseMatrix<double> sparsematrixDirichlet4 = Dirichlet4.sparsematrix;
	vector<double> bD4 = Dirichlet4.right_side_vector;
	// Velikost matice
	int numRows = right_side_vector0.size();
	int numCols = right_side_vector0.size();

	vector<double> right_side_vector(numRows, 0);

	for (int i = 0; i < right_side_vector0.size(); i++) {

		right_side_vector[i] = right_side_vector0[i] + bD1[i] + bD2[i] + bD3[i] + bD4[i];
		//cout << right_side_vector[i] << "   " << right_side_vector0[i] << "   " << right_side_vectorN1[i] << endl;
	}


	SparseMatrix<double> A = K+ sparsematrixDirichlet1 + sparsematrixDirichlet2 + sparsematrixDirichlet3 + sparsematrixDirichlet4;

	// Pravá strana soustavy b
	Eigen::VectorXd b = Eigen::Map<Eigen::VectorXd>(right_side_vector.data(), right_side_vector.size());
	//cout << "b" << b << endl;
	//cout << "bD1: "<<endl;
	//printVecd(bD1);
	//cout << endl;
	//cout << "bD2: " << endl;
	//printVecd(bD2);
	//cout << endl;
	//cout << "bD3: " << endl;
	//printVecd(bD3);
	//cout << endl;
	//cout << "bD4: " << endl;
	//printVecd(bD4);
	//cout << endl;
	//cout << "riggtsidevec0: " << endl;
	//printVecd(right_side_vector0);
	//cout << endl;
	//cout << "riggtsidevec: " << endl;
	//printVecd(right_side_vector);

	//cout << A << endl;

	//// Inicializace øešení
	Eigen::VectorXd x(numRows);

	ConjugateGradient<SparseMatrix<double>, Eigen::Upper> solver;
	x = solver.compute(A).solve(b);
	string meshFile2 = meshFile;

	// Odvození jmena vystupniho souboru
	size_t dotPos = meshFile2.find_last_of(".");

	if (dotPos != std::string::npos && dotPos != meshFile2.size() - 1) {
		// Nahrazení části řetězce začínající po tečce novým řetězcem
		meshFile2 = meshFile2.substr(0, dotPos) + "_2.vtk";
	}
	else {
		meshFile2 += "_2.vtk";
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
		UU[i] = Xi*Yi;
		UU[i + nNodes] = 0;
	}

	Eigen::VectorXd UUU = Eigen::Map<Eigen::VectorXd>(UU.data(), UU.size());
	Eigen::VectorXd d = x - UUU;
	string meshFilechyba = "ctverec_chyba2.vtk";

	VTKfile(meshFilechyba, meshFile, d);

}