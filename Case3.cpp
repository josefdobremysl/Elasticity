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
#include <chrono>	
#include <armadillo>


using namespace std;
using namespace Eigen;
using namespace arma;

using Time = std::chrono::steady_clock;
using ms = std::chrono::milliseconds;

void Case3() {

	// Definice typu pro řídkou matici s číselnými hodnotami typu double
	typedef Eigen::SparseMatrix<double> SparseMatrixd;

	std::string meshFile = "ctverec3.txt";


	const auto start1 = Time::now();
	SparseMatrixVectorResult stiffness_matrix = StiffnessMatrix(meshFile);
	const auto end1 = Time::now();
	const auto diff1 = std::chrono::duration_cast<ms>(end1 - start1).count();
	std::cout << " StiffnessMatrix time = " << diff1 << " ms" << "\n";
	SparseMatrix<double> K = stiffness_matrix.sparsematrix;
	vector<double> right_side_vector0 = stiffness_matrix.right_side_vector;
	int tagD1 = 100;

	const auto start2 = Time::now();
	SparseMatrixVectorResult Dirichlet1 = DirichletBoundaryCondition(meshFile, tagD1);
	const auto end2 = Time::now();
	const auto diff2 = std::chrono::duration_cast<ms>(end2 - start2).count();
	std::cout << " Dirichlet 1 time = " << diff2 << " ms" << "\n";

	SparseMatrix<double> sparsematrixDirichlet1 = Dirichlet1.sparsematrix;
	vector<double> bD1 = Dirichlet1.right_side_vector;

	int tagN2 = 200;
	vector<double> bN2 = NeumannBoundaryCondition(meshFile, tagN2);
	int tagD3 = 300;
	SparseMatrixVectorResult Dirichlet2 = DirichletBoundaryCondition(meshFile, tagD3);
	SparseMatrix<double> sparsematrixDirichlet3 = Dirichlet2.sparsematrix;
	vector<double> bD3 = Dirichlet2.right_side_vector;
	int tagD4 = 400;
	SparseMatrixVectorResult Dirichlet4 = DirichletBoundaryCondition(meshFile, tagD4);

	SparseMatrix<double> sparsematrixDirichlet4 = Dirichlet4.sparsematrix;
	vector<double> bD4 = Dirichlet4.right_side_vector;
	// Velikost matice
	int numRows = right_side_vector0.size();
	int numCols = right_side_vector0.size();

	vector<double> right_side_vector(numRows, 0);

	for (int i = 0; i < right_side_vector0.size(); i++) {

		right_side_vector[i] = right_side_vector0[i] + bD1[i] + bD3[i] + bN2[i] + bD4[i];
		//cout << right_side_vector[i] << "   " << right_side_vector0[i] << "   " << right_side_vectorN1[i] << endl;
	}

	SparseMatrix<double> A = K + sparsematrixDirichlet1 + sparsematrixDirichlet3 + sparsematrixDirichlet4;

	// Pravá strana soustavy b
	Eigen::VectorXd b = Eigen::Map<Eigen::VectorXd>(right_side_vector.data(), right_side_vector.size());
	//cout << "b" << b << endl;

	//////////////////////  ŘEŠENÍ SOUSTAVY  \\\\\\\\\\\\\\\\\\\\\\\\\\

	const auto start = Time::now();

	//// Inicializace øešení
	Eigen::VectorXd x(numRows);

	ConjugateGradient<SparseMatrix<double>, Eigen::Upper> solver;
	x = solver.compute(A).solve(b);
	string meshFile2 = meshFile;

	const auto end = Time::now();

	const auto diff = std::chrono::duration_cast<ms>(end - start).count();

	std::cout << "Solution time = " << diff << " ms" << "\n";

	/////////////////////////////////////////////

	//// Převod matice ze formátu Eigen sparse na formát Armadillo
	//sp_mat A_arm(A.rows(), A.cols()); // Vytvoření prázdné matice
	//for (int k = 0; k < A.outerSize(); ++k) {
	//	for (SparseMatrix<double>::InnerIterator it(A, k); it; ++it) {
	//		A_arm(it.row(), it.col()) = it.value(); // Nastavení hodnoty prvku
	//	}
	//}

	//// Paralelní řešení soustavy
	//
	//vec x2 = spsolve(A_arm, conv_to<vec>::from(right_side_vector));
	//

	/////////////////////////////////////////////

	const auto start3 = Time::now();

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

	const auto end3 = Time::now();
	const auto diff3 = std::chrono::duration_cast<ms>(end3 - start3).count();
	std::cout << "writing VTK file  time = " << diff3 << " ms" << "\n";

}