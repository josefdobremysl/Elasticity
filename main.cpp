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
#include "Case1.h"
#include "Case2.h"
#include "Case3.h"
#include "Case4.h"
#include <chrono>	

using Time = std::chrono::steady_clock;
using ms = std::chrono::milliseconds;

int main() {

	const auto start = Time::now();

	Case3();

	const auto end = Time::now();
	const auto diff = std::chrono::duration_cast<ms>(end - start).count();
	std::cout << " Total time = " << diff << " ms" << "\n";
}