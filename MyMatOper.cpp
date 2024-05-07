#include <Eigen/Sparse>
#include <armadillo>

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <tuple>
#include <cmath>

using namespace Eigen;
using namespace std;
using namespace arma;

void printVecd(vector<double>& u) {
	for (int i = 0; i < u.size(); ++i) {
		cout << u[i] << " ";
	}
}