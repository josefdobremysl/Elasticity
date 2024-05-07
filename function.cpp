#include "mesh_processor.h"
#include "function.h"
//#include "print.h"
//#include "quicksort.h"

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

const double pi = 3.14159265358979323846264338327950288;

double E = 100000;					// Materiálové vlastnosti
double sigma = 0.3;
double lambda = E * sigma / ((1 + sigma) * (1 - 2 * sigma));
double mu = E / (2 * (1 + sigma));

//double Function(double& X, double& Y) {
//	double f = 2 * pi * pi * (sin(pi * X) * sin(pi * Y));
//	return f;
//}

double f1;
double f2;

////////////////////////////////////////////////////////////////////////////////////
/////////////    Case1
//vector<double> FunctionV(double& X, double& Y) {
//	vector<double> fV(2, 0);
//	fV[0] = 0;
//	fV[1] = 0;
//	return fV;
//
//}
//
//
//
//double FunctionN1(double& X, double& Y, const int& tag, vector<double>& n) {
//	mat Tsigma = {
//	{ lambda + 2 * mu, 0},
//	{ 0, lambda}
//	};
//	double n_norm = (sqrt(n[0] * n[0] + n[1] * n[1]));
//	vec n_jednotkovy = {n[0]/n_norm,n[1] / n_norm };
//	vec ff1 = Tsigma * n_jednotkovy;
//	//ff1[0] = 0;
//	//ff1[1] = 0;
//	/*cout << "n: " << n[0]<<"   " << n[1] << endl;*/
//	if (tag == 100) {
//		f1 = ff1[0];
//		
//	}
//	if (tag == 200) {
//		f1 = ff1[0];
//	}
//	if (tag == 300) {
//		f1 = ff1[0];
//	}
//	if (tag == 400) {
//		f1 = 0;
//	}
//	return f1;
//}
//
//double FunctionN2(double& X, double& Y, const int& tag, vector<double>& n) {
//	mat Tsigma = {
//	{ lambda + 2 * mu, 0},
//	{ 0, lambda}
//	};
//	double n_norm = (sqrt(n[0] * n[0] + n[1] * n[1]));
//	vec n_jednotkovy = { n[0] / n_norm,n[1] / n_norm };
//	vec ff1 = Tsigma * n_jednotkovy;
//	//ff1[0] = 0;
//	//ff1[1] = 0;
//	
//	if (tag == 100) {
//		f2 = ff1[1];
//		//cout << "f2 " << f2 << endl;
//	}
//	if (tag == 200) {
//		f2 = ff1[1];
//	}
//	if (tag == 300) {
//		f2 = ff1[1];
//	}
//	if (tag == 400){
//		f2 = 0;
//	}
//	return f2;
//}
//////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////       case2
vector<double> FunctionV(double& X, double& Y) {
	vector<double> fV(2, 0);
	fV[0] = 0;// (mu + lambda);
	fV[1] =  -1 * (lambda + mu);
	return fV;

}



double FunctionN1(double& X, double& Y, const int& tag, vector<double>& n) {
	mat Tsigma = {
	{ lambda * (Y)+2 * mu * Y, mu * (X)},
	{ mu * (X), lambda * (Y) }
	};
	//cout << "X" << X << "Y<" << Y << endl;
	double n_norm = (sqrt(n[0] * n[0] + n[1] * n[1]));
	vec n_jednotkovy = { n[0] / n_norm,n[1] / n_norm };
	vec ff1 = Tsigma * n_jednotkovy;
	/*cout << "n: " << n[0]<<"   " << n[1] << endl;*/
	//ff1[0] = 0;
	//ff1[1] = 0;
	if (tag == 100) {
		f1 = 0;// ff1[0];

	}
	if (tag == 200) {
		f1 = ff1[0];
	}
	if (tag == 300) {
		f1 = Y*X;
	}
	if (tag == 400) {
		f1 = 0;
	}
	return f1;
}

double FunctionN2(double& X, double& Y, const int& tag, vector<double>& n) {
	mat Tsigma = {
	{ lambda * (Y ) + 2 * mu * Y, mu * ( X)},
	{ mu * (X), lambda*Y}
	};
	double n_norm = (sqrt(n[0] * n[0] + n[1] * n[1]));
	vec n_jednotkovy = { n[0] / n_norm,n[1] / n_norm };
	vec ff1 = Tsigma * n_jednotkovy;
	//ff1[0] = 0;
	//ff1[1] = 0;

	if (tag == 100) {
		f2 = 0;// ff1[1];
		//cout << "f2 " << f2 << endl;
	}
	if (tag == 200) {
		
		f2 = ff1[1];
	}
	if (tag == 300) {
		
		f2 =  0;
	}
	if (tag == 400) {
		f2 = 0;
	}
	return f2;
}
//////////////////////////////////////////////////////////////////////////////////


//bases functions
double fi1(double& x) {
	double fi1 = 1 - x;
	return fi1;
}
double fi2(double& x) {
	double fi2 = x;
	return fi2;
}

double dotprod(vector<double>& u, vector<double>& v) {
	double result = 0;
	for (int i = 0; i < u.size(); ++i) {
		result += u[i] * v[i];
	}
	return result;

}

//if (tag == 100) {			// Neumann
//	f1 = -pi * (sin(pi * X));
//}
//if (tag == 200) {		// Dirichlet
//	f1 = 0;
//}
//if (tag == 300) {		// Dirichlet
//	f1 = 0;
//}
//if (tag == 400) {		// Dirichlet
//	f1 = 0;
//}

