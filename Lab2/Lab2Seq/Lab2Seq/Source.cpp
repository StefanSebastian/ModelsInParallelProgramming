#include "Header.h";
#include<fstream>
#include<vector>
#include<iostream>
#include<cmath>
#include<random>
#include<chrono>

using namespace std;
using namespace std::chrono;

ifstream fin("in.txt");
ofstream fout("out.txt");

int n;
double a1, a2, a3, ast, aend, **X0;
double **X1, **X2;
int max_steps;
double min_error;
double *a, *b, *c;

void thomas(
	int n,
	double* a,
	double* b,
	double* c,
	double* u,
	double* d
) {
	vector<double> p(n);
	vector<double> q(n);

	// forward
	p[0] = -c[0] / b[0]; q[0] = d[0] / b[0];
	for (int i = 1; i < n; i++) {
		p[i] = -c[i] / (b[i] + a[i] * p[i - 1]);
		q[i] = (d[i] - a[i] * q[i - 1]) / (b[i] + a[i] * p[i - 1]);
	}

	// backward
	u[n - 1] = q[n - 1];
	for (int i = n - 2; i >= 0; i--) {
		u[i] = p[i] * u[i + 1] + q[i];
	}
}

double f(double p1, double p2, double p3, double p4, double p5) {
	return (p1 + p2 + p3 + p4 + p5) / 5;
}

double border(int k, double** mat) {
	return 1;
}

void generateThomasArrays() {
	cout << "generating thomas arrays" << endl;
	a = new double[n];
	b = new double[n];
	c = new double[n];
	a[0] = 0; b[0] = ast; c[0] = 0;
	for (int i = 1; i <= n - 2; i++) {
		a[i] = a1; b[i] = a2; c[i] = a3;
	}
	b[n - 1] = aend; c[n - 1] = 0; a[n - 1] = 0;
}

void readData() {
	cout << "Reading input data" << endl;
	fin >> n;
	X0 = new double*[n];
	X1 = new double*[n];
	X2 = new double*[n];
	for (int i = 0; i < n; i++) {
		X0[i] = new double[n];
		X1[i] = new double[n];
		X2[i] = new double[n];
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			fin >> X0[i][j];
		}
	}
	fin >> a1 >> a2 >> a3 >> ast >> aend;
	fin >> max_steps;
	fin >> min_error;
}

void generateData(uniform_real_distribution<double> unif, std::default_random_engine re) {
	X0 = new double*[n];
	X1 = new double*[n];
	X2 = new double*[n];
	for (int i = 0; i < n; i++) {
		X0[i] = new double[n];
		X1[i] = new double[n];
		X2[i] = new double[n];
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			// generate
			X0[i][j] = unif(re);
		}
	}

	a1 = unif(re);
	a2 = unif(re);
	a3 = unif(re);
	ast = unif(re);
	aend = unif(re);
}

void constructBForX1(int j, double* B) {
	B[0] = border(0, X0);
	for (int i = 1; i < n - 1; i++) {
		B[i] = f(X0[i][j], X0[i + 1][j], X0[i - 1][j], X0[i][j + 1], X0[i][j - 1]);
	}
	B[n - 1] = border(1, X0);
}

void solveX1() {
	cout << "Solving for X1" << endl;
	// solve for X1
	for (int j = 1; j < n - 1; j++) {
		// A x X(:,j) = Bj
		double* B = new double[n];
		constructBForX1(j, B);

		double* res = new double[n];
		thomas(n, a, b, c, res, B);

		for (int i = 0; i < n; i++) {
			X1[i][j] = res[i];
		}
	}
	for (int i = 0; i < n; i++) {
		X1[i][0] = X0[i][0];
		X1[i][n - 1] = X0[i][n - 1];
	}
}

void constructBForX2(int j, double* B) {
	B[0] = border(2, X1);
	for (int i = 1; i < n - 1; i++) {
		B[i] = f(X1[i][j], X1[i + 1][j], X1[i - 1][j], X1[i][j + 1], X1[i][j - 1]);
	}
	B[n - 1] = border(3, X1);
}

void solveX2() {
	cout << "Solving for X2 " << endl;
	// solve for X2
	for (int i = 1; i < n - 1; i++) {
		// A x X(i,:) = Bj
		double* B = new double[n];
		constructBForX2(i, B);

		double* res = new double[n];
		thomas(n, a, b, c, res, B);

		for (int j = 0; j < n; j++) {
			X2[i][j] = res[j];
		}
	}
	for (int i = 0; i < n; i++) {
		X2[0][i] = X1[0][i];
		X2[n - 1][i] = X1[n - 1][i];
	}
}

double computeError() {
	// error computation
	double error = std::numeric_limits<double>::min();
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			double curr = abs(X2[i][j] - X0[i][j]);
			if (error < curr) {
				error = curr;
			}
		}
	}
	cout << "Error is " << error << endl;
	return error;
}

void moveX2intoX0() {
	for (int i = 0; i < n; i++) {
		delete X0[i];
		X0[i] = X2[i];
		X2[i] = new double[n];
	}
}

void printMat() {
	cout << "----------------------X0-------------------" << endl;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			cout << X0[i][j] << " ";
		}
		cout << endl;
	}
	cout << "-----------------------X1------------------" << endl;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			cout << X1[i][j] << " ";
		}
		cout << endl;
	}
	cout << "------------------------X2-----------------" << endl;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			cout << X2[i][j] << " ";
		}
		cout << endl;
	}

	cout << "------------------------------------------" << endl;
}

void solveSystemsSequentially() {
	high_resolution_clock::time_point t1 = high_resolution_clock::now();

	cout << "Sequential solution " << endl;
	double error = std::numeric_limits<double>::max();
	int steps = 0;
	while (steps < max_steps && error > min_error) {
		steps += 1;
		cout << "Steps " << steps << endl;

		solveX1();
		solveX2();

		error = computeError();
		if (error > min_error) {
			moveX2intoX0();
		}
	}

	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	auto durationMilisec = duration_cast<milliseconds>(t2 - t1).count();
	cout << "Millisec " << durationMilisec << endl;
}

int parseInput(int argc, char* argv[]) {
	if (argc != 6) {
		cout << "invalid nr of args; expected <n> <min> <max> <maxsteps> <minerr>" << endl;
		return 1;
	}

	n = atoi(argv[1]);
	double min, max;
	sscanf(argv[2], "%lf", &min);
	sscanf(argv[3], "%lf", &max);
	max_steps = atoi(argv[4]);
	sscanf(argv[5], "%lf", &min_error);

	cout << "Running with n " << n << ", max_steps " << max_steps << ", min_error " << min_error << ", range " << min << ", " << max << endl;
	cout << "Sequential algorithm" << endl;

	std::uniform_real_distribution<double> unif(min, max);
	std::default_random_engine re;

	generateData(unif, re);

	return 0;
}

void printOutput() {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			fout << X2[i][j] << " ";
		}
		fout << endl;
	}
}


int main(int argc, char* argv[]) {

	cout << "Started..." << endl;
	//if (parseInput(argc, argv) != 0) {
	//	return 0;
	//}
	readData();
	generateThomasArrays();

	solveSystemsSequentially();
	
	printOutput();
	cout << "Done" << endl;

	return 0;
}

