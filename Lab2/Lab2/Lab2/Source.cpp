#include "Header.h";
#include<fstream>
#include<vector>
#include<iostream>
#include<cmath>
#include<random>

using namespace std;

ifstream fin("in.txt");
ofstream fout("out.txt");

int n;
double a1, a2, a3, ast, aend, **X0;
double **X1, **X2;
int max_steps;
double min_error;
vector<double> a, b, c;

void thomas(
	const vector<double>& a,
	const vector<double>& b,
	const vector<double>& c,
	vector<double>& u,
	const vector<double>& d) {

	int n = u.size();
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

	b.push_back(ast);
	a.push_back(0);
	c.push_back(0);
	for (int i = 1; i <= n - 2; i++) {
		b.push_back(a2);
		a.push_back(a1);
		c.push_back(a3);
	}
	b.push_back(aend);
	a.push_back(0);
	c.push_back(0);
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

void constructBForX1(int j, vector<double>& B) {
	B.push_back(border(0, X0));
	for (int i = 1; i < n - 1; i++) {
		B.push_back(f(X0[i][j], X0[i+1][j], X0[i-1][j], X0[i][j+1], X0[i][j-1]));
	}
	B.push_back(border(1, X0));
}

void solveX1() {
	// solve for X1
	for (int j = 1; j < n - 1; j++) {
		// A x X(:,j) = Bj
		vector<double> B;
		constructBForX1(j, B);

		vector<double> res(n);
		thomas(a, b, c, res, B);

		for (int i = 0; i < n; i++) {
			X1[i][j] = res[i];
		}
	}
	for (int i = 0; i < n; i++) {
		X1[i][0] = X0[i][0];
		X1[i][n - 1] = X0[i][n - 1];
	}
}

void constructBForX2(int j, vector<double>& B) {
	B.push_back(border(2, X1));
	for (int i = 1; i < n - 1; i++) {
		B.push_back(f(X1[i][j], X1[i + 1][j], X1[i - 1][j], X1[i][j + 1], X1[i][j - 1]));
	}
	B.push_back(border(3, X1));
}

void solveX2() {
	// solve for X2
	for (int i = 1; i < n - 1; i++) {
		// A x X(i,:) = Bj
		vector<double> B;
		constructBForX2(i, B);

		vector<double> res(n);
		thomas(a, b, c, res, B);

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
		for (int j = 0; j < n; j++) {
			X0[i][j] = X2[i][j];
		}
	}
}

void solveSystems() {
	
	double error = std::numeric_limits<double>::max();
	int steps = 0;
	while (steps < max_steps && error > min_error) {
		solveX1();
		solveX2();
		error = computeError();
		steps += 1;
		cout << "Steps " << steps << endl;
		moveX2intoX0();
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

int parseInput(int argc, char* argv[]) {
	// expect n, min, max, maxsteps, minerr
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
	
	std::uniform_real_distribution<double> unif(min, max);
	std::default_random_engine re;
	
	generateData(unif, re);

	return 0;
}

int main(int argc, char* argv[]) {
	cout << "Started..." << endl;
	if (parseInput(argc, argv) != 0) {
		return 0;
	}
	generateThomasArrays();

	solveSystems();

	cout << "Done" << endl;

	//printMat();

	return 0;
}