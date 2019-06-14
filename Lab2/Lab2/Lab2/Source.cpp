#include "Header.h";
#include<fstream>
#include<vector>
#include<iostream>
#include<cmath>
#include<random>
#include<chrono>

#include <mpi.h>
#include "stdio.h"
#include "stdlib.h"

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
		B[i] = f(X0[i][j], X0[i+1][j], X0[i-1][j], X0[i][j+1], X0[i][j-1]);
	}
	B[n - 1] = border(1, X0);
}

void constructBForX2(int j, double* B) {
	B[0] = border(2, X1);
	for (int i = 1; i < n - 1; i++) {
		B[i] = f(X1[i][j], X1[i + 1][j], X1[i - 1][j], X1[i][j + 1], X1[i][j - 1]);
	}
	B[n - 1] = border(3, X1);
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
	if (argc != 6) {
		cout << "invalid nr of args; expected <n> <min> <max> <maxsteps> <minerr> <soltype>" << endl;
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

void printOutput() {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			fout << X2[i][j] << " ";
		}
		fout << endl;
	}
}

void solveParallelX1() {
	cout << "Solving for X1" << endl;
	int size; MPI_Comm_size(MPI_COMM_WORLD, &size);

	// solve for X1
	for (int j = 1; j < n - 1; j++) {
		// A x X(:,j) = Bj
		double* B = new double[n];
		constructBForX1(j, B);
		int sendto = j % size;
		if (sendto == 0) {
			sendto += 1;
		}
		MPI_Send(B, n, MPI_DOUBLE, sendto, j, MPI_COMM_WORLD);
	}

	double* resp = new double[n];
	for (int j = 1; j < n - 1; j++) {
		int recvfr = j % size;
		if (recvfr == 0) {
			recvfr += 1;
		}
		MPI_Recv(resp, n, MPI_DOUBLE, recvfr, j, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		for (int i = 0; i < n; i++) {
			X1[i][j] = resp[i];
		}
	}
	
	for (int i = 0; i < n; i++) {
		X1[i][0] = X0[i][0];
		X1[i][n - 1] = X0[i][n - 1];
	}
}

void solveParallelX2() {
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

void broadcastData() {
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&a1, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&a2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&a3, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&ast, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&aend, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void worker_x1() {
	int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int size; MPI_Comm_size(MPI_COMM_WORLD, &size);

	cout << "Worker " << rank << " running " << endl;

	double* B = new double[n];
	double* res = new double[n];

	for (int j = 1; j < n - 1; j++) {
		int recvfr = j % size;
		if (recvfr == 0) {
			recvfr += 1;
		}
		if (recvfr == rank) {
			MPI_Recv(B, n, MPI_DOUBLE, 0, j, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			thomas(n, a, b, c, res, B);
			MPI_Send(res, n, MPI_DOUBLE, 0, j, MPI_COMM_WORLD);
		}
	}
}

void start_worker() {
	int stop = 0;
	while (stop == 0) {
		MPI_Recv(&stop, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		if (stop != 0) {
			break;
		}

		worker_x1();
		// worker x2
	}
}

int solveSystemWithMPI() {
	int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int size; MPI_Comm_size(MPI_COMM_WORLD, &size);

	if (rank == 0) {
		int stop = 0;
		double error = std::numeric_limits<double>::max();
		int steps = 0;
		while (steps < max_steps && error > min_error) {
			for (int i = 1; i < size; i++) {
				MPI_Send(&stop, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
			}

			steps += 1;
			cout << "Steps " << steps << endl;

			solveParallelX1();
			solveParallelX2();
			error = computeError();

			moveX2intoX0();
		}

		for (int i = 1; i < size; i++) {
			stop = 1;
			MPI_Send(&stop, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
		}
	}
	else {
		start_worker();
	}
	
	return 0;
}

int main(int argc, char* argv[]) {
	high_resolution_clock::time_point t1 = high_resolution_clock::now();

	// Initialize the MPI environment
	MPI_Init(&argc, &argv);

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	cout << "Started..." << endl;
	if (rank == 0) {
		//if (parseInput(argc, argv) != 0) {
		//	return 0;
		//}
		readData();
	}
	broadcastData();
	generateThomasArrays();

	solveSystemWithMPI();

	if (rank == 0) {
		printOutput();
		cout << "Done" << endl;
		high_resolution_clock::time_point t2 = high_resolution_clock::now();
		auto durationMilisec = duration_cast<milliseconds>(t2 - t1).count();
		cout << "Millisec " << durationMilisec << endl;
		
	}

	//printMat();
	MPI_Finalize();

	return 0;
}
