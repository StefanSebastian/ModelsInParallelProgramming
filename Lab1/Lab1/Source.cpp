#include<iostream>
#include<fstream>
#include<vector>
#include<chrono>
#include<cmath>
#include<omp.h>

#include "Header.h"

using namespace std;
using namespace std::chrono;

ifstream fin("in.txt");
ofstream fout("out.txt");

/*
a, b, c diagonalspow
u - the unknowns
d - rhs matrix
*/
void thomas(
	const vector<double>& a,
	const vector<double>& b,
	const vector<double>& c,
	vector<double>& u,
	const vector<double>& d) {

	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	
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

	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	auto durationMicrosec = duration_cast<microseconds>(t2 - t1).count();
	auto durationMilisec = duration_cast<milliseconds>(t2 - t1).count();
	cout << "Microsec " << durationMicrosec << endl;
	cout << "Millisec " << durationMilisec << endl;
}

/*
a, b, c diagonals
u - the unknowns
d - rhs matrix
*/
void cyclic_reduction(
	const vector<double>& a,
	const vector<double>& b,
	const vector<double>& c,
	vector<double>& u,
	const vector<double>& d) {

	int n = u.size();
	int r = log2(n);

	// memory allocation
	double *F = new double[n];
	double **A = new double*[n];
	for (int i = 0; i < n; i++) {
		A[i] = new double[n];
		for (int j = 0; j < n; j++) {
			A[i][j] = 0;
		}
	}
	for (int i = 0; i < n; i++) {
		F[i] = d[i];
		A[i][i] = b[i];
		if (i - 1 >= 0) {
			A[i][i - 1] = a[i];
		}
		if (i + 1 < n) {
			A[i][i + 1] = c[i];
		}
	}

	omp_set_num_threads(4);

	// cyclic reduction
	for (int i = 0; i < log2(n + 1) - 1; i++) {
		int kl = pow(2, i + 1);

		#pragma omp parallel 
		{
			#pragma omp for
			for (int j = kl - 1; j < n; j = j + kl) {
				int offset = pow(2, i);
				int idx1 = j - offset;
				int idx2 = j + offset;

				double alpha = A[j][idx1] / A[idx1][idx1];
				double gamma = A[j][idx2] / A[idx2][idx2];

				for (int k = 0; k < n; k++) {
					A[j][k] -= (alpha * A[idx1][k] + gamma * A[idx2][k]);
				}
				F[j] -= (alpha * F[idx1] + gamma * F[idx2]);
			}
		}
	}

	// back subst
	int idx = (n - 1) / 2;
	u[idx] = F[idx] / A[idx][idx];

	for (int i = log2(n + 1) - 2; i >= 0; i--) {
		for (int j = pow(2, i + 1) - 1; j < n; j += pow(2, i + 1)) {
			int offset = pow(2, i);
			int idx1 = j - offset;
			int idx2 = j + offset;

			u[idx1] = F[idx1];
			u[idx2] = F[idx2];
			for (int k = 0; k < n; k++) {
				if (k != idx1) {
					u[idx1] -= A[idx1][k] * u[k];
				}
				if (k != idx2) {
					u[idx2] -= A[idx2][k] * u[k];
				}
			}
			u[idx1] = u[idx1] / A[idx1][idx1];
			u[idx2] = u[idx2] / A[idx2][idx2];
		}
	}
}

void read_vector(vector<double>& v, int n) {
	for (int i = 0; i < n; i++) {
		fin >> v[i];
	}
}

void solve_from_file() {
	int n;
	fin >> n;
	vector<double> a(n);
	vector<double> b(n);
	vector<double> c(n);
	vector<double> d(n);
	read_vector(a, n);
	read_vector(b, n);
	read_vector(c, n);
	read_vector(d, n);

	vector<double> u(n); // solution
	cyclic_reduction(a, b, c, u, d);
	for (int i = 0; i < n; i++) {
		fout << u[i] << " ";
	}
}

void generate_benchmark() {
	int n = 1000000;
	vector<double> a(n);
	vector<double> b(n);
	vector<double> c(n);
	vector<double> d(n);
	vector<double> u(n);
	generate_thomas(n, a, b, c, u, d);
	fout << n << std::endl;
	for (int i = 0; i < n; i++) {
		fout << a[i] << " ";
	}
	fout << std::endl;
	for (int i = 0; i < n; i++) {
		fout << b[i] << " ";
	}
	fout << std::endl;
	for (int i = 0; i < n; i++) {
		fout << c[i] << " ";
	}
	fout << std::endl;
	for (int i = 0; i < n; i++) {
		fout << d[i] << " ";
	}
	fout << std::endl;
}

void generate_benchmark2() {
	int n = pow(2, 10) - 1;
	vector<double> a(n);
	vector<double> b(n);
	vector<double> c(n);
	vector<double> d(n);
	vector<double> u(n);
	generate_thomas(n, a, b, c, u, d);
	cyclic_reduction(a, b, c, u, d);
	for (int i = 0; i < u.size(); i++) {
		cout << u[i] << " ";
	}
}

int main() {
	//solve_from_file();
	//generate_benchmark();
	generate_benchmark2();
	return 0;
}
