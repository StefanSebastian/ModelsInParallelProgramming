#include<iostream>
#include<fstream>
#include<vector>
using namespace std;

ifstream fin("in.txt");
ofstream fout("out.txt");

/*
a, b, c diagonals
u - the unknowns
d - rhs matrix
*/
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
	thomas(a, b, c, u, d);
	for (int i = 0; i < n; i++) {
		fout << u[i] << " ";
	}
}

void main() {
	solve_from_file();
}
