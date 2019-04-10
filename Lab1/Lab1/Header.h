#pragma once

using namespace std;
#include<vector>
#include<string>

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
	const vector<double>& d);

void cyclic_reduction(
	vector<double>& a,
	vector<double>& b,
	vector<double>& c,
	vector<double>& u,
	vector<double>& d
);

void cyclic_reduction_omp(
	vector<double>& a,
	vector<double>& b,
	vector<double>& c,
	vector<double>& u,
	vector<double>& d,
	int num_threads
);

void cyclic_reduction_thr(
	vector<double>& a,
	vector<double>& b,
	vector<double>& c,
	vector<double>& u,
	vector<double>& d,
	int num_threads
);

void generate_thomas(int n,
	vector<double>& a,
	vector<double>& b,
	vector<double>& c,
	vector<double>& u,
	vector<double>& d);

void inner_loop_cyclic_red(int i, int n, double** A, double* F, int id, const vector<int>& indexes, int num_threads);

string algorithm_select(int alg);