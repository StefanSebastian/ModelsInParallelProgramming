#pragma once

using namespace std;
#include<vector>

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
	const vector<double>& a,
	const vector<double>& b,
	const vector<double>& c,
	vector<double>& u,
	const vector<double>& d
);

void parallel_cyclic_reduction(
	vector<double>& a,
	vector<double>& b,
	vector<double>& c,
	vector<double>& u,
	vector<double>& d
);

void generate_thomas(int n,
	vector<double>& a,
	vector<double>& b,
	vector<double>& c,
	vector<double>& u,
	vector<double>& d);