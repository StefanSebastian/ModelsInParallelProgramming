#pragma once

using namespace std;
#include<vector>

/*
a, b, c diagonals
u - the unknowns
d - rhs matrix
*/
void thomas(
	int n,
	double* a,
	double* b,
	double* c,
	double* u,
	double* d
);