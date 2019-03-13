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