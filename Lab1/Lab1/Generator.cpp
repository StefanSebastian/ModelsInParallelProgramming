#include "Header.h"

using namespace std;

/*
Generates a thomas system with size n
solution is array of 1, 1, 1 .. 
a, b, c are filled from 1 to n
*/
void generate_thomas(int n,
	vector<double>& a,
	vector<double>& b,
	vector<double>& c,
	vector<double>& u,
	vector<double>& d) {
	
	// sol vector
	for (int i = 0; i < n; i++) {
		u[i] = 1;
	}

	a[0] = 0;
	int ai = 1; int bi = 0; int ci = 0;
	for (int i = 1; i <= 4 + 3 * (n - 2); i++) { // 2 on the 1st row, 2 on the last, 3 on the others
		if (i % 3 == 1) {
			b[bi] = i; bi++;
		}
		else if (i % 3 == 2) {
			c[ci] = i; ci++;
		}
		else {
			a[ai] = i; ai++;
		}
	}
	c[ci] = 0;

	d[0] = 3; // first row
	int di = 1; // intermediary rows
	for (int i = 3; i <= 3 * (n - 2); i += 3) {
		d[di] = 3 * i + 3; di++;
	}
	d[di] += 3 * (n - 1) * 2 + 1; //last row
}