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
vector<double> a, b, c;


void writeData() {
	fout << n << endl;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			fout << X0[i][j] << " ";
		}
		fout << endl;
	}
	fout << a1 << " " << a2 << " " << a3 << " " << ast << " " << aend << endl;
	fout << max_steps << endl;
	fout << min_error << endl;
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
	writeData();

	return 0;
}

int main(int argc, char* argv[]) {

	parseInput(argc, argv);

	return 0;
}

