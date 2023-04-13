#include <iostream>
#include "matrix.h"
using namespace std;

int main() {
	//init(3);
	const int N[] = { 500,1000,2000,3000,5000 };

	for (int i = 0; i < sizeof(N) / sizeof(int); i++) {
		Matrix A = Matrix(N[i], N[i], 0.0);
		A.addDiagonal(12, MAIN_DIAGONAL);
		A.addDiagonal(-1.0, -1);
		A.addDiagonal(-1.0, -2);
		A.addDiagonal(-1.0, 1);
		A.addDiagonal(-1.0, 2);
		double* b = bVector(N[i]);

		cout << endl << "Matrix size: " << N[i] << endl;

		GaussSeidl(A, b);
		Jacobi(A, b);
		LUFactorization(A, b);
		A.freeNumbers();

	}
	return 0;
}