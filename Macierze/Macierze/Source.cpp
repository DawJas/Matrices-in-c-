#include <iostream>
#include "matrix.h"
using namespace std;

int main() {
	init(3);
	const int N = 933;
	Matrix A = Matrix(N, N, 0.0);
	A.addDiagonal(3,MAIN_DIAGONAL);
	A.addDiagonal(-1.0,-1);
	A.addDiagonal(-1.0,-2);
	A.addDiagonal(-1.0,1);
	A.addDiagonal(-1.0,2);
	double* b = bVector(N);

	GaussSeidl(A, b);
	Jacobi(A, b);



	A.freeNumbers();
	return 0;
}