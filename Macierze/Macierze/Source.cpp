#include <iostream>
#include <fstream>
#include <string>
#include "matrix.h"
using namespace std;

int main() {
	const int Size = 933;
	const int N[] = {100, 500, 1000, 2000, 3000};

	clearFiles();

	cout << "---------------ZADANIE B---------------" << endl;

	Matrix A = Matrix(Size, Size, 0.0);
	A.addDiagonal(12, MAIN_DIAGONAL);
	A.addDiagonal(-1.0, -1);
	A.addDiagonal(-1.0, -2);
	A.addDiagonal(-1.0, 1);
	A.addDiagonal(-1.0, 2);
	double* b = bVector(Size);

	GaussSeidl(A, b, 0);
	Jacobi(A, b, 0);
	A.freeNumbers();


	cout << endl << "---------------ZADANIE C---------------" << endl;



	A = Matrix(Size, Size, 0.0);
	A.addDiagonal(3, MAIN_DIAGONAL);
	A.addDiagonal(-1.0, -1);
	A.addDiagonal(-1.0, -2);
	A.addDiagonal(-1.0, 1);
	A.addDiagonal(-1.0, 2);

	GaussSeidl(A, b, 0);
	Jacobi(A, b, 0);


	cout << endl << "---------------ZADANIE D---------------" << endl;

	LUFactorization(A, b);

	for (int i = 0; i < sizeof(N) / sizeof(int); i++) {

		string FileName = "Rozmiary.csv";
		fstream FileStream;
		FileStream.open(FileName, ios_base::app);
		FileStream << N[i] << endl;
		FileStream.close();

		Matrix A = Matrix(N[i], N[i], 0.0);
		A.addDiagonal(12, MAIN_DIAGONAL);
		A.addDiagonal(-1.0, -1);
		A.addDiagonal(-1.0, -2);
		A.addDiagonal(-1.0, 1);
		A.addDiagonal(-1.0, 2);
		double* b = bVector(N[i]);

		cout << endl << "-------------------Matrix size: " << N[i] << "-------------------" << endl;

		GaussSeidl(A, b,0);
		Jacobi(A, b,0);
		LUFactorization(A, b);
		A.freeNumbers();
		delete[](b);

	}

	return 0;
}