#include <iostream>
#include "Matrix.h"
#include <iomanip>
#include <stdexcept>
#include <vector>
#include <fstream>
#include <string>
#include <chrono>

using namespace std;

Matrix::Matrix(int N, int M, double value = 0.0) {
	this->N = N;
	this->M = M;
	this->numbers = new double* [N];
	for (int i = 0; i < N; i++) {
		this->numbers[i] = new double[M];
		for (int j = 0; j < M; j++) {
			numbers[i][j] = value;
		}
	}
}

Matrix::Matrix(int N) {
	*this=Matrix(N, N);
}

Matrix::Matrix(int N, double value) {
	*this = Matrix(N, N);
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			(*this)[i][j] = value;
		}
	}
}

double* Matrix::operator[](int i) {
	return this->numbers[i];
}

double* Matrix::matrixMultiply(double* Second) {

	double* newMatrix = new double[this->getN()];
	int imax = this->getN();
	for (int i = 0; i < imax;i++) {
		newMatrix[i] = 0.0;
			for (int k = 0; k < imax; k++) {
				newMatrix[i] += (*this)[i][k] * Second[k];
			}
	}
	return newMatrix;
}

void Matrix::memoryMatrixMultiply(double* Second) {

	double* newMatrix = new double[this->getN()];
	int imax = this->getN();
	for (int i = 0; i < imax; i++) {
		newMatrix[i] = 0.0;
		for (int k = 0; k < imax; k++) {
			newMatrix[i] += (*this)[i][k] * Second[k];
		}
	}
	for (int i = 0; i < imax; i++) Second[i] = newMatrix[i];
	delete[] newMatrix;
}

Matrix trueMatrixMultiply(Matrix A, Matrix B) {

	Matrix newMatrix = Matrix(A.getN(), B.getM(), 0.0);
	for (int i = 0; i < A.getN(); i++) {
		for (int j = 0; j < B.getM(); j++) {
			for (int k = 0; k < A.getM(); k++) {
				newMatrix[i][j] += A[i][k] * B[k][j];
			}
		}
	}
	return newMatrix;
}

Matrix Matrix::matrixAdd(Matrix& Second) {

	Matrix newMatrix = Matrix(this->getN(), Second.getM(),0.0);
	for (int i = 0; i < newMatrix.getN(); i++) {
		for (int j = 0; j < newMatrix.getM(); j++) {
			newMatrix[i][j] = (*this)[i][j]+Second[i][j];
		}
	}
	return newMatrix;
}

Matrix Matrix::matrixSubstract(Matrix& Second) {

	Matrix newMatrix = Matrix(this->getN(), Second.getM(), 0.0);
	for (int i = 0; i < newMatrix.getN(); i++) {
		for (int j = 0; j < newMatrix.getM(); j++) {
			newMatrix[i][j] = (*this)[i][j] - Second[i][j];
		}
	}
	return newMatrix;
}

Matrix Matrix::Transpose() {
	Matrix newMatrix = Matrix(this->getM(), this->getN());
	for (int i = 0; i < newMatrix.getN(); i++) {
		for (int j = 0; j < newMatrix.getM(); j++) {
			newMatrix[i][j] = (*this)[j][i];
		}
	}
	return newMatrix;
}

double* Matrix::operator*(double* other) {
	//if (this->getM() != other.getN()) throw std::invalid_argument("Matrices dont have right sizes");
	return this->matrixMultiply(other);
}

Matrix Matrix::operator/(int other) {
	for (int i = 0; i < getN(); i++) numbers[i][i] = 1/numbers[i][i];
	return *this;
}

Matrix Matrix::operator+(Matrix& other) {
	if (this->getM() != other.getM() || this->getN() != other.getN()) throw std::invalid_argument("Matrices dont have right sizes");
	return this->matrixAdd(other);
}

Matrix Matrix::operator-(Matrix other) {
	if (this->getM() != other.getM() || this->getN() != other.getN()) throw std::invalid_argument("Matrices dont have right sizes");
	return this->matrixSubstract(other);
}


int Matrix::getN() {
	return this->N;
}

int Matrix::getM() {
	return this->M;
}

void print(Matrix matrix) {
	for (int i = 0; i < matrix.getN(); i++) {
		for (int j = 0; j < matrix.getM(); j++) {
			cout << matrix[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}

void clearFiles() {
	string FileName = "GaussCzas.csv";
	fstream FileStream;
	FileStream.open(FileName, ios::out);
	FileStream.close();

	FileName = "JacobiCzas.csv";
	FileStream;
	FileStream.open(FileName, ios::out);
	FileStream.close();

	FileName = "LUCzas.csv";
	FileStream;
	FileStream.open(FileName, ios::out);
	FileStream.close();

	FileName = "Rozmiary.csv";
	FileStream;
	FileStream.open(FileName, ios::out);
	FileStream.close();

	FileName = "resBGauss.csv";
	FileStream;
	FileStream.open(FileName, ios::out);
	FileStream.close();

	FileName = "resCGauss.csv";
	FileStream;
	FileStream.open(FileName, ios::out);
	FileStream.close();

	FileName = "resBJacobi.csv";
	FileStream;
	FileStream.open(FileName, ios::out);
	FileStream.close();

	FileName = "resCJacobi.csv";
	FileStream;
	FileStream.open(FileName, ios::out);
	FileStream.close();
}

//positive numbers shift to right, negative to the the left, 0 is the main diagonal
Matrix Diagonal(int N, double value, int diagonal) { 
	Matrix newMatrix = Matrix(N, N, 0.0);
	for (int i = 0; i < newMatrix.getN(); i++) {
		if(i-diagonal>=0 && i-diagonal < newMatrix.getN())newMatrix[i-diagonal][i] = value;
	}
	return newMatrix;
}

double* bVector(int N) {
	double* newMatrix = new double[N];
	for (int i = 0; i < N; i++) {
		newMatrix[i] = sin(i * (9));
	}
	return newMatrix;
}
void Matrix::addDiagonal(double value, int diagonal) {
	for (int i = 0; i < (*this).getN(); i++) {
		if (i - diagonal >= 0 && i - diagonal < (*this).getN())(*this)[i - diagonal][i] = value;
	}
}

Matrix getMainDiagonal(Matrix matrix) {
	Matrix newMatrix = Matrix(matrix.getN(), matrix.getM());
	for (int i = 0; i < newMatrix.getN(); i++) newMatrix[i][i] = matrix[i][i];
	return newMatrix;
}
Matrix Lower(Matrix matrix) {
	Matrix newMatrix = Matrix(matrix.getN(), matrix.getM());
	for (int i = 0; i < newMatrix.getN(); i++) {
		for (int j = 0; j < i; j++) {
			newMatrix[i][j] = matrix[i][j];
		}
	}
	return newMatrix;
}

Matrix Upper(Matrix matrix) {
	Matrix newMatrix = Matrix(matrix.getN(), matrix.getM());
	for (int i = 0; i < newMatrix.getN(); i++) {
		for (int j = matrix.getM(); j > i ; j--) {
			newMatrix[i][j] = matrix[i][j];
		}
	}
	return newMatrix;
}

void init(int precision) {
	cout << fixed;
	cout << setprecision(precision);
}
void Matrix::freeNumbers() {
	for (int i = 0; i < getN(); i++) {
		delete[](numbers[i]);
	}
	delete[](numbers);
}

double norm(double* matrix, int N) {
	double norm = 0.0;
	for (int i = 0; i < N; i++) {
			norm += matrix[i] * matrix[i];
	}
	norm = sqrt(norm);
	return norm;
}

void vectorSubstract(double* a, double* b, int N) {
	for (int i = 0; i < N; i++) a[i] = a[i] - b[i];
}

double* forwardSubstitution(Matrix L, double* b, bool clear) {
	double* y = new double[L.getN()];
	for (int i = 0; i < L.getN(); i++) y[i] = 1.0;
	double S = 0;
	for (int i = 0; i < L.getN(); i++) {
		S = 0;
		for (int j = 0; j < i; j++) {
			S += L[i][j] * y[j];
		}
		y[i] = (b[i] - S) / L[i][i];
	}
	if(clear)delete[](b);
	return y;
}

double* bakcwardSubstitution(Matrix U, double* b, bool clear) {
	int N = U.getN();
	double* y = new double[U.getN()];
	for (int i = N - 1; i >= 0; --i)
	{
		double S = 0;

		for (int j = i + 1; j < N; ++j) S += U[i][j] * y[j];

		y[i] = (b[i] - S) / U[i][i];
	}
	if(clear) delete[](b);
	return y;
}

Matrix copyMatrix(Matrix A) {
	Matrix newMatrix = Matrix(A.getN(), A.getM(), 0.0);
	for (int i = 0; i < A.getN(); i++) {
		for (int j = 0; j < A.getM(); j++) {
			newMatrix[i][j] = A[i][j];
		}
	}
	return newMatrix;
}
void Jacobi(Matrix A, double* b, int error, int time) {

	int N = A.getN();
	double* r = new double[A.getN()];
	for (int i = 0; i < A.getN(); i++) {
		r[i] = 1.0;
	}
	Matrix D = getMainDiagonal(A);
	Matrix L = Lower(A);
	Matrix U = Upper(A);
	Matrix Dinv = D / 1;
	Dinv = Dinv - Dinv - Dinv;
	double residual = 0.0;
	int iterations = 0;
	vector<double> residuals;
	cout << endl;
	auto start = chrono::high_resolution_clock::now();
	Matrix LU = L + U;
	do {
		r = (LU * r);
		vectorSubstract(r, b, N);
		Dinv.memoryMatrixMultiply(r);
		double* res = A * r;
		vectorSubstract(res, b, N);
		residual = norm(res, N);
		iterations++;
		residuals.push_back(residual);
		delete[] res;
	} while (residual > 1e-9 && residual < 1e+30);

	auto end = chrono::high_resolution_clock::now();
	auto difference = end - start;
	cout << "Czas obliczen Jacobi: " << chrono::duration<double, milli>(difference).count() <<  " ms" << endl;
	cout << "Liczba iteracji Jacobi : " << iterations << endl;
	cout << "Norma bledu rezydualnego Jacobi: " << residual << endl;

	string FileName;
	fstream FileStream;
	if (time) {
		FileName = "JacobiCzas.csv";
		FileStream;
		FileStream.open(FileName, ios_base::app);
		FileStream << chrono::duration<double, milli>(difference).count() << endl;
		FileStream.close();

	}
	if (error) {
		if (error == 1) {
			FileName = "resBJacobi.csv";
		}
		else FileName = "resCJacobi.csv";
		FileStream.open(FileName, ios_base::app);
		for (int i = 0; i < residuals.size(); i++) {
			FileStream << residuals[i] << endl;
		}
		FileStream.close();
	}


	delete[] (r);
}

void GaussSeidl(Matrix A, double* b, int error, int time) {
	int N = A.getN();
	double* r = new double[A.getN()];
	for (int i = 0; i < A.getN(); i++) {
		r[i] = 1.0;
	}
	Matrix D = getMainDiagonal(A);
	Matrix L = Lower(A);
	Matrix U = Upper(A);
	Matrix DL = D + L;
	DL = DL - DL - DL;
	vector<double> residuals;
	double residual = 0.0;
	int iterations = 0;
	cout << endl;
	auto start = chrono::high_resolution_clock::now();

	do {
		U.memoryMatrixMultiply(r);
		vectorSubstract(r, b, N);
		r = forwardSubstitution(DL, r,true);
		double* res = A * r;
		vectorSubstract(res, b, N);
		residual = norm(res, N);
		residuals.push_back(residual);
		iterations++;
		delete[] res;
	} while (residual > 1e-9 && residual < 1e+30);

	auto end = chrono::high_resolution_clock::now();
	auto difference = end - start;
	cout << "Czas obliczen Gauss-Seild: " << chrono::duration<double, milli>(difference).count() << " ms" << endl;
	cout << "Liczba iteracji Gauss-Seidl: " << iterations << endl;
	cout << "Norma bledu rezydualnego Gauss-Seidl: " << residual << endl;
	fstream FileStream;
	string FileName;

	if (time) {
		FileName = "GaussCzas.csv";
		FileStream;
		FileStream.open(FileName, ios_base::app);
		FileStream << chrono::duration<double, milli>(difference).count() << endl;
		FileStream.close();
	}
	if (error) {
		if (error == 1) {
			FileName = "resBGauss.csv";
		}
		else FileName = "resCGauss.csv";

		FileStream.open(FileName, ios_base::app);
		for (int i = 0; i < residuals.size(); i++) {
			FileStream << residuals[i] << endl;
		}
		FileStream.close();
	}

	L.freeNumbers();
	D.freeNumbers();
	U.freeNumbers();
	DL.freeNumbers();
	delete[] (r);
}

void LUFactorization(Matrix A, double* b, int time) {
	auto start = chrono::high_resolution_clock::now();
	int N = A.getN();
	Matrix U = copyMatrix(A);
	Matrix L = Matrix(A.getN(), A.getM(), 0.0);
	L.addDiagonal(1.0, MAIN_DIAGONAL);

	for (int i = 0; i < N - 1; ++i)		//columns
	{
		for (int j = i + 1; j < N; ++j)	//rows
		{
			L[j][i] = U[j][i] / U[i][i];

			for (int k = i; k < N; ++k)
				U[j][k] = U[j][k] - L[j][i] * U[i][k];
		}
	}

	double* y = forwardSubstitution(L, b, false);
	double* x = bakcwardSubstitution(U, y, true);
	double* temp = A * x;
	vectorSubstract(temp, b, N);

	auto end = chrono::high_resolution_clock::now();
	auto difference = end - start;
	cout << endl << "Czas obliczen LU: " << chrono::duration<double, milli>(difference).count() << " ms" << endl;
	cout << "Norma bledu rezydualnego LU: " << norm(temp, N) << endl;
	if (time) {
		string FileName = "LUCzas.csv";
		fstream FileStream;
		FileStream.open(FileName, ios_base::app);
		FileStream << chrono::duration<double, milli>(difference).count() << endl;
		FileStream.close();
	}
	delete[] (x);
	delete[] (temp);
}
