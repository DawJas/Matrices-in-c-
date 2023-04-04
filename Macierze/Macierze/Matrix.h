#pragma once
#define MAIN_DIAGONAL 0
class Matrix {
	// N - [x][]
	// M - [][x]
private:
	int N;
	int M;
	double** numbers;
public:
	Matrix(int N); //square Matrix
	Matrix(int N, double value); //square Matrix with values
	Matrix(int N, int M,double value); //N x M Matrix with walues
	int getN();
	int getM();
	double* operator[](int); //allows M[i][j] notation instead of M.numbers[i][j];
	double* operator*(double* other);
	Matrix operator/(int empty); // ONLY FOR DIGONALS!!!
	Matrix operator+(Matrix& other);
	Matrix operator-(Matrix other);
	double* matrixMultiply(double* Second);
	Matrix matrixAdd(Matrix& Second);
	Matrix matrixSubstract(Matrix& Second);
	Matrix Transpose();
	void addDiagonal(double value, int diagonal);
	void freeNumbers();

};
Matrix Diagonal(int N, double value, int diagonal);
double* bVector(int N);
Matrix getMainDiagonal(Matrix matrix);
Matrix Lower(Matrix matrix);
Matrix Upper(Matrix matix);
void print(Matrix matrix);
void init(int precision);
double norm(double* matrix, int N);
void vectorSubstract(double* a, double* b, int N);
void Jacobi(Matrix A, double* b);
void GaussSeidl(Matrix A, double* b);
void LUFactorization(Matrix A, double* b);