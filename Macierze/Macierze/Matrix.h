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
	double* operator[](int);
	Matrix operator*(Matrix& other);
	Matrix operator+(Matrix& other);
	Matrix matrixMultiply(Matrix& Second);
	Matrix matrixAdd(Matrix& Second);
	Matrix Transpose();
	void addDiagonal(double value, int diagonal);

};
Matrix Diagonal(int N, double value, int diagonal);
void print(Matrix matrix);