#pragma once
class Matrix {

private:
	int N;
	int M;
	double** numbers;
public:
	Matrix(int N); //square Matrix
	Matrix(int N, double value); //square Matrix with values
	Matrix(int N,int M); //N x M Matrix
	Matrix(int N, int M,double value); //N x M Matrix with walues
	int getN();
	int getM();
	double* operator[](int);
	Matrix operator*(Matrix& other);
	Matrix operator+(Matrix& other);
	Matrix& matrixMultiply(Matrix& Second);
	Matrix& matrixAdd(Matrix& Second);
};
