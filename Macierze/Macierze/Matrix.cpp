#include <iostream>
#include "Matrix.h"
#include <stdexcept>
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

Matrix Matrix::matrixMultiply( Matrix& Second) {

	Matrix newMatrix = Matrix(this->getN(), Second.getM(),0);
	for (int i = 0; i < newMatrix.getN();i++) {
		for (int j = 0; j < newMatrix.getM(); j++) {
			for (int k = 0; k < Second.getN(); k++) {
				newMatrix[i][j] += (*this)[i][k] * Second[k][j];
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

Matrix Matrix::Transpose() {
	Matrix newMatrix = Matrix(this->getM(), this->getN());
	for (int i = 0; i < newMatrix.getN(); i++) {
		for (int j = 0; j < newMatrix.getM(); j++) {
			newMatrix[i][j] = (*this)[j][i];
		}
	}
	return newMatrix;
}

Matrix Matrix::operator*(Matrix& other) {
	if (this->getM() != other.getN()) throw std::invalid_argument("Matrices dont have right sizes");
	return this->matrixMultiply(other);
}

Matrix Matrix::operator+(Matrix& other) {
	if (this->getM() != other.getM() || this->getN() != other.getN()) throw std::invalid_argument("Matrices dont have right sizes");
	return this->matrixAdd(other);
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
}

//positive numbers shift to right, negative to the the left, 0 is the main diagonal
Matrix Diagonal(int N, double value, int diagonal) { 
	Matrix newMatrix = Matrix(N, N, 0.0);
	for (int i = 0; i < newMatrix.getN(); i++) {
		if(i-diagonal>=0 && i-diagonal < newMatrix.getN())newMatrix[i-diagonal][i] = value;
	}
	return newMatrix;
}

void Matrix::addDiagonal(double value, int diagonal) {
	for (int i = 0; i < (*this).getN(); i++) {
		if (i - diagonal >= 0 && i - diagonal < (*this).getN())(*this)[i - diagonal][i] = value;
	}
}
