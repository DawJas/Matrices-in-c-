#include "Matrix.h"
#include <stdexcept>

Matrix::Matrix(int N, int M) {
	this->N = N;
	this->M = M;
	this->numbers = new double* [N];
	for (int i = 0; i < N; i++) {
		this->numbers[i] = new double[M];
	}
}

Matrix::Matrix(int N) {
	*this=Matrix(N, N);
}

Matrix::Matrix(int N, int M,double value) {
	*this = Matrix(N, M);
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) {
			(*this)[i][j] = value;
		}
	}
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

Matrix& Matrix::matrixMultiply( Matrix& Second) {

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

Matrix& Matrix::matrixAdd(Matrix& Second) {

	Matrix newMatrix = Matrix(this->getN(), Second.getM(),0.0);
	for (int i = 0; i < newMatrix.getN(); i++) {
		for (int j = 0; j < newMatrix.getM(); j++) {
			newMatrix[i][j] = (*this)[i][j]+Second[i][j];
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