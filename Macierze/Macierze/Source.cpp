#include <iostream>
#include "matrix.h"
using namespace std;

int main() {

	Matrix Mac = Matrix(2,1.0);
	Matrix Mac2 = Matrix(2);

	for (int i = 0; i < Mac2.getN(); i++) {
		for (int j = 0; j < Mac2.getM(); j++) {
		Mac2[i][j] =6.0;
		}

	}

	Matrix Mac3 = Mac * Mac2;

	for (int i = 0; i < Mac3.getN(); i++) {
		for (int j = 0; j < Mac3.getM(); j++) {
			cout << Mac3[i][j] << " ";
		}
		cout << endl;
	}
	return 0;
}