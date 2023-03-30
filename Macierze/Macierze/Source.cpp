#include <iostream>
#include <iomanip>
#include "matrix.h"
using namespace std;

int main() {
	//used to set output stream to output 2 decimal places
	cout << fixed;
	cout << setprecision(2);
	

	cout << "mam 1" << endl;
	Matrix Mac =  Matrix(100,100,1.0);
	cout << "mam 2" << endl;
	Matrix Mac2 = Matrix(100,100, 4);

	cout << "zaczynam mnozyc" << endl;
	Matrix Mac3 = Mac * Mac2;
	Matrix Mac4 = Diagonal(5, 1.0, MAIN_DIAGONAL);
	Mac4.addDiagonal(2.98, -1);
	print(Mac4);
	cout << "Done";
	return 0;
}