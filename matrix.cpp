/* matrix.cpp */
#include "matrix.h"
#include "iostream"
#include "fstream"
#include <cstdio>
#include <cmath>
#include "image2D.h"
#include "complex"
#include <iterator>
   
void test1() {
	typedef complex<int> cx;
	Matrix<double> C(4,1);
	for (int i=0; i<4; i++) C(i) = i;
	Matrix< complex<double> > D = fft(C);
	Matrix< complex<double> > E = ifft(D);
	cout << C << endl;
	cout << D << endl;
	cout << E << endl;
}

void test2() {
	typedef complex<int> cx;
	ifstream ifs("dat");
	Matrix<double> C(ifs);
	ifs.close();
	cout << C << endl;
	cout << fft(C) << endl;
	cout << fft(fft(C),-1) << endl;
}

void test3() {
	bool a[] = { 0, 1, 0,  1, 0, 0,  1, 1, 1 };
	Matrix<bool> A(a, 3,3);
	bool b[] = { 0, 0, 0,  0, 1, 0,  0, 1, 0 };
	Matrix<bool> B(b, 3,3);
	Matrix<bool> C = erosion(A, B);
	cout << A << endl;
	cout << B << endl;
	cout << C << endl;

	bool x = 0;
	bool y = 1;
	bool z = x && y;
	cout << z;
}

void test4() {
	double a[] = { 25.5, 1, 5,  50, 80, 100,  0, 3, 200 };
	Matrix<double> A(a, 3,3);
	//savePNM(A, "test.pgm", true);
	cout << A << endl;
	Matrix<int> B = A;
	cout << B << endl;
}

int main(int argc, char* argv[]) {
	//test1();
	//test2();
	//test3();
	test4();

	return 0;
}

