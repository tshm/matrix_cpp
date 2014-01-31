// testMatrix.cpp
#include "matrix.h"
#include <unittest++/UnitTest++.h>
#include <unittest++/TestReporterStdout.h>

using namespace matrix;
struct fixture {
	Matrix<double> A, B, C;
	Matrix<double> x, y, z;

	fixture() {
		// Matrix Fixtures
		double a[] = {0.0, 1.5, 2.0, 1.3};
		double b[] = {0.0, 1.5, 2.0, 1.3};
		double c[] = {0.0, 1.5, 2.0, 1.3};
		A = Matrix<double>(a,2,2);
		B = Matrix<double>(b,2,2);
		C = Matrix<double>(c,2,2);
		// Vector Fixtures
		double xx[] = {1.0, 2.0};
		double yy[] = {1.0, 2.0};
		double zz[] = {1.0, 2.0};
		x = Matrix<double>(xx,2,1);
		y = Matrix<double>(yy,2,1);
		z = Matrix<double>(zz,2,1);
	}
	~fixture() {}
};

SUITE(base) {
	TEST_FIXTURE(fixture,testBase) {
		CHECK_EQUAL(A, B);
		CHECK_EQUAL(B, C);
	}
}

TEST(external_allocation) {
	double tmp[] = {1,2,3,4};
	double *x[2];
	x[0] = &(tmp[0]);
	x[1] = &(tmp[2]);
	Matrix<double> X(x, 2,2, true);  //reference delegation
	Matrix<double> Y(X, true);  //reference delegation
	Matrix<double> Z(X);  //reference delegation
	CHECK_EQUAL(X(0,0), 1);
	CHECK_EQUAL(Y(0,0), 1);
	CHECK_EQUAL(Z(0,0), 1);
	tmp[0] = 5;
	CHECK_EQUAL(X(0,0), 5);
	CHECK_EQUAL(Y(0,0), 5);
	CHECK_EQUAL(Z(0,0), 1);
}

TEST_FIXTURE(fixture,pointer_retriever) {
	double orig = A(0,0);
	double **a = A.getPointer();
	CHECK_EQUAL(**a, orig);
	**a = 3.0;
	CHECK_EQUAL(**a, 3.0);
	CHECK_EQUAL(A(0,0), 3.0);
}

SUITE(matrix_math) {
	TEST_FIXTURE(fixture,sum) {
		Matrix<double> TMP = A + B;
		CHECK_EQUAL(TMP(0,0), 0.0);
		CHECK_EQUAL(TMP(0,1), 3.0);
		CHECK_EQUAL(TMP(1,0), 4.0);
		CHECK_EQUAL(TMP(1,1), 2.6);
	}
	TEST_FIXTURE(fixture,subtract) {
		Matrix<double> TMP = A - B;
		CHECK_EQUAL(TMP(0,0), 0.0);
		CHECK_EQUAL(TMP(0,1), 0.0);
		CHECK_EQUAL(TMP(1,0), 0.0);
		CHECK_EQUAL(TMP(1,1), 0.0);
	}
	TEST_FIXTURE(fixture,transpose) {
		Matrix<double> TMP = A.t();
		CHECK_EQUAL(TMP(0,0), 0.0);
		CHECK_EQUAL(TMP(0,1), 2.0);
		CHECK_EQUAL(TMP(1,0), 1.5);
		CHECK_EQUAL(TMP(1,1), 1.3);
	}
}

SUITE(vector_math) {
	TEST_FIXTURE(fixture,multiply) {
		Matrix<double> TMP = A * x;
		CHECK_EQUAL(TMP(0), 3.0);
		CHECK_EQUAL(TMP(1), 4.6);
	}
	TEST_FIXTURE(fixture,sum) {
		Matrix<double> TMP = x + y;
		CHECK_EQUAL(TMP(0), 2.0);
		CHECK_EQUAL(TMP(1), 4.0);
	}
	TEST_FIXTURE(fixture,subtract) {
		Matrix<double> TMP = x - y;
		CHECK_EQUAL(TMP(0), 0.0);
		CHECK_EQUAL(TMP(1), 0.0);
	}
}

SUITE(matrix_analysis) {
	TEST_FIXTURE(fixture,determinant) {
		CHECK_EQUAL(A.det(), 3.0);
		A(0,0) = A(0,1) = A(1,0) = A(1,1) = 3;  //singular matrix
		CHECK_CLOSE(A.det(), 0.0, 0.00001);
	}
}

int main(int argc, char* const argv[]) {
	if (argc==1) {
		return UnitTest::RunAllTests();
	} else {
		UnitTest::TestReporterStdout reporter;
		return UnitTest::RunAllTests(reporter,
				UnitTest::Test::GetTestList(), argv[1], 0);
	}
}
