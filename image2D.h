/** image2D.h
 * @author Toshihide Shimayama <tshm@csc.jp>
 */
#ifndef _IMAGE2DH_

#define _IMAGE2DH_
#include "matrix.h"
#include <complex>
#include <fstream>
using namespace std;
using namespace matrix;

// Image processing module
struct rgb {
	char r;
	char g;
	char b;
};

/** PGM save function
 * @param A input image
 * @param fname output filename
 * @param binary format flag
 */
void savePNM(const Matrix<char>& A, const char* fname, bool binary=false) {
	ofstream ofs(fname);
	if (binary) {
		ofs << "P5\n";
		ofs << A.cols() << " " << A.rows() << "\n255" << endl;
		for (int i=0; i<A.rows(); i++) {
			for (int j=0; j<A.cols(); j++) {
				ofs.put(A(i,j));
			}
		}
	} else {
		ofs << "P2\n";
		ofs << A.cols() << " " << A.rows() << "\n255" << endl;
		for (int i=0; i<A.rows(); i++) {
			for (int j=0; j<A.cols(); j++) {
				ofs << A(i,j) << " ";
			}
			ofs << "\n";
		}
	}
}

/** convolution filtering
 * @param A Original Image.
 * @param B Convolution kernel.
 */
template<class T> const Matrix<T> conv(const Matrix<T>& A, const Matrix<T>& B) {
	int Arow = A.rows();
	int Acol = A.cols();
	int Brow = B.rows();
	int Bcol = B.cols();
	Matrix<T> res(0.0, Arow, Acol);
	int p0 = (int)(Brow / 2);
	int q0 = (int)(Bcol / 2);
	for (int i=0; i<Arow; i++) {
		for (int j=0; j<Acol; j++) {
			T tmp=0;
			for (int p=0; p<Brow; p++) {
				int ii = i + p - p0;
				if (ii<0 || ii>=Arow) continue;
				for (int q=0; q<Bcol; q++) {
					int jj = j + q - q0;
					if (jj<0 || jj>=Acol) continue;
					tmp += A(ii,jj) * B(p,q);
				}
			}
			res(i,j) = tmp;
		}
	}
	return res;
}

/** bit reverse function (used by fft)
 */
inline unsigned int __bitReverse(unsigned int x, int nbits) {
  int n = 0;
  for (int i=0; i<nbits; i++) {
    n <<= 1;
    n |= (x & 1);
    x >>= 1;
  }
  return n;
}
/** Fourier Transform
 * @param org input image
 * @param rev direction rev=1(default) does forward FFT, rev=-1 does inverse FFT.
 * @return FFT(A)
 */
template<class T>
const Matrix< complex<double> > fft(const Matrix<T>& org, int rev=1) {
	const double PI = 3.1415926536;
	rev = rev>0 ? 1 : -1;
	int rowbits = (int)ceil( log2((double)org.rows()) );
	int colbits = (int)ceil( log2((double)org.cols()) );
	typedef complex<double> cx;
  unsigned int nrows = 1 << rowbits;
  unsigned int ncols = 1 << colbits;
	Matrix<cx> res(nrows, ncols);
  const cx I(0, rev);
	for (int i=0; i < org.rows(); ++i) {
		for (int j=0; j < org.cols(); ++j) {
			res(__bitReverse(i, rowbits), __bitReverse(j, colbits)) = org(i, j);
		}
	}
	// rows (the first dimension)
	for (int jj=0; jj<(int)ncols; jj++) {
		for (int s=1; s <= rowbits; ++s) {
			unsigned int m = 1 << s;
			unsigned int m2 = m >> 1;
			cx w(1, 0);
			cx wm = exp(-I * (PI / m2));
			for (unsigned int j=0; j < m2; ++j) {
				for (unsigned int k=j; k < nrows; k += m) {
					cx t = w * res(k + m2, jj);
					cx u = res(k, jj);
					res(k, jj) = u + t;
					res(k + m2, jj) = u - t;
				}
				w *= wm;
			}
		}
	}
	// cols (the second dimension)
	for (int jj=0; jj<(int)nrows; jj++) {
		for (int s=1; s <= colbits; ++s) {
			unsigned int m = 1 << s;
			unsigned int m2 = m >> 1;
			cx w(1, 0);
			cx wm = exp(-I * (PI / m2));
			for (unsigned int j=0; j < m2; ++j) {
				for (unsigned int k=j; k < ncols; k += m) {
					cx t = w * res(jj, k + m2);
					cx u = res(jj, k);
					res(jj, k) = u + t;
					res(jj, k + m2) = u - t;
				}
				w *= wm;
			}
		}
	}
	return rev==1 ? res : res/cx(ncols*nrows,0);
}

/** Inverse Fourier Transform
 * @param org input image
 * @param rev direction rev=1(default) does inverse FFT, rev=-1 does forward FFT.
 * @return IFFT(A)
 */
template<class T>
const Matrix< complex<double> > ifft(const Matrix<T>& org, int rev=1) {
	return fft(org, -1);
}
/** dilation
 * @param A input bool image
 * @param B kernel
 */
const Matrix<bool> dilation(const Matrix<bool>& A, const Matrix<bool>& B) {
	int Arow = A.rows();
	int Acol = A.cols();
	int Brow = B.rows();
	int Bcol = B.cols();
	Matrix<bool> res(Arow, Acol);
	int p0 = (int)(Brow / 2);
	int q0 = (int)(Bcol / 2);
	for (int i=0; i<Arow; i++) {
		for (int j=0; j<Acol; j++) {
			if (!A(i,j)) continue;
			for (int p=0; p<Brow; p++) {
				int ii = i + p - p0;
				if (ii<0 || ii>=Arow) continue;
				for (int q=0; q<Bcol; q++) {
					int jj = j + q - q0;
					if (jj<0 || jj>=Acol) continue;
					res(ii,jj) |= A(ii,jj) || B(p,q);
				}
			}
		}
	}
	return res;
}

/** erosion
 * @param A input bool image
 * @param B kernel
 */
const Matrix<bool> erosion(const Matrix<bool>& A, const Matrix<bool>& B) {
	int Arow = A.rows();
	int Acol = A.cols();
	int Brow = B.rows();
	int Bcol = B.cols();
	Matrix<bool> res(1, Arow, Acol);
	int p0 = (int)(Brow / 2);
	int q0 = (int)(Bcol / 2);
	for (int i=0; i<Arow; i++) {
		for (int j=0; j<Acol; j++) {
			if (A(i,j)) continue;
			for (int p=0; p<Brow; p++) {
				int ii = i + p - p0;
				if (ii<0 || ii>=Arow) continue;
				for (int q=0; q<Bcol; q++) {
					int jj = j + q - q0;
					if (jj<0 || jj>=Acol) continue;
					res(ii,jj) &= A(ii,jj) && !B(p,q);
				}
			}
		}
	}
	return res;
}
/** opening
 * @param A input bool image
 * @param B kernel
 */
const Matrix<bool> opening(const Matrix<bool>&A, const Matrix<bool>& B) {
	return dilation(erosion(A, B), B);
}
/** closing
 * @param A input bool image
 * @param B kernel
 */
const Matrix<bool> closing(const Matrix<bool>&A, const Matrix<bool>& B) {
	return erosion(dilation(A, B), B);
}
#endif /* _IMAGE2DH_ */
// vim:syntax=cpp.doxygen:
