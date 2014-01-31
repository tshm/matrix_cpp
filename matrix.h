/** @file matrix.h
 * @author Toshihide Shimayama <tshm@csc.jp>
 */
#ifndef _MATRIX_
#define _MATRIX_

#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <cstdarg>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <typeinfo>

namespace matrix {
/** boundary check for accessors. */
const int    _IO_BUFFER_LENGTH_=50000;    // fgets() buffer for file reading
const double _INFINITESIMAL_SIZE_=1.0e-30;// threshold size for == operator
/** rassert : runtime_assert macro.
 * Runtime error checking followed by abortion.
 */
#define rassert(cnd,msg) if (!(cnd)) {fprintf(stderr, msg); abort();}
using namespace std;

/** T valued fixed dimenstion Matrix class.
 * Matrix class of T values.  Matrix related operations should 
 *  be implemented within this class.
 * @note CAUTION: indices start from 0, not 1.
 * @todo move the implementation of inner classes into separate cpp file.
 */
template<class T> class Matrix {
	protected:
		T **m;  ///< 1D T array
		int row;  ///< number of rows
		int col;  ///< number of cols
		bool alloc; ///< external allocation?
	private:
		/** memory allocation.
		 */
		void allocate() {
			int vsize = sizeof(T)*row;
			m = (T**)malloc(sizeof(T*)*col);
			assert(m);
			for (int i=0; i<col; i++) {
				m[i] = (T*)malloc(vsize);
				assert(m[i]);
			}
		}
		void de_allocate() {
			for (int i=0; i<col; i++) free(m[i]);
			free(m);
		}
		inline bool bchk(int i, int j=0) {
			return i>=0 && i<col && j>=0 && j<row;
		}
	public:
		/** T array valued constructor.  */
		Matrix(
				T* data, ///< 1D T array.
				///< (i,j)-th value is accessed with arr[ i*col + j ];
				int row,  ///< row length.
				int col   ///< column length.
				): row(row), col(col), alloc(false) {
			allocate();
			for (int i=0; i<row; i++)
				for (int j=0; j<col; j++) operator()(i,j) = data[i*col + j];
		}
		/** T 2D-array valued constructor.
		 * @note: treated as a collection of column-vectors.
		 * It means i-th row and j-th column is accessed by arr[j][i].
		 */
		Matrix(
				T** data, ///< 2D T array.
				///< (i,j)-th value is accessed with arr[j][i].
				int row,  ///< row length.
				int col,  ///< column length.
				bool alloc=false
				): row(row), col(col), alloc(alloc) {
			if (alloc) {
				m = data;
			} else {
				allocate();
				for (int i=0; i<row; i++)
					for (int j=0; j<col; j++) operator()(i,j) = data[j][i];
			}
		}

		/** Sequence constructor.
		 * Create integer range column vector.
		 * e.g. Matrix(3) gives [0, 1, 2]'
		 */
		Matrix(int n): row(n), col(1), alloc(false) {
			allocate();
			for (int i=0; i<n; i++) operator()(i)=i;
		}

		/** single valued (row x col) matrix.
		 * uniform matrix of value v with (row x col) is constructed.
		 * @param val T value assigned uniformly.
		 * @param row row length.
		 * @param col column length.
		 */
		Matrix(T val,int row,int col): row(row), col(col), alloc(false) {
			allocate();
			for(int i=0; i<row; i++) 
				for(int j=0; j<col; j++) m[j][i]=val;
		}

		/** zero valued (row x col) matrix.
		 * same as calling Matrix(0.0, row, col)
		 * @see Matrix(T, int, int)
		 * @param row row length.
		 * @param col column length.
		 */
		Matrix(int row,int col): row(row), col(col), alloc(false) {
			allocate();
			for (int i=0; i<col; i++)
				memset(m[i], 0, sizeof(T)*col);
		}

		/** compilation-of-vectors constructor.
		 * Vectors must be row or column (1-dimensional) vectors.
		 * Compilation direction is depends on
		 * the dimension of the FIRST vector. (i.e. v[0])
		 * CAUTION: USE WITH CARE!
		 * THIS METHOD DOES NOT CHECK IF ALL VECTORS HAVE 
		 * THE SAME DIMENSION! 
		 * @param v array of vectors.
		 * @param num number of vectors.
		 */
		Matrix(const Matrix v[], int num): alloc(false) {
			int dim = v[0].length();
			if (v[0].row == 1) {  // case of column vector
				row = num;
				col = dim;
				allocate();
				for (int i=0; i<row; i++)
					for (int j=0; j<col; j++) operator()(i,j) = v[i](j);
			} else if (v[0].col == 1) { // case of row vector
				row = dim;
				col = num;
				allocate();
				for (int j=0; j<col; j++)
					for (int i=0; i<row; i++) operator()(i,j) = v[j](i);
			} else {
				rassert(0, "elements are not 1D-vectors\n");
			}
		}

		/** copy constructor */
		Matrix(const Matrix& a, bool alloc=false): row(a.row), col(a.col), alloc(alloc) {
			if (alloc) {
				m = a.m;
			} else {
				allocate();
				for (int i=0; i<row; i++)
					for (int j=0; j<col; j++) operator()(i,j) = a(i,j);
			}
		}
		/** copy constructor (castable) */
		template<class TT> Matrix(const Matrix<TT>& a): row(a.rows()), col(a.cols()), alloc(false) {
			allocate();
			for (int i=0; i<row; i++)
				for (int j=0; j<col; j++) operator()(i,j) = (T)a(i,j);
		}

		/** default constructor.
		 * empty 0x0 matrix is created.
		 */
		Matrix(): m(NULL), row(0), col(0), alloc(false) {};
		/// destructor.
		virtual ~Matrix() {
			if (alloc) return;
			if (m) de_allocate();
		}

		/** pointer retriever.
		 * Use with extreme care.
		 */
		T** getPointer() {
			return m;
		}
		//template <class TT> Matrix& operator =(const Matrix<TT>& a) { }
		/// Matrix-Matrix assignment
		Matrix& operator =(const Matrix& b) {
			if (this == &b) return *this;
			de_allocate();
			row = b.row;
			col = b.col;
			allocate();
			for(int i=0; i<row; i++)
				for(int j=0; j<col; j++) operator()(i,j) = b(i,j);
			return *this;
		}
		/// T-Matrix assignment (uniform)
		Matrix& operator =(const T& a) {
			for(int i=0; i<row; i++)
				for(int j=0; j<col; j++) operator()(i,j) = a;
			return *this;
		}
		/** @name logical operators
		 * @{
		 */
		/** bool operation to find which element is finite (not Nan / infinity).
		 * @return bool matrix
		 */
		const Matrix finite() const {
			Matrix<T> res(row,col);
			for (int i=0; i<row; i++)
				for (int j=0; j<col; j++) res(i,j) = -isfinite(operator()(i,j));
			return res;
		}
		/// comparison operator
		friend const Matrix<bool> operator ==(const Matrix& a, const T& b) {
			Matrix<bool> res(a.row, a.col);
			for (int i=0; i<a.row; i++)
				for (int j=0; j<a.col; j++) res(i,j) = a(i,j)==b;
			return res;
		}
		/// comparison operator
		friend const Matrix<bool> operator !=(const Matrix& a, const T& b) {
			Matrix<bool> res(a.row, a.col);
			for (int i=0; i<a.row; i++)
				for (int j=0; j<a.col; j++) res(i,j) = a(i,j)!=b;
			return res;
		}
		/// comparison operator
		friend const Matrix<bool> operator <(const Matrix& a, const T& b) {
			Matrix<bool> res(a.row, a.col);
			for (int i=0; i<a.row; i++)
				for (int j=0; j<a.col; j++) res(i,j) = a(i,j)<b;
			return res;
		}
		/// comparison operator
		friend const Matrix<bool> operator >(const Matrix& a, const T& b) {
			Matrix<bool> res(a.row, a.col);
			for (int i=0; i<a.row; i++)
				for (int j=0; j<a.col; j++) res(i,j) = a(i,j)>b;
			return res;
		}
		/// comparison operator
		/*
		friend const Matrix<bool> operator ==(const Matrix& a, const Matrix& b) {
			Matrix<bool> res(a.row, a.col);
			for (int i=0; i<a.row; i++)
				for (int j=0; j<a.col; j++) res(i,j) = a(i,j)==b(i,j);
			return res;
		}*/
		friend bool operator ==(const Matrix& a, const Matrix& b) {
			for (int i=0; i<a.row; i++)
				for (int j=0; j<a.col; j++)
					if (a(i,j)!=b(i,j)) return false;
			return true;
		}
		/// comparison operator
		friend const Matrix<bool> operator !=(const Matrix& a, const Matrix& b) {
			Matrix<bool> res(a.row, a.col);
			for (int i=0; i<a.row; i++)
				for (int j=0; j<a.col; j++) res(i,j) = a(i,j)!=b(i,j);
			return res;
		}
		/// comparison operator
		friend const Matrix<bool> operator <(const Matrix& a, const Matrix& b) {
			Matrix<bool> res(a.row, a.col);
			for (int i=0; i<a.row; i++)
				for (int j=0; j<a.col; j++) res(i,j) = a(i,j)<b(i,j);
			return res;
		}
		/// comparison operator
		friend const Matrix<bool> operator >(const Matrix& a, const Matrix& b) {
			Matrix<bool> res(a.row, a.col);
			for (int i=0; i<a.row; i++)
				for (int j=0; j<a.col; j++) res(i,j) = a(i,j)>b(i,j);
			return res;
		}
		/// logical "or" operator
		friend const Matrix<bool> operator ||(const Matrix& a, const Matrix& b) {
			Matrix<bool> res(a.row, a.col);
			for (int i=0; i<a.row; i++)
				for (int j=0; j<a.col; j++) res(i,j) = a(i,j)||b(i,j);
			return res;
		}
		/// logical "and" operator
		friend const Matrix<bool> operator &&(const Matrix& a, const Matrix& b) {
			Matrix<bool> res(a.row, a.col);
			for (int i=0; i<a.row; i++)
				for (int j=0; j<a.col; j++) res(i,j) = a(i,j)&&b(i,j);
			return res;
		}
		const Matrix operator !() const {/// boolean "not" prepend ope.
			Matrix res(row,col);
			for (int i=0; i<row; i++)
				for (int j=0; j<col; j++) res(i,j) = -!operator()(i,j);
			return res;
		}
		//@}
		/** @name Arithmetic Operators
		 * Simple addition/subtraction/multiplication/division arithmetic
		 * for both Matrix-Matrix operations and Matrix-Scalar operations.
		 * @{
		 */
		/// sumation operator
		friend const Matrix operator +(const T& a, const Matrix& b) { return operator+(b,a); }
		friend const Matrix operator +(const Matrix& a, const T& b) {
			Matrix c(a.row, a.col);
			for (int i=0; i<a.row; i++)
				for (int j=0; j<a.col; j++) c(i,j) = a(i,j) + b;
			return c;
		}
		/// sumation operator
		friend const Matrix operator +(const Matrix& a, const Matrix& b) {
			Matrix c(a.row, a.col);
			rassert(a.row == b.row && a.col == b.col, "operator+: dimension mismatch.\n");
			for (int i=0; i<a.row; i++)
				for (int j=0; j<a.col; j++) c(i,j) = a(i,j) + b(i,j);
			return c;
		}
		/// subtraction operator
		friend const Matrix operator -(const T& a, const Matrix& b) { return operator+(b,-a); }
		friend const Matrix operator -(const Matrix& a, const T& b) {
			Matrix c(a.row, a.col);
			for (int i=0; i<a.row; i++)
				for (int j=0; j<a.col; j++) c(i,j) = a(i,j) - b;
			return c;
		}
		/// subtraction operator
		friend const Matrix operator -(const Matrix& a,const Matrix& b) {
			Matrix c(a.row, a.col);
			rassert(a.row == b.row && a.col == b.col, "operator-: dimension mismatch.\n");
			for (int i=0; i<a.row; i++)
				for (int j=0; j<a.col; j++) c(i,j) = a(i,j) - b(i,j);
			return c;
		}
		/// matrix multplication operator
		friend const Matrix operator *(const Matrix& a, const T& coef) { return operator*(coef,a); }
		friend const Matrix operator *(const T& coef, const Matrix& a) {
			Matrix c(a.row, a.col);
			for (int i=0; i<a.row; i++)
				for (int j=0; j<a.col; j++) c(i,j) = a(i,j) * coef;
			return c;
		}
		/// matrix multplication operator
		friend const Matrix operator *(const Matrix& a, const Matrix& b) {
			Matrix c(a.row, b.col);
			rassert(a.col == b.row, "operator*: dimension mismatch.\n");
			for(int i=0; i<a.row; i++) {
				for(int j=0; j<b.col; j++) {
					T sum = 0;
					for (int k=0; k<a.col; k++) sum += a(i,k) * b(k,j);
					c(i,j) = sum;
				}
			}
			return c;
		}
		/// elementwise division operator
		friend const Matrix operator /(const Matrix& a, const T& b) {
			Matrix c(a.row, a.col);
			for (int i=0; i<a.row; i++)
				for (int j=0; j<a.col; j++) c(i,j) = a(i,j) / b;
			return c;
		}
		/// elementwise division operator
		friend const Matrix operator /(const Matrix& a, const Matrix& b) {
			Matrix c(a.row, a.col);
			rassert(a.row == b.row && a.col == b.col, "operator/: dimension mismatch.\n");
			for (int i=0; i<a.row; i++) {
				for (int j=0; j<a.col; j++) {
					rassert(0!=b(i,j), "operator/: division by zero.\n");
					c(i,j) = a(i,j) / b(i,j);
				}
			}
			return c;
		}

		const Matrix operator -() const {/// minus prepend ope.
			Matrix res(row,col);
			for (int i=0; i<row; i++)
				for (int j=0; j<col; j++) res(i,j) = -operator()(i,j);
			return res;
		}
		Matrix operator +=(const T& c) {/// method in group1
			for (int i=0; i<row; i++)
				for (int j=0; j<col; j++) operator()(i,j) += c;
			return *this;
		}
		Matrix operator +=(const Matrix& c) {
			rassert(row == c.row && col == c.col, "operator+=: dimension mismatch.\n");
			for (int i=0; i<row; i++)
				for (int j=0; j<col; j++) operator()(i,j) += c(i,j);
			return *this;
		}
		Matrix operator -=(const T& c) {
			for (int i=0; i<row; i++)
				for (int j=0; j<col; j++) operator()(i,j) -= c;
			return *this;
		}
		Matrix operator -=(const Matrix& c) {
			rassert(row==c.row && col==c.col ,"operator-=: dimension mismatch.\n");
			for (int i=0; i<row; i++)
				for (int j=0; j<col; j++) operator()(i,j) -= c(i,j);
			return *this;
		}
		Matrix operator *=(const T& c) {
			for (int i=0; i<row; i++)
				for (int j=0; j<col; j++) operator()(i,j) *= c;
			return *this;
		}
		Matrix operator *=(const Matrix& c) {
			Matrix res(row, c.col);
			rassert(col==c.row ,"operator*=: dimension mismatch.\n");
			for(int i=0; i<row; i++) {
				for(int j=0; j<c.col; j++) {
					T sum = 0;
					for (int k=0; k<col; k++) sum += operator()(i,k) * c(k,j);
					res(i,j) = sum;
				}
			}
			operator=(res);
			return *this;
		}
		Matrix operator /=(const T& c) {
			for (int i=0; i<row; i++)
				for (int j=0; j<col; j++) operator()(i,j) /= c;
			return *this;
		}
		Matrix operator /=(const Matrix& c) {
			for (int i=0; i<row; i++)
				for (int j=0; j<col; j++) operator()(i,j) /= c(i,j);
			return *this;
		}
		//@}

		/** @name Accessors
		 * @{
		 */
		/** 2D accessor.
		 * Index starts from 0 not 1.
		 * @note no boundary checking applied here.  Simply don't overrun.
		 * @param i row index
		 * @param j col index
		 * @return value at i-th row and j-th column.
		 */
		virtual const T operator()(int i,int j) const {assert(bchk(i,j)); return m[j][i];}
		/// 2D accessor.
		virtual T& operator()(int i,int j) {assert(bchk(i,j)); return m[j][i];}
		/// 1D accessor.
		virtual const T operator()(int i) const {assert(bchk(i,0)); return m[0][i];}
		/// 1D accessor.
		virtual T& operator()(int i) {assert(bchk(i,0)); return m[0][i];}
		//@}
		/** @name Attribute accessors
		 * @{
		 */
		/// "number of rows" accessor
		const int rows() const { return row; }
		/// "number of cols" accessor
		const int cols() const { return col; }
		/** size / dimension accessor.
		 * @param dir returns row if zero, column otherwise.
		 * @return ( dir==0 ? row : col )
		 */
		int size(int dir) const { return dir==0 ? row : col; }
		/** size / dimension accessor.
		 * returns 2D dimension with a vector.
		 * @return column vector of [row; col].
		 */
		const Matrix size() const {
			Matrix s(2,1);
			s(0,0) = row;
			s(1,0) = col;
			return s;
		}
		/** length accessor.
		 * @return longer dimension
		 */
		int length() const { return row > col ? row : col; }
		/** width accessor.
		 * @return shorter dimension
		 */
		int width() const  { return row < col ? row : col; }
		/** diagonal value accessor.
		 * @return column vector with diagonal values
		 */
		const Matrix diag() const {
			rassert(row == col ,"diag(): dimension mismatch.\n");
			Matrix res(row,1);
			for (int i=0; i<row; i++) res(i,0) = operator()(i,i);
			return res;
		}
		//@}
		/** @name Matrix Composition/Decomposition
		 * @{
		 */
		/** Right-Left matrix concatenation.
		 * Pseudo Code: @code
		 * cout << A;
		 *  #> [ 1 2 ]
		 * cout << B;
		 *  #> [ 3 4 ]
		 * cout << (A | B);
		 *  #> [ 1 2 3 4 ]
		 * @endcode
		 */
		const Matrix operator |(const Matrix& b) {  // left-right matrix concatenation.
			rassert(row == b.row, "operator|: unmatched #row.\n");
			Matrix c(row, col + b.col);
			for (int j=0; j<col; j++)
				for (int i=0; i<row; i++) c(i,j) = operator()(i,j);
			for (int j=0; j<b.col; j++)
				for (int i=0; i<row; i++) c(i,j+col) = b(i,j);
			return c;
		}
		/** Top-Bottom matrix concatenation.
		 * Pseudo Code: @code
		 * cout << A;
		 *  #> [ 1 2 ]
		 * cout << B;
		 *  #> [ 3 4 ]
		 * cout << (A ^ B);
		 *  #> [ 1 2 ]
		 *  #> [ 3 4 ]
		 * @endcode
		 */
		const Matrix operator ^(const Matrix& b) {  // top-bottom matrix concatenation.
			rassert(col == b.col, "operator^: unmatched #col.\n");
			Matrix c(row + b.row, col);
			for (int i=0; i<row; i++)
				for (int j=0; j<col; j++) c(i,j) = operator()(i,j);
			for (int i=0; i<b.row; i++)
				for (int j=0; j<col; j++) c(i+row,j) = b(i,j);
			return c;
		}
		/** Matrix decompose.
		 * Decompile the matrix into array of vectors.
		 * Pseudo Code: @code
		 * cout << X;
		 *  #> [1 2]
		 *  #> [3 4]
		 * Matrix *x = new Matrix[2];
		 * X.decompose(x);
		 * cout << x[0];
		 *  #> [1 2]
		 * cout << x[1];
		 *  #> [3 4]
		 * @endcode
		 * @param v matrix array for storage.
		 *           CAUTION!: Array size must be large enough.
		 * @param row returns row vectors if true, column(default) vectors otherwise.
		 */
		const void decompose(Matrix *v, bool row=false) const {//top-bottom matrix de-compilation.
			if (row) {
				for (int i=0; i<row; i++)
					v[i] = vector_row(i);
			} else {
				for (int i=0; i<col; i++)
					v[i] = vector_col(i);
			}
		}
		/** partial matrix operator.
		 * Negative indexing(index from the end) is supported for
		 *  all the partial operations.
		 * @param i0 starting row index
		 * @param i1 ending row index
		 * @param j0 starting col index
		 * @param j1 ending row index
		 * @return partial matrix of the range specified by params.
		 *
		 * Pseudo sample code:
		 * @code
		 * cout << A;
		 *  #> [ 1 2 3 ]
		 *  #> [ 4 5 6 ]
		 * cout << A(0,0, 1,-1);
		 *  #> [ 2 3 ]
		 * @endcode
		 */
		const Matrix partial(int i0, int i1, int j0, int j1) const {
			if (i0<0) i0+=row;
			if (i1<0) i1+=row;
			if (j0<0) j0+=col;
			if (j1<0) j1+=col;
			// must comply; 0<=i0<=i1<row, 0<=j0<=j1<col (%d,%d,%d,%d)
			rassert(0<=i0 && i0<=i1 && i1<row && 0<=j0 && j0<=j1 && j1<col ,"operator(): out-of-range.\n");
			int row1 = i1 - i0 + 1;
			int col1 = j1 - j0 + 1;
			Matrix res(row1,col1);
			for (int i=i0; i<i1+1; i++)
				for (int j=j0; j<j1+1; j++) res(i-i0,j-j0) = operator()(i,j);
			return res;
		}
		/** partial matrix operator.
		 * @param b a bool vector to specify which row(s) or col(s) to be chosen.
		 *          Must be either row or column vector.
		 * @return partial matrix logically-indexed with param.
		 */
		const Matrix partial(const Matrix<bool>& b) const {
			if (b.size(0)==row && b.size(1)==1) { //column vector
				int size = 0;
				for (int i=0; i<row; i++) { if (b(i)!=0) size++; }
				Matrix res(size,col);
				int k = 0;
				for (int i=0; i<row; i++) { 
					if (b(i)!=0) {
						for (int j=0; j<col; j++) { 
							res(k,j) = operator()(i,j);
						}
						k++;
					}
				}
				return res;
			} else if (b.size(0)==1 && b.size(1)==col) { //row vector
				int size = 0;
				for (int i=0; i<col; i++) { if (b(i)!=0) size++; }
				Matrix res(row,size);
				int k = 0;
				for (int i=0; i<col; i++) { 
					if (b(i)!=0) {
						for (int j=0; j<row; j++) { 
							res(j,k) = operator()(j,i);
						}
						k++;
					}
				}
				return res;
			} else {
				rassert(0, "partial(): argument must be column or row vector.\n");
			}
			return Matrix(); // dummy
		}
		/** partial matrix operator.
		 * @param b a bool vector to specify which row(s) to be chosen.
		 * @param c a bool vector to specify which col(s) to be chosen.
		 * @return partial matrix logically-indexed with params.
		 */
		const Matrix partial(const Matrix<bool>& b, const Matrix<bool>& c) const {
			int rowsize = 0, colsize = 0;
			for (int i=0; i<row; i++) { if (b(i)!=0) rowsize++; }
			for (int i=0; i<col; i++) { if (c(i)!=0) colsize++; }
			Matrix res(rowsize,colsize);
			int k = 0;
			for (int i=0; i<row; i++) { 
				if (b(i)) {
					int l=0;
					for (int j=0; j<col; j++) { 
						if (c(j)) res(k,l++) = operator()(i,j);
					}
					k++;
				}
			}
			return res;
		}
		/** matrix permutation.
		 * @param index A vector of row indexes.
		 *          Can be row or column vector.
		 *          Must be a vector.
		 * @return permuted matrix.
		 *
		 * It works like... (pseudo-code):
		 * @code
		 * cout << A;
		 *  #> [10 20 30 40]
		 * cout << B;
		 *  #> [1 0 3 2]
		 * cout << A.permute(B);
		 *  #> [20 10 40 30]
		 * @endcode
		 */
		const Matrix permute(const Matrix& index) const {
			if ( index.size(0) > index.size(1) ) {
				Matrix res(index.size(0), col);
				for (int i=0; i<index.size(0); i++)
					for (int j=0; j<col; j++)  res(i,j) = operator()((int)index(i),j);
				return res;
			} else {
				Matrix res(row, index.size(1));
				for (int i=0; i<row; i++)
					for (int j=0; j<index.size(1); j++)  res(i,j) = operator()(i,(int)index(j));
				return res;
			}
		}
		/** matrix permutation.
		 * @param srow A vector of row indexes.
		 *              Can be row or column vector.
		 *              Must be a vector.
		 *              Can be a null vector. i.e. Matrix()
		 * @param scol a vector of col indexes.
		 * @return permuted matrix.
		 * @see permute(const Matrix&, const Matrix&) const
		 */
		const Matrix permute(const Matrix& srow, const Matrix& scol) const {
			int nrow = srow.size(0) > srow.size(1) ? srow.size(0) : srow.size(1);
			int ncol = scol.size(0) > scol.size(1) ? scol.size(0) : scol.size(1);
			if (nrow == 0) {
				Matrix res(row, ncol);
				for (int i=0; i<row; i++)
					for (int j=0; j<ncol; j++) res(i,j) = operator()(i, (int)scol(j));
				return res;
			} else if (ncol == 0) {
				Matrix res(nrow, col);
				for (int i=0; i<nrow; i++)
					for (int j=0; j<col; j++) res(i,j) = operator()((int)srow(i), j);
				return res;
			} else {
				Matrix res(nrow, ncol);
				for (int i=0; i<nrow; i++)
					for (int j=0; j<ncol; j++) res(i,j) = operator()((int)srow(i), (int)scol(j));
				return res;
			}
			return Matrix();
		}

		/** row vector extractor.
		 * @param i index specifier
		 * @return vector of i-th row
		 *
		 * Pseudo sample code:
		 * @code
		 * cout << A;
		 *  #> [1 2]
		 *  #> [3 4]
		 * cout << A.vector_row(1);
		 *  #> [3 4]
		 * @endcode
		 * @see vector_col(int) const
		 */
		const Matrix vector_row(int i) const {
			rassert(0<=i && i<row ,"vector_row(): index out of range.\n");
			Matrix vec(1, col);
			for (int j=0; j<col; j++) vec(0, j)=m[j][i];
			return vec;
		}
		/** column vector extractor.
		 * @param j index specifier
		 * @return vector of j-th column
		 * @see vector_row(int) const
		 */
		const Matrix vector_col(int j) const {
			rassert(0<=j && j<col ,"vector_col(): index out of range.\n");
			Matrix vec(m[j], row, 1);
			return vec;
		}
		//@}
		/** @name Basic linear algebra operations.
		 * Basic stuff you can find from text books.
		 * @{
		 */
		/** inner product
		 */
		friend T dot(const Matrix& a, const Matrix& b) {
			rassert(a.row == 1 || a.col == 1 ,"A must be a vector.\n");
			rassert(b.row == 1 || b.col == 1 ,"B must be a vector.\n");
			rassert(a.row * a.col == b.row * b.col ,"dimension mismatch.\n");
			T sum = 0;
			for (int i=0; i<a.row*a.col; i++) sum += a(i) * b(i);
			return sum;
		}

		/** Transpose operation which simply flips row and column.
		 * @return transposed matrix
		 */
		const Matrix t() const {
			Matrix v = Matrix(col,row);
			for (int i=0; i<row; i++)
				for (int j=0; j<col; j++) v(j,i) = operator()(i,j);
			return v;
		}
		T det(void) const{///determinant.
			rassert(row == col, "Matrix must be square.\n");
			Matrix<T> D = triangle();

			T determinant = 1;
			for (int i=0; i<col; i++) {
				determinant *= D(i,i);
			}
			return determinant;
		}
		/** Triangulation.
		 * achieve upper-triangle matrix
		 *
		 * @param[in] min: Minimum eigenvalue limit.
		 *       Number of eigen values/vectors are limited by this threshold.
		 * @return upper-triangle matrix
		 */
		const Matrix triangle(double min=0) const {
			Matrix<T> *V = new Matrix[row];
			for (int i=0; i<col; i++) {
				V[i] = Matrix(m[i], row, 1);
			}
			Matrix<T> D(row,row);

			for (int i=0; i<col; i++) {
				for (int j=0; j<i; j++) {
					T tmp = dot(V[i],V[j]);
					V[i] = V[i] - V[j] * tmp;
					D(i,j) = tmp;
				}
				T tmp = V[i].norm();
				if (tmp > min) V[i] /= tmp;
				D(i,i) = tmp;
			}
			delete[] V;
			return D;
		}
		/** Eigen value decomposition / diagonalization.
		 * if "V = A.eig(L)" is called, then  A V = V L  or  A = V L V' .
		 * Where I = V' V = V V', L is diagonal eigen value matrix.
		 * The eigen values are sorted from highest values.
		 *
		 * @param[in]  L: Square matrix or 1D-vector.
		 *       If the dimension of the matrix/vector is n, then n-th largest
		 *       eigen values/vectors will be returned.
		 * @param[out] L: Either diagonal matrix or 1D-vector of eigen-values.
		 * @param[in] min: Minimum eigenvalue limit.
		 *       Number of eigen values/vectors are limited by this threshold.
		 * @return P matrix of column eigen vectors
		 */
		const Matrix eig(Matrix& L, double min=0) const {
			rassert(row == col ,"eig(): dimension mismatch.\n");
			// check whether it is symmetric
			for (int i=0; i<row; i++) {
				for (int j=0; j<col; j++) {
					rassert(fabs(operator()(i,j)-operator()(j,i))<_INFINITESIMAL_SIZE_, "eig(): not symmetric.\n");
				}
			}
			Matrix V(*this);

			int n = row;
			double *d = new double[n];
			double *e = new double[n];

			for (int j=0; j<n; j++) { d[j] = V(n-1,j); }
			for (int i=n-1; i>0; i--) {
				// Scale to avoid under/overflow.
				double scale = 0.0;
				double h = 0.0;
				for (int k=0; k<i; k++) {
					scale = scale + fabs(d[k]);
				}
				if (scale == 0.0) {
					e[i] = d[i-1];
					for (int j = 0; j < i; j++) {
						d[j] = V(i-1,j);
						V(i,j) = 0.0;
						V(j,i) = 0.0;
					}
				} else {
					// Generate Householder vector.
					for (int k = 0; k < i; k++) {
						d[k] /= scale;
						h += d[k] * d[k];
					}
					double f = d[i-1];
					double g = ::sqrt(h);
					if (f > 0) {
						g = -g;
					}
					e[i] = scale * g;
					h = h - f * g;
					d[i-1] = f - g;
					for (int j = 0; j < i; j++) {
						e[j] = 0.0;
					}

					// Apply similarity transformation to remaining columns.
					for (int j = 0; j < i; j++) {
						f = d[j];
						V(j,i) = f;
						g = e[j] + V(j,j) * f;
						for (int k = j+1; k <= i-1; k++) {
							g += V(k,j) * d[k];
							e[k] += V(k,j) * f;
						}
						e[j] = g;
					}
					f = 0.0;
					for (int j = 0; j < i; j++) {
						e[j] /= h;
						f += e[j] * d[j];
					}
					double hh = f / (h + h);
					for (int j = 0; j < i; j++) {
						e[j] -= hh * d[j];
					}
					for (int j = 0; j < i; j++) {
						f = d[j];
						g = e[j];
						for (int k = j; k <= i-1; k++) {
							V(k,j) -= (f * e[k] + g * d[k]);
						}
						d[j] = V(i-1,j);
						V(i,j) = 0.0;
					}
				}
				d[i] = h;
			}

			// Accumulate transformations.
			for (int i = 0; i < n-1; i++) {
				V(n-1,i) = V(i,i);
				V(i,i) = 1.0;
				double h = d[i+1];
				if (h != 0.0) {
					for (int k = 0; k <= i; k++) {
						d[k] = V(k,i+1) / h;
					}
					for (int j = 0; j <= i; j++) {
						double g = 0.0;
						for (int k = 0; k <= i; k++) {
							g += V(k,i+1) * V(k,j);
						}
						for (int k = 0; k <= i; k++) {
							V(k,j) -= g * d[k];
						}
					}
				}
				for (int k = 0; k <= i; k++) {
					V(k,i+1) = 0.0;
				}
			}
			for (int j = 0; j < n; j++) {
				d[j] = V(n-1,j);
				V(n-1,j) = 0.0;
			}
			V(n-1,n-1) = 1.0;
			e[0] = 0.0;
			// Symmetric tridiagonal QL algorithm.
			for (int i = 1; i < n; i++) { e[i-1] = e[i]; }
			e[n-1] = 0.0;

			double f = 0.0;
			double tst1 = 0.0;
			double eps = pow(2.0,-52.0);
			for (int l = 0; l < n; l++) {
				// Find small subdiagonal element
				{
					double tmp = fabs(d[l])+fabs(e[l]);
					tst1 = ( tst1 > tmp ? tst1 : tmp );
				}
				int m = l;
				while (m < n) {
					if (fabs(e[m]) <= eps*tst1) {
						break;
					}
					m++;
				}

				// If m == l, d[l] is an eigenvalue,
				// otherwise, iterate.
				if (m > l) {
					int iter = 0;
					do {
						iter = iter + 1;  // (Could check iteration count here.)

						// Compute implicit shift
						double g = d[l];
						double p = (d[l+1] - g) / (2.0 * e[l]);
						double r = ::sqrt(p*p + 1.0);
						if (p < 0) {
							r = -r;
						}
						d[l] = e[l] / (p + r);
						d[l+1] = e[l] * (p + r);
						double dl1 = d[l+1];
						double h = g - d[l];
						for (int i = l+2; i < n; i++) {
							d[i] -= h;
						}
						f = f + h;

						// Implicit QL transformation.

						p = d[m];
						double c = 1.0;
						double c2 = c;
						double c3 = c;
						double el1 = e[l+1];
						double s = 0.0;
						double s2 = 0.0;
						for (int i = m-1; i >= l; i--) {
							c3 = c2;
							c2 = c;
							s2 = s;
							g = c * e[i];
							h = c * p;
							r = ::sqrt(p*p + e[i]*e[i]);
							e[i+1] = s * r;
							s = e[i] / r;
							c = p / r;
							p = c * d[i] - s * g;
							d[i+1] = h + s * (c * g + s * d[i]);

							// Accumulate transformation.

							for (int k = 0; k < n; k++) {
								h = V(k,i+1);
								V(k,i+1) = s * V(k,i) + c * h;
								V(k,i) = c * V(k,i) - s * h;
							}
						}
						p = -s * s2 * c3 * el1 * e[l] / dl1;
						e[l] = s * p;
						d[l] = c * p;

						// Check for convergence.

					} while (fabs(e[l]) > eps*tst1);
				}
				d[l] = d[l] + f;
				e[l] = 0.0;
			}
			// Sort eigenvalues and corresponding vectors.
			for (int i = 0; i < n-1; i++) {
				int k = i;
				double p = d[i];
				for (int j = i+1; j < n; j++) {
					if (::fabs(d[j]) > ::fabs(p)) {
						k = j;
						p = d[j];
					}
				}
				if (k != i) {
					d[k] = d[i];
					d[i] = p;
					for (int j = 0; j < n; j++) {
						p = V(j,i);
						V(j,i) = V(j,k);
						V(j,k) = p;
					}
				}
			}
			// construction of Lambda(eigenvalue-diagonal) matrix
			int Lw=L.size(0);
			int Lc=L.size(1);
			if (min>0) {
				for (int i=1; i<n; i++) {
					if (::fabs(d[i])<min) {
						n = i;
						break;
					}
				}
			}
			if ( Lw==1 || Lc==1 ) {
				if      (Lw > n) { L=Matrix(n,1); Lw=n; }
				else if (Lc > n) { L=Matrix(1,n); Lc=n; }
				for (int i=0; i<Lw*Lc && i<n; i++)  L(i) = d[i];
			} else if ( Lw == Lc ) {
				if ( Lw + Lc == 0 || Lw > n ) { L=Matrix(n,n); Lw=Lc=n; }
				for (int i=0; i<Lw && i<n; i++)  L(i,i) = d[i];
			} else {
				rassert(0, "eig(L): L must be a vector or square matrix.\n");
			}
			delete[] d;
			delete[] e;
			return V.partial(0,-1, 0,L.length()-1);
		}
		/** apply SVD (Singular Value Decomposition).
		 * Singular Value Decomposition of the matrix A(=U S V').
		 * The decomposition is accomplished  by \n
		 *  (1) get U and S by diagonalizing A*A', and \n
		 *  (2) get V with A, S and V as follows. \n
		 * The matrix must be skinny(i.e. \#row >= \#col) for this
		 *  method to run. (Transpose it as needed.)
		 * @code
		 * (A * A') * U = A * V * S^2
		 * V = A' * U * S^-1
		 * @endcode
		 * @param[out] U: left-singular matrix (m x m).
		 * @param[out] V: right-singular matrix (n x n).
		 * @return singular value matrix S (m x n).
		 */
		const Matrix SVD(Matrix& U, Matrix& V) const {
			rassert(row>=col, "SVD(): matrix must be skinny.\n");
			Matrix M(m, row,col);
			Matrix L(row,1);
			U = (M * M.t()).eig(L);
			Matrix S(0.0, row,col);
			Matrix Si(0.0, row,col);
			for (int i=0; i<col; i++) {
				S(i,i) = ::sqrt(L(i));
				if ( L(i) > _INFINITESIMAL_SIZE_ )  Si(i,i) = 1.0/S(i,i);
			}
			V = M.t() * U * Si;
			return S;
		}
		/** Matrix inversion.
		 *
		 * i.e. It satisfies; I = A.inv() * A = A * A.inv()
		 *
		 * eig() method is used for the inversion.
		 * @return inverted matrix
		 * @see eig(Matrix&, double) const
		 */
		const Matrix inv() const {
			rassert(row == col ,"inv(): dimension mismatch.\n");
			for (int i=0; i<row; i++) { // symmetricity check
				for (int j=0; j<col; j++) {
					if (fabs(operator()(i,j)-operator()(j,i))> _INFINITESIMAL_SIZE_) {
						Matrix A(*this);
						return ( A.t() * A ).inv() * A.t();
					}
				}
			}
			Matrix L(row,row);
			Matrix V = eig(L);
			Matrix Linv(row,row);
			for (int i=0; i<row; i++) Linv(i,i) = 1.0/L(i,i);
			return ( V * Linv * V.t() );
		}
		/** covariance matrix calculator.
		 * i.e. A.cov() = A' A
		 *
		 * Where A' = [ x_1 x_2 ... x_n ]
		 * @note no column should include extremely large or small values.
		 *       Otherwise result in incorrect values.
		 * @return covariance matrix
		 */
		const Matrix cov() const {
			Matrix<double> u = mean();
			Matrix<double> S(col,col);
			double *x = new double[col];
			rassert(x, "cov(): running out memory for temporary allocation..");
			for (int i=0; i<row; i++) {
				for (int j=0; j<col; j++) x[j] = operator()(i,j) - u(j);
				for (int j=0; j<col; j++) // top half
					for (int k=j; k<col; k++)
						S(j,k) += x[j] * x[k];
			}
			for (int j=0; j<col; j++) // bottom half reflecting from top part
				for (int k=0; k<j; k++)  S(j,k) = S(k,j) /= (double)(row-1);
			delete[] x;
			return S;
		}
		//@}
		/** @name vector-wise / element-wise operations
		 * special methods to treat matrix as a collection of 
		 * vectors or even simply indivisual elements.
		 * @{
		 */
		/// element-wise *func.
		const Matrix<double> apply(double (*func)(T)) const {
			Matrix<double> res(row,col);
			for (int i=0; i<row; i++) {
				for (int j=0; j<col; j++) {
					res(i,j) = func(operator()(i,j));
				}
			}
			return res;
		}
		/// element-wise square root.
		const Matrix sqrt() const { return apply(::sqrt); }
		/// element-wise absolute values.
		const Matrix abs() const { return apply(::fabs); }
		/** row-direction(vertical) summation.
		 * @return row vector of vertical sum.
		 */
		const Matrix sum() const {
			Matrix res(1,col);
			for (int i=0; i<row; i++) {
				for (int j=0; j<col; j++) {
					res(j) += operator()(i,j);
				}
			}
			return res;
		}
		/** row-direction(vertical) average.
		 * @return row vector of vertical average.
		 */
		const Matrix mean() const {
			Matrix<double> res = sum();
			for (int j=0; j<col; j++) {
				res(j) /= (double)row;
			}
			return res;
		}
		/** row-direction(vertical) variance.
		 * @return row vector of vertical variance.
		 */
		const Matrix var() const {
			Matrix<double> m = mean();
			Matrix<double> res(1,col);
			if (row<=1) return res;
			for (int i=0; i<row; i++) {
				for (int j=0; j<col; j++) {
					double dev = operator()(i,j) - m(j);
					res(j) += dev * dev;
				}
			}
			for (int j=0; j<col; j++) {
				res(j) /= (double)(row-1);
			}
			return res;
		}
		/** row-direction(vertical) standard deviation.
		 * @return row vector of vertical std.
		 */
		const Matrix std() const { return var().sqrt(); }
		/** row-direction(vertical) max finder.
		 * @return row vector of vertical maximum values.
		 */
		const Matrix max() const {
			Matrix<T> res = vector_row(0);
			if (row==1) return res;
			for (int i=1; i<row; i++)
				for (int j=0; j<col; j++) if ( res(j) < operator()(i,j) ) res(j)=operator()(i,j);
			return res;
		}
		/** overall max finder.
		 * @return max value among all elements.
		 */
		const T max1() const {
			T res = operator()(0);
			for (int i=1; i<row*col; i++) {
				if (res < operator()(i)) res = operator()(i);
			}
			return res;
		}
		/** row-direction(vertical) min finder.
		 * @return row vector of vertical minimum values.
		 */
		const Matrix min() const {
			Matrix<T> res = vector_row(0);
			if (row==1) return res;
			for (int i=1; i<row; i++)
				for (int j=0; j<col; j++) if ( res(j) > operator()(i,j) ) res(j)=operator()(i,j);
			return res;
		}
		/** overall min finder.
		 * @return min value among all elements.
		 */
		const T min1() const {
			T res = operator()(0);
			for (int i=1; i<row*col; i++) {
				if (res > operator()(i)) res = operator()(i);
			}
			return res;
		}
		/** row-wise sorter with a key column (returns index matrix).
		 * Sort the matrix with a key values indicated by key_index.
		 * Sorting is done along row-wise(vertical) direction.
		 * @param index column index of a sorting key.
		 * @return vector with sorted index numbers.
		 *
		 * Pseudo sample code:
		 * @code
		 * cout << A;
		 *  #> [ -1.0 3.0 ]
		 *  #> [ -3.0 1.0 ]
		 *  #> [ -8.0 2.0 ]
		 * cout << A.sort(1);
		 *  #> [ 1 ]
		 *  #> [ 3 ]
		 *  #> [ 0 ]
		 * @endcode
		 */
		const Matrix sort_index(int index) const {
			class SortElem { // inner object for simplicity
				public:
					int ind;
					T val;

					static void qsort(SortElem a[], int l, int r) {
						static SortElem m;
						static int j;
						int i;
						if (r > l) {
							m = a[r];  i = l-1;  j = r;
							for (;;) {
								while (a[++i].val < m.val && i<r);
								while (a[--j].val > m.val && l<j);
								if (i>=j) break;
								swap(a[i], a[j]);
							}
							swap(a[i],a[r]);
							qsort(a,l,i-1);
							qsort(a,i+1,r);
						}
					}
					static void swap(SortElem& a, SortElem& b) {
						SortElem tmp = a;
						a = b;
						b = tmp;
					}
			};
			SortElem *se = new SortElem[row];
			for (int i=0; i<row; i++) {
				se[i].ind = i;
				se[i].val = operator()(i,index);
			}
			SortElem::qsort(se, 0, row-1);
			Matrix permute_index(row, 1);
			for (int i=0; i<row; i++) { permute_index(i) = se[i].ind; }
			delete[] se;
			return permute_index;
		}
		/** row-wise sorter with a key column (returns sorted matrix).
		 * Sort the matrix with a key values indicated by key_index.
		 * Sorting is done along row-wise(vertical) direction.
		 * @param key_index column index of a sorting key.
		 * @return sorted matrix
		 *
		 * Pseudo sample code:
		 * @code
		 * cout << A;
		 *  #> [ -1.0 3.0 ]
		 *  #> [ -3.0 1.0 ]
		 *  #> [ -8.0 2.0 ]
		 * cout << A.sort(1);
		 *  #> [ -3.0 1.0 ]
		 *  #> [ -8.0 2.0 ]
		 *  #> [ -1.0 3.0 ]
		 * @endcode
		 */
		const Matrix sort(int key_index) const {
			return permute(sort_index(key_index));
		}
		/** Row-wise swap.
		 * Swap x-th row and y-th row.
		 */
		void swap(int x, int y) {
			T tmp;
			for (int i=0; i<col; i++) {
				tmp = operator()(x,i);
				operator()(x,i) = operator()(y,i);
				operator()(y,i) = tmp;
			}
		}
		void reset(void) {///resetting all values to ZERO.
			for (int i=0; i<row; i++)
				for (int j=0; j<col; j++) operator()(i,j) = 0.0;
		}
		double norm() const {///vector norm.
			rassert(row==1 || col==1 ,"norm(): dimension mismatch.\n");
			int dim = ( col < row ? row : col );
			double res = 0;
			for (int i=0; i<dim; i++) {
				res += operator()(i) * operator()(i);
			}
			return ::sqrt(res);
		}
		T tot() const {///element-wise total sum.
			T res = 0;
			for (int i=0; i<row; i++) {
				for (int j=0; j<col; j++) {
					res += operator()(i,j);
				}
			}
			return res;
		}
		//@}
		/** @name Input/Output operations
		 * @{
		 */
		/** content dumper (human readable format).
		 * When dumped into a ascii-file.  Designed to be used to
		 * quickly view the contents.  Values can also be loaded with
		 * Matrix(FILE*) constructor with ascii-mode turned on.
		 * If the precision must be kept,
		 * use save(FILE*) / Matrix(FILE*) methods in binary mode instead.
		 * @see Matrix(istream&, bool)
		 * @see save(ostream&, bool) const
		 */
		friend ostream& operator <<(ostream& os, const Matrix& a) {
			if (a.row==0 && a.col==0) { os << "[]\n"; }
			for(int i=0; i<a.row; i++) {
				os << "[ ";
				for(int j=0; j<a.col; j++)
					os << setw(9) << setprecision(4) << a(i,j) << " ";
				os << "]\n";
			}
			return os;
		}
		/** content save.
		 * Save the matrix content(dimension and element values) into FILE pointer.
		 * Can be later loaded with constructor for reproduction.
		 * @param os ostream.
		 * @param ascii optional parameter specifying whether ascii/binary.
		 * @see Matrix(istream&, bool)
		 */
		void save(ostream& os, bool ascii=true) const {
			if (ascii) {
				for (int i=0; i<row; i++) {
					for (int j=0; j<col; j++)
						os << operator()(i,j) << " ";
					os << "\n";
				}
			} else {
				os.write((char*)&row, sizeof(int));
				os.write((char*)&col, sizeof(int));
				os.write((char*)m, sizeof(T)*row*col);
			}
		}
		/** stream loading constructor. 
		 * load from the save()'d stream.
		 * Following CSV format is supported: 
		 * @li @c The value separator is /[ ,\\t]+/.
		 * @li @c The line(row) separator is "\n".
		 *
		 * @param is istream.
		 * @param ascii optional parameter specifies whether istream is ascii/binary.
		 * @see save(ostream&, bool) const
		 */
		Matrix(istream& is, bool ascii=true) {
			if (ascii) {
				const char DELIM[] = "Ee+-1234567890.";
				char buf[_IO_BUFFER_LENGTH_];
				char *ptr = buf;
				{ // peeking the matrix dimension
					int num=0, dim=0;
					do {
						rassert(is.getline(buf, _IO_BUFFER_LENGTH_), "Matrix(): file format error.\n");
					} while (strlen(buf)<2 || *buf=='#');
					while (NULL!=(ptr = strpbrk(ptr, DELIM))) {
						ptr += strspn(ptr, DELIM);
						dim++;
					}
					is.seekg(0);
					while (is.getline(buf, _IO_BUFFER_LENGTH_).good()) {
						if (strlen(buf)>1 && *buf!='#') num++; 
					}
					is.clear();
					is.seekg(0);
					row = num;
					col = dim;
				}
				allocate();
				for (int i=0; i<row; i++) {
					ptr = buf;
					if (NULL == is.getline(ptr, _IO_BUFFER_LENGTH_)) break;
					if (strlen(ptr)<2 || *buf=='#') { i--; continue; }
					for (int k=0; k<col; k++) {
						ptr = strpbrk(ptr, DELIM);
						m[k][i] = atof(ptr);
						ptr += strspn(ptr, DELIM);
					}
				}
			} else {
				is.read((char*)&row, sizeof(int));
				is.read((char*)&col, sizeof(int));
				allocate();
				is.read((char*)m, sizeof(T)*row*col);
			}
		}
		//@}
		/** Boolean indexer Matrix.
		 * Used implicitly for boolean indexing of various leftvalue access.
		 * Should not be used directly.
		 * @note not heavily tested.
		 */
		class BoolIndexer { 
			private:
				Matrix<T> *target;
				Matrix<bool> *mask;
			public:
				int row, col;
				BoolIndexer(Matrix<T> &m, const Matrix<bool> &msk):
					target(&m), mask(new Matrix<bool>(msk)), row(m.size(0)), col(m.size(1)) {}
				virtual ~BoolIndexer() { delete mask; }
				BoolIndexer& operator =(const T& val) {
					for (int i=0; i<row; i++)
						for (int j=0; j<col; j++)
							if ((*mask)(i,j)) (*target)(i,j) = val;
					return *this;
				}
		};
		/** Leftvalued logical(boolean) indexed accessor.
		 * @param mask boolean mask matrix.
		 * @return masked leftvalue access(reference).
		 *
		 * Pseudo sample code:
		 * @code
		 * cout << A;
		 *  #> [1 3]
		 *  #> [3 2]
		 * A(A==3) = 0;  // find 3 and change them to 0.
		 * cout << A;
		 *  #> [1 0]
		 *  #> [0 2]
		 * @endcode
		 */
		BoolIndexer operator ()(const Matrix<bool>& mask) { return BoolIndexer(*this, mask); }
		/** Leftvalued logical(boolean) indexed accessor.
		 * @param rmask boolean row mask-vector.
		 * @param cmask boolean col mask-vector.
		 * @return masked leftvalue access(reference).
		 * @see operator()(const Matrix<bool>&)
		 */
		BoolIndexer operator ()(const Matrix<bool>& rmask, const Matrix<bool>& cmask) {
			Matrix<bool> mask(row, col);
			for (int i=0; i<row; i++)
				for (int j=0; j<col; j++) mask(i,j) = (rmask(i) && cmask(j));
			return BoolIndexer(*this, mask);
		}
		/** Leftvalued logical(boolean) indexed accessor.
		 * @param rmask boolean row mask-vector.
		 * @param c_idx column index.
		 * @return masked leftvalue access(reference).
		 * @see operator()(const Matrix<bool>&)
		 */
		BoolIndexer operator ()(const Matrix<bool>& rmask, const int& c_idx) {
			Matrix<bool> mask(row, col);
			for (int i=0; i<row; i++)  mask(i,c_idx) = rmask(i);
			return BoolIndexer(*this, mask);
		}
		/** Leftvalued logical(boolean) indexed accessor.
		 * @param r_idx raw index.
		 * @param cmask boolean col mask-vector.
		 * @return masked leftvalue access(reference).
		 * @see operator()(const Matrix<bool>&)
		 */
		BoolIndexer operator ()(const int& r_idx, const Matrix<bool>& cmask) {
			Matrix<bool> mask(row, col);
			for (int j=0; j<col; j++) mask(r_idx,j) = cmask(j);
			return BoolIndexer(*this, mask);
		}
#ifdef matrix2222

		/** Partial/Permuted indexer Matrix.
		 * Used for implicit permutation/partial indexing of various
		 *  leftvalue access.
		 * Should not be used directly.
		 * @note not heavily tested.
		 */
		class PermIndexer {
			private:
				Matrix *parent;
				int *rows;
				int *cols;
			public:
				int row, col;
				PermIndexer(Matrix &m, Matrix& rows, Matrix& cols) {
					row = rows.size(0)>rows.size(1) ? rows.size(0):rows.size(1);
					col = cols.size(0)>cols.size(1) ? cols.size(0):cols.size(1);
					PermIndexer::rows = new int[row];
					PermIndexer::cols = new int[col];
					for (int i=0; i<row; i++)  PermIndexer::rows[i] = (int)rows(i);
					for (int i=0; i<col; i++)  PermIndexer::cols[i] = (int)cols(i);
					parent = &m;
				}
				PermIndexer& operator =(const Matrix& b) {
					for (int i=0; i<row; i++)
						for (int j=0; j<col; j++) operator()(i,j) = b(i,j);
					return *this;
				}
				PermIndexer& operator =(const T& b) {
					for (int i=0; i<row; i++)
						for (int j=0; j<col; j++) operator()(i,j) = b;
					return *this;
				}
				virtual ~PermIndexer() { delete[] rows;  delete[] cols; }
				virtual       T& operator()(int i,int j)       {return parent->operator()(rows[i],cols[j]);}
				virtual const T  operator()(int i,int j) const {return parent->operator()(rows[i],cols[j]);}
				virtual       T& operator()(int i)       {return parent->operator()(rows[i]);}
				virtual const T  operator()(int i) const {return parent->operator()(rows[i]);}
		};
		/** Leftvalued partial accessor.
		 * @param i0 starting row index
		 * @param i1 ending row index
		 * @param j0 starting col index
		 * @param j1 ending col index
		 * @return masked leftvalue access(reference).
		 *
		 * For example (Pseudo code):
		 * @code
		 * cout << A;
		 *  #> [0 1 2 3 4]
		 * A.L(0,0, 1,2) = [-1 -2];
		 * cout << A;
		 *  #> [0 -1 -2 3 4]
		 * @endcode
		 */
		PermIndexer L(int i0, int i1, int j0, int j1) {
			if (i1<0) i1+=row;
			if (j1<0) j1+=col;
			Matrix rows(i1-i0+1,1);
			Matrix cols(j1-j0+1,1);
			for (int i=0; i<i1-i0+1; i++) rows(i) = i0 + i;
			for (int j=0; j<j1-j0+1; j++) cols(j) = j0 + j;
			return PermIndexer(*this, rows, cols);
		}
		/** Leftvalued permuted indexed accessor.
		 * @param rows vector of row permutation.
		 * @param cols vector of col permutation.
		 * @return permutated leftvalue access(reference).
		 *
		 * Pseudo sample:
		 * @code
		 * cout << A;
		 *  #> [1 2 3 4 5 6]
		 * A.L([],[0 2 4]) = [10 20 30];
		 * cout << A;
		 *  #> [10 2 20 3 30 6]
		 * @endcode
		 */
		PermIndexer L(Matrix& rows, Matrix& cols) { return PermIndexer(*this, rows, cols); }

		/** Matrix constructor from PermIndexer.
		 * This constructor is intended to be called implicitly and
		 *  should not be used directly.
		 * @param b permutation matrix
		 * @return permutated matrix
		 */
		Matrix(const PermIndexer& b) {
			row = b.row;
			col = b.col;
			allocate();
			for (int i=0; i<row; i++)
				for (int j=0; j<col; j++) operator()(i,j) = b(i,j);
		}
#endif //matrix2222
};

#ifdef matrix11111
/** Single constant value & uniform matrix.
 * It has single value for all the elements.
 */
template<class T> class Matrix1 : public Matrix<T> {
	public:
		/// uniform val of (row x col).
		Matrix1(T val, int row, int col): Matrix<T>(val,1,1) {
			Matrix<T>::row = row;
			Matrix<T>::col = col;
		}
		virtual       T& operator()(int i,int j) { 
			cerr << "Matrix1 cannot be a leftvalue.\n";  exit(1);
		}
		virtual const T  operator()(int i,int j) const { return *m; }
		virtual       T& operator()(int i) { 
			cerr << "Matrix1 cannot be a leftvalue.\n";  exit(1);
		}
		virtual const T  operator()(int i) const { return *m; }
	private:
		Matrix1(const Matrix1&);
};

/** Eye matrix with constant scale.
 * It is a constant diagonal matrix.
 */
class MatrixEye : public Matrix {
	public:
		/// Identity matrix of (dim x dim).
		MatrixEye(int dim): Matrix(1.0,1,1) {
			Matrix::row = Matrix::col = dim;
		}
		/// scale matrix of factor val.
		MatrixEye(T val, int dim): Matrix(val,1,1) {
			Matrix::row = Matrix::col = dim;
		}
		virtual       T& operator()(int i,int j) { 
			fprintf(stderr, "MatrixEye cannot be a leftvalue.\n");  exit(1);
		}
		virtual const T  operator()(int i,int j) const { return *m * (i==j); }
		virtual       T& operator()(int i) { 
			fprintf(stderr, "MatrixEye cannot be a leftvalue.\n");  exit(1);
		}
		virtual const T  operator()(int i) const { return *m; }
	private:
		MatrixEye(const Matrix&);
};
#endif //matrix1111
} //namespace matrix

#endif /* _MATRIX_ */
// vim:syntax=cpp.doxygen:
