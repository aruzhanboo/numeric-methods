#include <iostream>
#include <iomanip>
#include <ctime>
#include <cmath>
#include <vector>

using namespace std;

template <class T>
class Matrix {
	vector<vector <T> > mat;
	int rows;
	int cols;
public:
	Matrix();
	Matrix(int);
	Matrix(int, int);
	Matrix(const Matrix&);
	
	template<class F>friend Matrix<F> operator- (const Matrix<F>&, const Matrix<F>&);
	template<class F>friend Matrix<F> operator* (const Matrix<F>&, const Matrix<F>&);
	template<class F>friend Matrix<F> operator* (const Matrix<F>&, const F&);
	template<class F>friend vector<F> operator* (const Matrix<F>&, const vector<F>&);
	Matrix<T> operator= (const Matrix<T>&);
	T& operator() (int, int);
	const T& operator() (int, int) const;
	vector<T>& operator[] (int);
	const vector<T>& operator[] (int) const;
	int row() const;
	int col() const;

};

//constructors

template<class T>Matrix<T>::Matrix(): rows(0), cols(0) {}
template<class T>Matrix<T>::Matrix(int m, int n): rows(m), cols(n) {
	vector<T> v(cols);
	for(int i = 0; i < rows; i++) mat.push_back(v);
}
template<class T>Matrix<T>::Matrix(int n): rows(n), cols(n) {
	vector<T> v(cols);
	for(int i = 0; i < rows; i++) mat.push_back(v);
}
template<class T>Matrix<T>::Matrix(const Matrix& M): rows(M.rows), cols(M.cols) {
	vector<T> v(cols);
	for(int i = 0; i < rows; i++) mat.push_back(v);
	for(int i = 0; i < rows; i++){
		for(int j = 0; j < cols; j++) mat[i][j] = M(i, j);
	}
}


template<class F> Matrix<F> operator- (const Matrix<F>& A, const Matrix<F>& B){
	
	Matrix<F> Res(A.rows, A.cols);
	for(int i = 0; i < Res.rows; i++){
		for(int j = 0; j < Res.cols; j++) Res(i,j) = A(i,j) - B(i,j);
	}
	return Res;
}
template<class F> Matrix<F> operator* (const Matrix<F>& A, const Matrix<F>& B){
	
	Matrix<F> Res(A.rows, B.cols);
	for(int i = 0; i < A.rows; i++){
		for(int j = 0; j < B.cols; j++){
			Res(i, j) = 0;
			for(int s = 0; s < A.cols; s++) Res(i,j) += A(i,s) * B(s,j);
		}
	}
	return Res;
}
template<class F> vector<F> operator* (const Matrix<F>& A, const vector<F>& v){
	
	vector<F> res(A.rows);
	for(int i = 0; i < A.rows; i++){
			res[i] = 0;
			for(int j = 0; j < A.cols; j++) res[i] += A(i,j) * v[j];
	}
	return res;
}
template<class F> Matrix<F> operator* (const Matrix<F>& A, const F& a){
	Matrix<F> Res(A.rows, A.cols);
	for(int i = 0; i < Res.rows; i++){
		for(int j = 0; j < Res.cols; j++) Res(i,j) = A(i,j) * a;
	}
	return Res;
}
template<class T> Matrix<T> Matrix<T>::operator= (const Matrix<T>& B){
	for(int i = 0; i < rows; i++){
		for(int j = 0; j < cols; j++) mat[i][j] = B(i,j);
	}
	return *this;
}
template<class T> T& Matrix<T>::operator() (int i, int j) {return mat[i][j];}
template<class T> const T& Matrix<T>::operator() (int i, int j) const {return mat[i][j];}
template<class T> vector<T>& Matrix<T>::operator[] (int i) {return mat[i];}
template<class T> const vector<T>& Matrix<T>::operator[] (int i) const {return mat[i];}



template<class T> int Matrix<T>::row() const {return rows;}
template<class T> int Matrix<T>::col() const {return cols;}


template<class T> vector<T> operator+ (const vector<T>& a, const vector<T>& b){
	int n = a.size();	
	vector<T> res(n);
	for(int i = 0; i < n; i++) res[i] = a[i] + b[i];
	return res;
}
template<class T> ostream& operator<< (ostream& os, const vector<T>& v){
	int n = v.size();
	for(int i = 0; i < n; i++) os << setprecision(2) << setw(10) << v[i];
	return os;
}
vector<double> linspace(double start, double end, int num){
	vector<double> linspaced;
	if (!num) {return linspaced;}
	if (num == 1){
		linspaced.push_back(start);
		return linspaced;
	}
	double delta = (end - start)/(num - 1);
	for(int i = 0; i < num - 1; i++) linspaced.push_back(start + delta * i);
	linspaced.push_back(end);
	return linspaced;
}
