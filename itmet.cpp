#include "Matrix.h"
#include <fstream>
#include <chrono>
double u_a(double x, double y) 
{ 	return sin(M_PI * x) * cos(M_PI * y);}
double fun(double x, double y) 
{	return 2 * M_PI * M_PI * sin(M_PI * x) * cos(M_PI * y);}
double north(double x) 
{	return u_a(x, 1.0);}
double south(double x) 
{	return u_a(x, 0.0);}
double east(double y) 
{	return u_a(1.0, y);}
double west(double y) 
{	return u_a(0.0, y);}
void jacobi(const Matrix<double>& U1, Matrix<double>& U2, const Matrix<double>& f, double h2) 
{
	int n = U2.row() - 1;
	for(int i = 1; i < n; ++i) 
		for(int j = 1; j < n; ++j) 
			U2(i,j) = (U1(i-1,j) + U1(i+1,j) + U1(i,j-1) + U1(i,j+1) + h2 * f(i,j))/4;
}
void seidel(Matrix<double>& U, const Matrix<double>& f, double h2) 
{
	int n = U.row() - 1;
	for(int i = 1; i < n; ++i) 
		for(int j = 1; j < n; ++j) 
			U(i,j) = (U(i-1,j) + U(i+1,j) + U(i,j-1) + U(i,j+1) + h2 * f(i,j))/4;
}
void sor(Matrix<double>& U, const Matrix<double>& f, double h2, double omega) 
{
	int n = U.row() - 1;
	for(int i = 1; i < n; ++i) 
	{
		for(int j = 1; j < n; ++j) 
		{
			double temp = (U(i-1,j) + U(i+1,j) + U(i,j-1) + U(i,j+1) + h2 * f(i,j))/4 - U(i,j);
			temp *= omega;
			U(i,j) += temp;
		}
	}
}



double C_norm(const Matrix<double>& A)
{
	int n = A.row();
	double norm = 0.0;
	for(int i = 0; i < n; ++i) 
	{
		for(int j = 0; j < n; ++j)
		{
			double max = abs(A(i,j));
			if(norm < max)
				norm = max;
		}
	}
	return norm;
}

double matrix_min(const Matrix<double>& A) {
	int n = A.row();
	double min = A(0,0);
	for(int i = 0; i < n; ++i) {
		for(int j = 0; j < n; ++j){
			if(min > A(i,j)) {
				min = A(i,j);
			}
		}
	}
	return min;
}
double matrix_max(const Matrix<double>& A){
	int n = A.row();
	double max = A(0,0);
	for(int i = 0; i < n; ++i) {
		for(int j = 0; j < n; ++j){
			if(max < A(i,j)) {
				max = A(i,j);
			}
		}
	}
	return max;
}
double approx(const Matrix<double>& u, const Matrix<double>& f, double h2)
{
	int n = u.row() - 1;
	double norm = 0.0;
	for(int i = 1; i < n; i++)
	{
		for(int j = 1; j < n; j++)
		{
			double max = abs((u(i-1,j) + u(i+1,j) + u(i,j-1) + u(i,j+1) - 4 * u(i,j)) / h2 + f(i,j));
			if(max > norm) 
				norm = max;
		}
	}
	return norm;
}

int main()
{
	int n = 64;
	int np = n + 1;
	int p = 3;
	double h = 1.0/n;
	double h2 = h * h;
	vector<double> x = linspace(0.0, 1.0, np);
	vector<double> y = linspace(0.0, 1.0, np);
	Matrix<double> f(n);

	Matrix<double> U1(np), U2(np);
	for(int i = 1; i < n; i++)
		for(int j = 1; j < n; j++)
			f(i,j) += fun(x[i], y[j]);
	Matrix<double> ua(np);
	for(int i = 0; i < np; i++)
		for(int j = 0; j < np; j++)
			ua(i,j) = u_a(x[i],y[j]);
	for(int i = 0; i < np; i++)
	{
		U1(i,0) = U2(i,0) = south(x[i]);
		U1(0,i) = U2(0,i) = west(y[i]);
		U1(i,n) = U2(i,n) = north(x[i]);
		U1(n,i) = U2(n,i) = east(y[i]);
	}
	double c1 = approx(U1, f, h2);
	double c4 = C_norm(U1 - ua);
	double rs, r, tmp, omega;
	char option;
	int m;
	cout << "Метод:   ";
	cin >> option;
	switch(option) 
	{
		case '1':
			m = ceil(2 * p * log(10) / (M_PI * M_PI * h2));
			rs = cos(M_PI * h);
			break;
		case '2':
			m = ceil(p * log(10) / (M_PI * M_PI * h2));
			rs = 1.0 - M_PI * M_PI * h2;
			break;
		case '3':
			m = ceil(2 * p * log(10) / (M_PI * h));
			r = cos(M_PI * h);
			omega = 2.0/(1.0 + sqrt(1.0 - r * r));
			rs = omega - 1.0; 
			break;
		default:
			return false;
	}
	cout << endl;
	cout << "Число итераций: " << m << endl;
	cout << "Спектральный радиус: " << rs << endl;
	cout << endl;
	cout << "k\t        |F-AUk|    |F-AUk|/|F_AUo|\t    |Uk-aU|   |Uk-aU|/|Uo-aU|\t |Uk-U(k-1)|   pogr.\t     rs" << endl;
	auto begin = std::chrono::steady_clock::now();
	for(int k = 1; k <= m; ++k) 
	{
		switch(option) 
		{
			case '1':
				jacobi(U1, U2, f, h2);
				break;
			case '2':
				seidel(U2,f,h2);
				break;
			case '3':
				sor(U2,f,h2,omega);
				break;
		}
		double c5 = C_norm(U2 - U1);
		double c2 = approx(U2,f,h2);
		double c3 = C_norm(U2 - ua);
		if(option == '3') 
		{
			if(k%10 == 0 || k == m) 
			{
				cout << k;
				if(k >= 100) 
				{
					printf("%14.10lf", c2);
				}
				else 
				{
					printf("%15.7lf", c2);
				}
				printf("%20.15lf %12.7lf %12.7lf %22.12lf %14.10lf %13.5lf\n", c2/c1, c3, c3 / c4, c5, rs * c5/(1.0-rs), c5 / tmp);
			}
		}
		else if(k%100 == 0 || k == m) 
		{
			cout << k;
			if(k >= 1000) 
			{
				printf("%14.10lf", c2);
			}
			else 
			{
				printf("%15.7lf", c2);
			}
			printf("%20.15lf %12.7lf %12.7lf %22.12lf %14.10lf %13.5lf\n", c2/c1, c3, c3 / c4, c5, rs * c5/(1.0-rs), c5 / tmp);
		}
		tmp = c5;
		U1 = U2;
	}
	auto end = std::chrono::steady_clock::now();
	auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
			cout <<"Time =     " << elapsed_ms.count() << " ms\n";
	ofstream out("info.txt");
	for(int i = 0; i < n + 1; i++)
		for(int j = 0; j < n + 1; j++)
        	out << x[i] <<"		" << y[j] << "	" << U1(i, j) <<endl;
	out.close();
	return 0;
}


















/*void new_seidel(Matrix<double>& U, const Matrix<double>& f, double h2) {
	int n = U.row() - 1;
	for(int i = 1; i < n; ++i) {
		for(int j = 1; j < n; ++j) {
			if((i + j)%2 == 0){
				U(i,j) = (U(i-1,j) + U(i+1,j) + U(i,j-1) + U(i,j+1) + h2 * f(i,j))/4;
			}

		}
	}
	for(int i = 1; i < n; ++i) {
		for(int j = 1; j < n; ++j) {
			if((i + j)%2 == 1){
				U(i,j) = (U(i-1,j) + U(i+1,j) + U(i,j-1) + U(i,j+1) + h2 * f(i,j))/4;
			}

		}
	}
}*/