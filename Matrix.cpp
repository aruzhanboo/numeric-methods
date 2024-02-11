#include <cmath>
#include <fstream>
#include "Matrix.h"
#include <chrono>
Matrix<double> inver(Matrix<double>);
Matrix<double> sweep(const Matrix<double> &, const Matrix<double> &);

double u_a(double x, double y)
{   return exp(2 * x) * sin(2*y);}//sin(M_PI * x) * cos(M_PI * y);   }  //x*x+y*y;}
double f(double x, double y)
{   return 0;}//2 * M_PI * M_PI * sin(M_PI * x) * cos(M_PI * y);   } //-4.0 ;} 
double west(double y)
{   return u_a(0.0, y); }
double east(double y)
{   return u_a(1.0, y); }
double south(double x)
{   return u_a(x, 0.0); }
double north(double x)
{   return u_a(x, 1.0); }
int main(void)
{
	int n = 64, nm1 = n - 1, nm2 = n - 2, np1 = n + 1;
	double h = 1.0 / n, h2 = h * h;
	vector<double> x = linspace(0.0, 1.0, np1), y = linspace(0.0, 1.0, np1);
	Matrix<double> F(nm1);
	for (int j = 0; j < nm1; j++)
	{
		for (int i = 0; i < nm1; i++)
		{
			F(j, i) += f(x[i + 1], y[j + 1]) * h2;
			if (!j)
				F(0, i) += south(x[i + 1]);
			else if (j == nm2)
				F(nm2, i) += north(x[i + 1]);
		}
		F(j, 0) += west(y[j + 1]);
		F(j, nm2) += east(y[j + 1]);
	}
	Matrix<double> c(nm1);
	for (int i = 0; i < nm1; i++)
	{
		c(i, i) = 4.0;
		if (i > 0)
			c(i - 1, i) = -1.0;
		if (i < nm1 - 1)
			c(i + 1, i) = -1.0;
	}
	auto begin = std::chrono::steady_clock::now();
	Matrix<double> U = sweep(c, F), U_a(nm1);
	auto end = std::chrono::steady_clock::now();
	auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
	double dlt = 0.0;
	int I, J;
	for (int j = 0; j < nm1; j++)
	{
		for (int i = 0; i < nm1; i++)
		{
			U_a(j, i) = u_a(x[i + 1], y[j + 1]);
			double t = fabs(U(j, i) - U_a(j, i));
			if (t > dlt)
			{
				I = i;
				J = j;
				dlt = t;
			}
		}
	}
	cout << "Погрешность = " << dlt << "	x[i]= " << x[I] << "	y[j]= " << y[J] << "	time = " << elapsed_ms.count() << " ms\n";
	ofstream ofs("m.txt");
	for (int j = 0; j < nm1; j++)
	{
		for (int i = 0; i < nm1; i++)
		{
			ofs << x[i] << ' ' << y[j] << ' ' << U(j, i) << endl;
		}
	}
	ofs.close();
	return 0;
}

Matrix<double> inver(Matrix<double> a)
{
	int n = a.row();
	Matrix<double> e(n);
	for (int i = 0; i < n; i++)
		e(i, i) = 1;
	double eps = 1.0e-3;
	for (int k = 0; k < n - 1; k++)
	{
		double max = fabs(a(k, k));
		int pos = k;
		for (int t = k + 1; t < n; t++)
		{
			if (fabs(a(t, k)) > max)
			{
				max = fabs(a(t, k));
				pos = t;
			}
		}
		if (pos != k)
		{
			for (int j = k; j < n; j++)
				swap(a(pos, j), a(k, j));
			for (int s = 0; s < n; s++)
				swap(e(pos, s), e(k, s));
		}
		for (int i = k + 1; i < n; i++)
		{
			double tmp = a(i, k) / a(k, k);
			a(i, k) = 0;
			for (int j = k + 1; j < n; j++)
				a(i, j) -= tmp * a(k, j);
			for (int s = 0; s < n; s++)
				e(i, s) -= tmp * e(k, s);
		}
	}
	Matrix<double> x(n);
	for (int k = 0; k < n; k++)
	{
		for (int i = n - 1; i >= 0; i--)
		{
			x(i, k) = e(i, k);
			double sum = 0;
			for (int j = n - 1; j > i; j--)
				sum += a(i, j) * x(j, k);
			x(i, k) = (x(i, k) - sum) / a(i, i);
		}
	}
	return x;
}

Matrix<double> sweep(const Matrix<double> &c, const Matrix<double> &F)
{
	int n = c.row();
	vector<Matrix<double>> alph;
	Matrix<double> x(n), bet(n);
	alph.push_back(inver(c));//прямой ход
	bet[0] = alph[0] * F[0];
	for (int i = 0; i < n - 1; i++)
	{
		alph.push_back(inver(c - alph[i]));
		bet[i + 1] = alph[i + 1] * (bet[i] + F[i + 1]);
	}
	x[n - 1] = bet[n - 1]; //обратный ход
	for (int i = n - 2; i >= 0; i--)
		x[i] = alph[i] * x[i + 1] + bet[i];
	return x;
}