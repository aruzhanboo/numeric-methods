#include <iostream>
#include <math.h>
#include <fstream>
#include "Matrix.h"
#include <vector>
#include <tuple>
#include <unistd.h>
#include <chrono>
using namespace std;
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

double dot(int m, Matrix<double> a, Matrix<double> b)
{
    double dt = 0.0;
    for (int j = 0; j < m; j++)
        for (int i = 0; i < m; i++)
            dt += a(j, i) * b(j, i);
    return dt;
}

double norm(int m, Matrix<double> a)
{   return sqrt(dot(m, a, a));  }

Matrix<double> Ax(int m, Matrix<double> x)
{
    Matrix<double> t(m);
    int j = 0;
    for (int i = 0; i < m; i++)
    {
        t(j, i) = 4 * x(j, i) - x(j + 1, i);
        if (i > 0)
            t(j, i) += -x(j, i - 1);
        if (i < m - 1)
            t(j, i) += -x(j, i + 1);
    }
    for (j = 1; j < m - 1; j++)
    {
        for (int i = 0; i < m; i++)
        {
            t(j, i) = 4 * x(j, i) - x(j - 1, i) - x(j + 1, i);
            if (i > 0)
                t(j, i) += -x(j, i - 1);
            if (i < m - 1)
                t(j, i) += -x(j, i + 1);
        }
    }
    j = m - 1;
    for (int i = 0; i < m; i++)
    {
        t(j, i) = 4 * x(j, i) - x(j - 1, i);
        if (i > 0)
            t(j, i) += -x(j, i - 1);
        if (i < m - 1)
            t(j, i) += -x(j, i + 1);
    }
    return t;
}

Matrix<double> xaddcy(int m, Matrix<double> x, double c, Matrix<double> y)
{
    Matrix<double> t(m);
    for (int j = 0; j < m; j++)
        for (int i = 0; i < m; i++)
            t(j, i) = x(j, i) + c * y(j, i);
    return t;
}

tuple<int, double, int, Matrix<double>> CG(int m, Matrix<double> b, int max_iter, double tol)
{
    Matrix<double> x(m), r(m), p(m), z(m);
    r =  p = b;                     
    double nb = norm(m, r), mu, v, r2; 
    if (nb == 0)
        nb = 1.0;
    double resid = norm(m, r) / nb;
    if (resid <= tol)
    {
        tol = resid;
        max_iter = 0;
        return make_tuple(0, tol, max_iter, x);
    }
    double r1 = dot(m, r, r);
    for (int j = 1; j < max_iter; j++)
    {
        z = Ax(m, p);
        v = r1 / dot(m, p, z);
        x = xaddcy(m, x, v, p);
        r = xaddcy(m, r,  -v, z);
        r2 = dot(m, r, r);
        resid = norm(m, r) / nb;
        if (resid <= tol)
        {
            tol = resid;
            max_iter = j;
            return make_tuple(0, tol, max_iter, x);
        }
        mu = r2 / r1;
        p = xaddcy(m, r, mu, p);
        r1 = r2;
    }
    tol = resid;
    return make_tuple(1, tol, max_iter, x);
}

int main()
{
    int n = 64, nm1 = n - 1, nm2 = n - 2, np1 = n + 1, mx_it = 100;
    double h = 1.0 / n, h2 = h * h;
    vector<double> x = linspace(0.0, 1.0, np1), y = linspace(0.0, 1.0, np1);
    Matrix<double> F(nm1), U(nm1),  U_a(nm1);                                          
    for (int j = 0; j < nm1; j++)
    {
        for (int i = 0; i < nm1; i++)
        {
            F(j, i) += f(x[i + 1], y[j + 1]) * h2;
            if (j == 0)
                F(0, i) += south(x[i + 1]);
            else if (j == nm2)
                F(nm2, i) += north(x[i + 1]);
        }
        F(j, 0) += west(y[j + 1]);
        F(j, nm2) += east(y[j + 1]);
    }
    double eps = 1.0e-3, tol; //
    int ok, it;
    auto begin = std::chrono::steady_clock::now();
    tie(ok, tol, it, U) = CG(nm1, F, mx_it, eps);
    auto end = std::chrono::steady_clock::now();
	auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
    cout << "ok = " << ok << " tol = " << tol << " Число итераций = " << it << endl;
    double dlt = 0.0, t;
    int I, J;
    for (int j = 0; j < nm1; j++)
    {
        for (int i = 0; i < nm1; i++)
        {
            U_a(j, i) = u_a(x[i + 1], y[j + 1]);
            t = fabs(U(j, i) - U_a(j, i));
            if (t > dlt)
            {
                I = i;
                J = j;
                dlt = t;
            }
        }
    }
    cout << "Погрешность = " << dlt << " x[i] = " << x[I] << " y[j] =  " << y[J] << "	time = " << elapsed_ms.count() << " ms\n";;
    ofstream ofs("cg.txt");
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
