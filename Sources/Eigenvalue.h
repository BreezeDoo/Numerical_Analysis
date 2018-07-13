//Eigenvalue.h头文件，用于计算方阵及压缩矩阵的特征值，包括幂法，反幂法，QR分解法等。
#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <assert.h>
#include "Linear_Equations.h"
#include "MyVector.h"
#include "Matrix.h"
#include <math.h>


using std::vector;
using std::string;
using std::cout;
using std::cin;
using std::istream;
using std::ostream;
using std::ifstream;
using std::ofstream;
using std::istringstream;
using std::cerr;
using std::endl;
using namespace Linear_Equations;

namespace Eigenvalue
{
	double powerMethod_algorithm(const Matrix& m_input, const vector<double>& u0, const double& epsilon, const int& iterations);
	double powerMethod(const Matrix& m_input);
	double powerMethod_for_bandMatrix_algorithm(const Matrix& m_input, const vector<double>& u0, const double& epsilon, const int& iterations, const int &s, const int &r);
	double powerMethod_for_bandMatrix(const Matrix& m_input, const int &s, const int &r);
	double inv_pow_algorithm(const Matrix& m_input, vector<double> &eigenvector, const vector<double>& u0, const double& epsilon, const int& iterations);
	double inv_pow(const Matrix& m_input, vector<double> &eigenvector);
	void QR_composition(const Matrix &A, Matrix &R, Matrix &Q);//QR分解
	Matrix hessen(const Matrix &A);//矩阵的拟上三角化
	Matrix double_displacement_QRmethod(const Matrix &m_input, const double &epsilon, const int &MAX_ITERATIONS);
}


double Eigenvalue::powerMethod_algorithm(const Matrix& m_input, const vector<double>& u0, const double& epsilon, const int& iterations)
{
	int size = m_input.cols();
	vector<double> u(u0);
	double eta;
	vector<double> y(size, 0.);
	double beta, beta_tmp;
	beta_tmp = 0;

	for (int k = 1; k < iterations; k++)
	{
		eta = sqrt(u * u);
		y = u / eta;
		u = m_input * y;
		beta = y * u;
		if ( (abs(beta - beta_tmp) / abs(beta)) < epsilon )
		{
			//cout << "The iterations right now are : " << k << endl;
			return beta;
			break;
		}
		beta_tmp = beta;
	}
	cout << "Failed: the iterations are not enough." << endl;
}


double Eigenvalue::powerMethod(const Matrix& m_input)
{
	double given_epsilon;
	given_epsilon = 1e-12;
//	cout << "Please input the required error Epsilon: " << endl;
//	cin >> given_epsilon;
	vector<double> u0(m_input.cols(), 1);

	const int MAX_INTERATIONS = 1000;
	return powerMethod_algorithm(m_input, u0, given_epsilon, MAX_INTERATIONS);
}


double Eigenvalue::powerMethod_for_bandMatrix_algorithm(const Matrix& m_input, const vector<double>& u0, const double& epsilon, const int& iterations, const int &s, const int &r)
{
	int size = m_input.cols();
	vector<double> u(u0);
	double eta;
	vector<double> y(size, 0.);
	double beta, beta_tmp;
	beta_tmp = 0;

	for (int k = 1; k < iterations; k++)
	{
		eta = sqrt(u * u);
		y = u / eta;
		u = bandMatrix_multiply_vector(m_input, s, r, y);
		beta = y * u;
		if ((abs(beta - beta_tmp) / abs(beta)) < epsilon)
		{
			//cout << "The iterations right now are : " << k << endl;
			return beta;
			break;
		}
		beta_tmp = beta;
	}
	cout << "Failed: the iterations are not enough." << endl;
}


double Eigenvalue::powerMethod_for_bandMatrix(const Matrix& m_input, const int &s, const int &r)
{
	double given_epsilon;
	given_epsilon = 1e-12;
	//	cout << "Please input the required error Epsilon: " << endl;
	//	cin >> given_epsilon;
	vector<double> u0(m_input.cols(), 1);

	const int MAX_INTERATIONS = 10000;
	return powerMethod_for_bandMatrix_algorithm(m_input, u0, given_epsilon, MAX_INTERATIONS, s, r);
}

double Eigenvalue::inv_pow_algorithm(const Matrix& m_input, vector<double> &eigenvector, const vector<double>& u0, const double& epsilon, const int& iterations)
{
	int size = m_input.cols();
	vector<double> u(u0);
	double eta;
	eigenvector.resize(size);
	double beta, beta_tmp;
	beta_tmp = 0;
	Matrix m = LU_sequential_resolve(m_input);
	for (int k = 1; k < iterations; k++)
	{
		eta = sqrt(u * u);
		eigenvector = u / eta;
		u = Substition_Doolittle(m, eigenvector);
		beta = eigenvector * u;
		if ( (fabs(1 / beta - 1 / beta_tmp) / fabs(1 / beta)) < epsilon )
		{
			//cout << "The iterations right now are : " << k << endl;
			return ( 1 / beta);//返回beta的倒数为特征值
			break;
		}
		beta_tmp = beta;
	}
	cout << "Failed: the iterations are not enough." << endl;
}

double Eigenvalue::inv_pow(const Matrix& m_input, vector<double> &eigenvector)
{
	double given_epsilon;
	given_epsilon = 1e-12;
	vector<double> u0(m_input.cols(), 1);
	const int MAX_INTERATIONS = 10000;
	return inv_pow_algorithm(m_input, eigenvector, u0, given_epsilon, MAX_INTERATIONS);
}

int sgn(const double &num)
{
	if (num <= 0)
		return -1;
	else
		return 1;
}

void Eigenvalue::QR_composition(const Matrix &A, Matrix &R, Matrix &Q)
{
	assert(A.square());
	int n = A.cols();
	Q.resize(n, n);
	for (int i = 0; i < n; i++)
	{
		Q[i][i] = 1;
	}
	R = A;
	double d;
	double c;
	double h;
	vector<double> u(n);
	vector<double> omega(n);
	vector<double> p(n);
	for (int r = 1; r < n; r++)
	{
		//检测当前步是否需要将第 r 列集中到对角线位置
		d = 0;
		for (int i = r + 1; i <= n; i++)
		{
			d += R.at(i, r) * R.at(i, r);
		}
		if (d == 0)
		{
			continue;
		}
		//下面为主算法
		else
		{
			//中间变量初始化
			d += R.at(r, r) * R.at(r, r);
			d = sqrt(d);
			c = -sgn(R.at(r, r)) * d;
			h = c * c - c * R.at(r, r);

			int i = 1;
			for (; i <= r - 1; i++)
			{
				u[i - 1] = 0;
			}
			u[i - 1] = R.at(r, r) - c;
			i++;
			for (; i <= n; i++)
			{
				u[i - 1] = R.at(i, r);
			}

			omega = Q * u;
			Q = Q - vector_multiply(omega, u / h);
			p = trans(R) * (u / h);
			R -= vector_multiply(u, p);
		}
	}
}

//拟上三角化过程，将初始矩阵化为hessen矩阵
Matrix Eigenvalue::hessen(const Matrix &A)
{
	assert(A.square());
	int n = A.cols();
	Matrix m(A);
	double d;
	double c;
	double h;
	vector<double> u(n);
	vector<double> p(n);
	vector<double> q(n);
	double t;
	vector<double> omega(n);

	for (int r = 1; r < n - 1; r++)
	{
		//检测当前步是否需要将第 r 列利用集中到对角线及对角线下一列的位置
		d = 0;
		for (int i = r + 2; i <= n; i++)
		{
			d += m.at(i, r) * m.at(i, r);
		}
		if (d == 0)
		{
			continue;
		}
		//下面为主算法
		else
		{
			//中间变量初始化
			d += m.at(r + 1, r) * m.at(r + 1, r);
			d = sqrt(d);
			c = -sgn(m.at(r + 1, r)) * d;
			h = c * c - c * m.at(r + 1, r);

			int i = 1;
			for (; i <= r; i++)
			{
				u[i - 1] = 0;
			}
			u[i - 1] = m.at(r + 1, r) - c;
			i++;
			for (; i <= n; i++)
			{
				u[i - 1] = m.at(i, r);
			}

			p = trans(m) * (u / h);
			q = m * (u / h);
			t = p * u / h;
			omega = q - u * t;
			m -= (vector_multiply(omega, u) + vector_multiply(u, p));
		}
	}
	return m;
}

void solve_quadratic_equation(const double &a, const double &b, const double &c, double &s1_re, double &s1_im, double &s2_re, double &s2_im)
{
	double delta = b * b - 4 * a * c;
	if (delta >= 0)
	{
		s1_re = (-b + sqrt(delta)) / (2 * a);
		s2_re = (-b - sqrt(delta)) / (2 * a);
		s1_im = 0;
		s2_im = 0;
	}
	else
	{
		s1_re = -b / (2 * a);
		s2_re = s1_re;
		s1_im = sqrt(-delta) / (2 * a);
		s2_im = -s1_im;
	}
}

//带双步位移的QR方法
Matrix Eigenvalue::double_displacement_QRmethod(const Matrix &m_input, const double &epsilon, const int &MAX_ITERATIONS)
{
	assert(m_input.square());
	int n = m_input.cols();
	Matrix A = hessen(m_input);
	cout << "矩阵 A 经过拟上三角化得到的矩阵 A^(n-1) 为：" << endl;
	printMatrix(A);

	int k = 0;
	int m = n;
	Matrix eigenvalue(2, n);
	while (1)
	{
		if (m == 1)
		{
			eigenvalue.at(1, 1) = A.at(1, 1);
			eigenvalue.at(2, 1) = 0;
			break;
		}
		if (fabs(A.at(m, m - 1)) <= epsilon)
		{
			eigenvalue.at(1, m) = A.at(m, m);
			eigenvalue.at(2, m) = 0;
			m--;
			continue;
		}
		if (m == 2 || fabs(A.at(m - 1, m - 2)) <= epsilon)
		{
			solve_quadratic_equation(1, -(A.at(m - 1, m - 1) + A.at(m, m)), 
				(A.at(m, m) * A.at(m - 1, m - 1) - A.at(m, m - 1) * A.at(m - 1, m)), 
				eigenvalue.at(1, m), eigenvalue.at(2, m), eigenvalue.at(1, m - 1), eigenvalue.at(2, m - 1));
			m = m - 2;
			if (m == 0)
				break;
			continue;
		}

		if (k > MAX_ITERATIONS)
		{
			cout << "The iterations are not enough.";
			break;
		}
		if (A.cols() != m)
		{
			A = submatrix(A, 1, m, 1, m);
		}

		double s = A.at(m - 1, m - 1) + A.at(m, m);
		double t = A.at(m - 1, m - 1) * A.at(m, m) - A.at(m, m - 1) * A.at(m - 1, m);
		Matrix M = A * A - s * A;
		for (int i = 0; i < m; i++)
		{
			M[i][i] += t;
		}

		Matrix B(M), C(A);

		for (int r = 1; r <= m - 1; r++)
		{
			bool flag = 0;
			for (int i = r + 1; i <= m; i++)
			{
				if (B.at(i, r) != 0)
				{
					flag = 1;
					break;
				}
			}
			if (flag == 0)
				continue;
			double d = 0;
			//中间变量初始化
			for (int i = r; i <= m; i++)
			{
				d += B.at(i, r) * B.at(i, r);
			}
			d = sqrt(d);
			double c = -sgn(B.at(r, r)) * d;
			double h = c * c - c * B.at(r, r);

			vector<double> u(m);
			u[r - 1] = B.at(r, r) - c;
			for (int i = r + 1; i <= m; i++)
			{
				u[i - 1] = B.at(i, r);
			}

			vector<double>v(trans(B) * u / h);
			B = B - vector_multiply(u, v);
			vector<double>p(trans(C) * (u / h));
			vector<double>q(C * (u / h));
			t = p * u / h;
			vector<double>omega(q - u * t);
			C = C - vector_multiply(omega, u) - vector_multiply(u, p);
		}
		A = C;
		k++;
	}
	return eigenvalue;
}