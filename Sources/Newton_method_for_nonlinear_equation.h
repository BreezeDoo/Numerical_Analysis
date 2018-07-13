#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <assert.h>
#include "Matrix.h"
#include "MyVector.h"
#include "Linear_Equations.h"
#include <math.h>
#include <algorithm>

using std::vector;
using std::string;
using std::cout;
using std::cin;
using std::ios;
using std::setprecision;
using std::setiosflags;
using std::dec;
using std::cerr;
using std::endl;
using namespace Linear_Equations;

double fun(int n, const double &x, const double &y, const double &t, const double &u, const double &v, const double &w)
{
	if (n == 0)
		return 0.5 * cos(t) + u + v + w - x - 2.67;
	if (n == 1)
		return t + 0.5 * sin(u) + v + w - y - 1.07;
	if (n == 2)
		return 0.5 * t + u + cos(v) + w - x - 3.74;
	if (n == 3)
		return t + 0.5 * u + v + sin(w) - y - 0.79;
}

//离散Newton法求解非线性方程组
void solveNonlinearEquations(const int &dec, vector<double> &var, const double &epsilon, const int &M)
{
	double &x = var[0];
	double &y = var[1];
	double &t = var[2];
	double &u = var[3];
	double &v = var[4];
	double &w = var[5];
	double h;
	Matrix J(4, 4);
	vector<double> sol(4);
	vector<double> F(4);
	//已知t, u求解x, y
	if (dec == 0)
	{
		sol[0] = x;
		sol[1] = y;
		sol[2] = v;
		sol[3] = w;
		h = 0.1;
		for (int k = 0; k < M; k++)
		{
			for (int i = 0; i < 4; i++)
			{
				F[i] = fun(i, x, y, t, u, v, w);
			}
			for (int j = 0; j < 4; j++)
			{
				J[j][0] = (fun(j, x + h, y, t, u, v, w) - fun(j, x, y, t, u, v, w)) / h;
			}
			for (int j = 0; j < 4; j++)
			{
				J[j][1] = (fun(j, x, y + h, t, u, v, w) - fun(j, x, y, t, u, v, w)) / h;
			}
			for (int j = 0; j < 4; j++)
			{
				J[j][2] = (fun(j, x, y, t, u, v + h, w) - fun(j, x, y, t, u, v, w)) / h;
			}
			for (int j = 0; j < 4; j++)
			{
				J[j][3] = (fun(j, x, y, t, u, v, w + h) - fun(j, x, y, t, u, v, w)) / h;
			}

			vector<double> delta_sol = solve_pivot_Gauss(J, -F);
			sol = sol + delta_sol;
			x = sol[0];
			y = sol[1];
			v = sol[2];
			w = sol[3];
			if (abs(delta_sol) / abs(sol) <= epsilon)
				return;
			if (k >= M)
				cout << "Fail to solve the equaitons : the iterations are not enough." << endl;
			h = 1e-11;
			//h = abs(F);
		}
	}

	//已知x, y求解t, u
	if (dec == 1)
	{
		sol[0] = t;
		sol[1] = u;
		sol[2] = v;
		sol[3] = w;
		h = 0.1;
		for (int k = 0; k < M; k++)
		{
			for (int i = 0; i < 4; i++)
			{
				F[i] = fun(i, x, y, t, u, v, w);
			}
			for (int j = 0; j < 4; j++)
			{
				J[j][0] = (fun(j, x, y, t + h, u, v, w) - fun(j, x, y, t, u, v, w)) / h;
			}
			for (int j = 0; j < 4; j++)
			{
				J[j][1] = (fun(j, x, y, t, u + h, v, w) - fun(j, x, y, t, u, v, w)) / h;
			}
			for (int j = 0; j < 4; j++)
			{
				J[j][2] = (fun(j, x, y, t, u, v + h, w) - fun(j, x, y, t, u, v, w)) / h;
			}
			for (int j = 0; j < 4; j++)
			{
				J[j][3] = (fun(j, x, y, t, u, v, w + h) - fun(j, x, y, t, u, v, w)) / h;
			}

			vector<double> delta_sol(solve_pivot_Gauss(J, -F));
			sol = sol + delta_sol;
			t = sol[0];
			u = sol[1];
			v = sol[2];
			w = sol[3];
			if (abs(delta_sol) / abs(sol) <= epsilon)
				return;
			if (k >= M)
				cout << "Fail to solve the equaitons : the iterations are not enough." << endl;
			h = 1e-11;
			//h = abs(F);
		}
	}
}



//分片二次代数插值求函数g(t, u)的值以和f(x, y)做比较
double lx(const int &i, const int &k, const double &x)
{
	double res = 1;
	for (int t = i - 1; t <= i + 1; t++)
	{
		if (t == k)
			continue;
		res *= (x - 0.2 * t);
		res /= (0.2 * k - 0.2 * t);
	}
	return res;
}

double ly(const int &j, const int &r, const double &y)
{
	double res = 1;
	for (int t = j - 1; t <= j + 1; t++)
	{
		if (t == r)
			continue;
		res *= (y - 0.4 * t);
		res /= (0.4 * r - 0.4 * t);
	}
	return res;
}

double fun_tu(const double &t, const double &u)
{
	//生成U矩阵
	Matrix U(6, 6);
	U[0][0] = -0.5;  U[0][1] = -0.34; U[0][2] = 0.14;  U[0][3] = 0.94;  U[0][4] = 2.06;  U[0][5] = 3.5;
	U[1][0] = -0.42; U[1][1] = -0.5;  U[1][2] = -0.26; U[1][3] = 0.3;   U[1][4] = 1.18;  U[1][5] = 2.38;
	U[2][0] = -0.18; U[2][1] = -0.5;  U[2][2] = -0.5;  U[2][3] = -0.18; U[2][4] = 0.46;  U[2][5] = 1.42;
	U[3][0] = 0.22;  U[3][1] = -0.34; U[3][2] = -0.58; U[3][3] = -0.5;  U[3][4] = -0.1;  U[3][5] = 0.62;
	U[4][0] = 0.78;  U[4][1] = -0.02; U[4][2] = -0.5;  U[4][3] = -0.66; U[4][4] = -0.5;  U[4][5] = -0.02;
	U[5][0] = 1.5;   U[5][1] = 0.46;  U[5][2] = -0.26; U[5][3] = -0.66; U[5][4] = -0.74; U[5][5] = -0.5;

	//判断输入的t, u在哪个片区，决定i, j的值
	int i, j;
	if (t <= 0.3)
	{
		i = 1;
	}
	if (t > 0.7)
	{
		i = 4;
	}
	if (u <= 0.6)
	{
		j = 1;
	}
	if (u > 1.4)
	{
		j = 4;
	}
	for (int p = 2; p <= 3; p++)
	{
		if (t <= 0.2 * p + 0.1 && t > 0.2 * p - 0.1)
		{
			i = p;
		}
	}
	for (int q = 2; q <= 3; q++)
	{
		if (u <= 0.4 * q + 0.2 && u > 0.4 * q - 0.2)
		{
			j = q;
		}
	}

	double res = 0;
	for (int k = i - 1; k <= i + 1; k++)
	{
		for (int r = j - 1; r <= j + 1; r++)
		{
			res += lx(i, k, t) * ly(j, r, u) * U[k][r];
		}
	}
	return res;
}


double fun_p_xy(const double &x, const double &y, const int &k, const Matrix &C)
{
	double res = 0;
	for (int r = 0; r < k + 1; r++)
	{
		for (int s = 0; s < k + 1; s++)
		{
			res += C[r][s] * pow(x, r) * pow(y, s);
		}
	}
	return res;
}


void algorithm_forThesSecondQuesion()
{
	//算法初始化
	vector<double> var(6);
	double &x = var[0];
	double &y = var[1];
	double &t = var[2];
	double &u = var[3];
	double &v = var[4];
	double &w = var[5];
	double epsilon = 1e-12;
	const int MAX_INTERATIONS = 100000;

	double k_opt;
	Matrix C_opt;

	//利用曲面拟合迭代求解最小的 k 值满足sigma <= 10^(-7)
	Matrix U(11, 21);
	//k 值最大为7，否则为插值型多项式。
	for (int k = 1; k < 7; k++)
	{
		//矩阵 U 的初始化
		for (int i = 0; i <= 10; i++)
		{
			for (int j = 0; j <= 20; j++)
			{
				//迭代求解非线性方程组的输入矩阵初始化
				x = 0.08 * i;
				y = 0.5 + 0.05 * j;
				t = 0.5;
				u = 1;
				v = 1;
				w = 1;

				solveNonlinearEquations(1, var, epsilon, MAX_INTERATIONS);
				U[i][j] = fun_tu(t, u);
			}
		}

		//矩阵 B 的初始化
		//选取正交多项式系为{x^i}i = 0 ~ n + 1型。
		Matrix B(11, k + 1);
		for (int i = 0; i < 11; i++)
		{
			for (int r = 0; r < k + 1; r++)
			{
				//计算phai_r(x_i)
				B[i][r] = pow(0.08 * i, r);
			}
		}

		//矩阵 G 的初始化
		//选取正交多项式系为{y^i}i = 0 ~ n + 1型。
		Matrix G(21, k + 1);
		for (int j = 0; j < 21; j++)
		{
			for (int s = 0; s < k + 1; s++)
			{
				//计算pusai_s(y_j)
				G[j][s] = pow((0.5 + 0.05 * j), s);
			}
		}

		Matrix C(inverse(trans(B) * B) * trans(B) * U * G * inverse(trans(G) * G));

		double sigma = 0;

		for (int i = 0; i < 11; i++)
		{
			for (int j = 0; j < 21; j++)
			{
				//计算 sigma 值
				sigma += pow((U[i][j] - fun_p_xy(0.08 * i, (0.5 + 0.05 * j), k, C)), 2) ;
			}
		}

		cout << "k = " << std::fixed << k << ", sigma(" << k << ") = " << std::setprecision(12) << std::scientific << sigma << endl;

		if (sigma <= 1e-7)
		{
			cout << "此时 k 值达到要求，k 和 sigma 值如上。" << endl;
			cout << endl;
			//输出数表（x_i, y_j, f(x_i, y_j)）
			cout << "数表（x_i, y_j, f(x_i, y_j)），(i = 0...10, j = 0...20) 为：" << endl << endl;
			
			for (int i = 0; i < 11; i++)
			{
				for (int j = 0; j < 21; j++)
				{
					cout <<  "( x_" << i << std::setw(2) << " = " << std::fixed << std::setprecision(2) << 0.08 * i << ", y_" << std::setw(2) << j << " = " << (0.5 + 0.05 * j);
					cout << ", f(x_" << std::setw(2) << i << ", y_" << std::setw(2) << j << " ) = " << std::setprecision(12) << std::scientific << U[i][j] << " )" << endl;
				}
				cout << endl;
			}
			//此时算法符合要求，将 k 的值 和矩阵 C 保留下来。
			C_opt = C;
			k_opt = k;
			cout << "矩阵 C 为：" << endl << endl;
			printMatrix(C);
			cout << endl;
			break;
		}
	}


	//利用x, y的值反解t, u
	cout << "数表（xi*, yj*, f(xi*, yj*), p(xi*, yj*) )，(i = 1...8, j = 1...5) 为：" << endl << endl;
	for (int i = 1; i <= 8; i++)
	{
		for (int j = 1; j <= 5; j++)
		{
			//迭代求解非线性方程组的输入矩阵初始化
			x = 0.1 * i;
			y = 0.5 + 0.2 * j;
			t = 0.5;
			u = 1;
			v = 1;
			w = 1;

			solveNonlinearEquations(1, var, epsilon, MAX_INTERATIONS);
			//输出数表（xi*, yj*, f(xi*, yj*), p(xi*, yj*) )
			cout << "( x_" << i << "* = " << std::fixed  << std::setprecision(2) << 0.1 * i << ", y_" << j << "* = " << (0.5 + 0.2 * j);
			cout << ", f(x_" << i << "*, y_" << j << "* ) = " << std::setprecision(12) << std::scientific << fun_tu(t, u);
			cout << ", p(x_" << i << "*, y_" << j << "* ) = ";
			cout << fun_p_xy(0.1 * i, (0.5 + 0.2 * j), k_opt, C_opt) << " )" << endl;
		}
		cout << endl;
	}
}