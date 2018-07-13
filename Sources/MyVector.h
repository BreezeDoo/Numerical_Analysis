//MyVector.h，提供和vector相关的各种运算，目前包括向量点成，矩阵乘向量，压缩矩阵乘向量等。
#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <assert.h>
#include "Matrix.h"

using std::vector;
using std::string;
using std::cout;
using std::cin;
using std::cerr;
using std::endl;

template <typename Type>
vector<Type> operator+(const vector<Type>& x, const vector<Type>& y)
{
	assert(x.size() == y.size());
	int size = x.size();
	vector<Type> result(size);
	for (int i = 0; i < size; i++)
	{
		result[i] = x[i] + y[i];
	}
	return result;
}

template <typename Type>
vector<Type> operator-(const vector<Type>& x)
{
	int size = x.size();
	vector<Type> y(size);
	for (int i = 0; i < size; i++)
	{
		y[i] = -x[i];
	}
	return y;
}

template <typename Type>
vector<Type> operator-(const vector<Type>& x, const vector<Type>& y)
{
	assert(x.size() == y.size());
	int size = x.size();
	vector<Type> result(size);
	for (int i = 0; i < size; i++)
	{
		result[i] = x[i] - y[i];
	}
	return result;
}

template <typename Type1, typename Type2>
vector<Type1> operator*(const vector<Type1>& x, const Type2& n)
{
	vector<Type1> y(x.size());
	for (size_t i = 0; i < x.size(); ++i)
		y[i] = x[i] * n;
	return y;
}

template <typename Type>
Type operator*(const vector<Type>& x, const vector<Type>& y)
{
	assert(x.size() == y.size());
	int size = x.size();
	Type sum = 0.;
	for (size_t i = 0; i < size; ++i)
		sum += x[i] * y[i];
	return sum;
}

template <typename Type>
vector<Type> operator*(const MATRIX<Type>& m, const vector<Type>& v)
{
	assert(m.cols() == v.size());
	int rows = m.rows();
	int size = v.size();
	vector<Type> result(size, 0);
	for (int j = 0; j < rows; j++)
	{
		for (int i = 0; i < size; i++)
		{
			result[j] += m[j][i] * v[i];
		}
	}
	return result;
}


vector<double> operator/(const vector<double>& v, const double& x)
{
	int size = v.size();
	vector <double> result(size, 0);
		for (int i = 0; i < size; i++)
		{
			result[i] = v[i] / x;
		}
	return result;
}

int bandMatrix_trans_i(const int &input_i, const int &input_j, const int &s)
{
	return input_i - input_j + s + 1;
}

vector<double> bandMatrix_multiply_vector(const Matrix &m_input, const int &s, const int &r, const vector<double> &v_input)
{
	assert(m_input.cols() == v_input.size());
	int size = v_input.size();
	vector<double> sol(size, 0.);
	for (int i = 1; i <= r + 1; i++)
	{
		for (int j = 1; j <= s + i; j++)
		{
			sol[i - 1] += m_input[bandMatrix_trans_i(i, j, s) - 1][j - 1] * v_input[j - 1];
		}
	}
	for (int i = r + 2; i <= size - s; i++)
	{
		for (int j = i - r; j <= i + s; j++)
		{
			sol[i - 1] += m_input[bandMatrix_trans_i(i, j, s) - 1][j - 1] * v_input[j - 1];
		}
	}
	for (int i = size - s + 1; i <= size; i++)
	{
		for (int j = i - r; j <= size; j++)
		{
			sol[i - 1] += m_input[bandMatrix_trans_i(i, j, s) - 1][j - 1] * v_input[j - 1];
		}
	}
	return sol;
}


Matrix vector_multiply(const vector<double> &v1, const vector<double> &v2)
{
	int rows = v1.size();
	int cols = v2.size();
	Matrix m(rows, cols);
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			m[i][j] = v1[i] * v2[j];
		}
	}
	return m;
}