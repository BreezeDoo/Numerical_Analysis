//Matrix.h头文件，包含矩阵的基本运算。
#pragma once
#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iomanip>  //用于设置输出格式

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

// 任意类型矩阵类
template <typename Object>
class MATRIX
{
public:
	explicit MATRIX() : array(0) {}

	MATRIX(int rows, int cols) :array(rows)
	{
		for (int i = 0; i < rows; ++i)
		{
			array[i].resize(cols);
		}
	}

	MATRIX(const MATRIX<Object>& m) { *this = m; }

	void resize(int rows, int cols);           // 改变当前矩阵大小
	bool push_back(const vector<Object>& v);   // 在矩阵末尾添加一行数据
	void swap_row(int row1, int row2);         // 交换两行的数据

	int  rows() const { return array.size(); }
	int  cols() const { return rows() ? (array[0].size()) : 0; }
	bool empty() const { return rows() == 0; }        // 是否为空
	bool square() const { return (!(empty()) && rows() == cols()); }  // 是否为方阵


	const vector<Object>& operator[](int row) const { return array[row]; } //[]操作符重载
	vector<Object>& operator[](int row) { return array[row]; }

	const Object& at(int row, int col) const { return array[row - 1][col - 1]; }
	Object& at(int row, int col) { return array[row - 1][col - 1]; }

protected:
	vector< vector<Object> > array;
};

// 改变当前矩阵大小
template <typename Object>
void MATRIX<Object>::resize(int rows, int cols)
{
	int rs = this->rows();
	int cs = this->cols();

	if (rows == rs && cols == cs)
	{
		return;
	}
	else if (rows == rs && cols != cs)
	{
		for (int i = 0; i < rows; ++i)
		{
			array[i].resize(cols);
		}
	}
	else if (rows != rs && cols == cs)
	{
		array.resize(rows);
		for (int i = rs; i < rows; ++i)
		{
			array[i].resize(cols);
		}
	}
	else
	{
		array.resize(rows);
		for (int i = 0; i < rows; ++i)
		{
			array[i].resize(cols);
		}
	}
}

// 在矩阵末尾添加一行
template <typename Object>
bool MATRIX<Object>::push_back(const vector<Object>& v)
{
	if (rows() == 0 || cols() == (int)v.size())
	{
		array.push_back(v);
	}
	else
	{
		return false;
	}

	return true;
}

// 交换两行
template <typename Object>
void MATRIX<Object>::swap_row(int row1, int row2)
{
	if (row1 != row2 && row1 >= 0 &&
		row1 < rows() && row2 >= 0 && row2 < rows())
	{
		vector<Object>& v1 = array[row1];
		vector<Object>& v2 = array[row2];
		vector<Object> tmp = v1;
		v1 = v2;
		v2 = tmp;
	}
}

// 矩阵转置
template <typename Object>
const MATRIX<Object> trans(const MATRIX<Object>& m)
{
	MATRIX<Object> ret;
	if (m.empty()) return ret;

	int row = m.cols();
	int col = m.rows();
	ret.resize(row, col);

	for (int i = 0; i < row; ++i)
	{
		for (int j = 0; j < col; ++j)
		{
			ret[i][j] = m[j][i];
		}
	}

	return ret;
}

//////////////////////////////////////////////////////////
// double类型矩阵类，用于科学计算
// 继承自MATRIX类
// 实现常用操作符重载，并实现计算矩阵的行列式、逆以及LU分解
class Matrix :public MATRIX<double>
{
public:
	Matrix() :MATRIX<double>() {}
	Matrix(int c, int r) :MATRIX<double>(c, r) {}
	Matrix(const Matrix& m) { *this = m; }

	const Matrix& operator+=(const Matrix& m);
	const Matrix& operator-=(const Matrix& m);
	const Matrix& operator*=(const Matrix& m);
	const Matrix& operator/=(const Matrix& m);
};

bool  operator==(const Matrix& lhs, const Matrix& rhs);        // 重载操作符==
bool  operator!=(const Matrix& lhs, const Matrix& rhs);        // 重载操作符!=
const Matrix operator+(const Matrix& lhs, const Matrix& rhs);  // 重载操作符+
const Matrix operator-(const Matrix& lhs, const Matrix& rhs);  // 重载操作符-
const Matrix operator*(const Matrix& lhs, const Matrix& rhs);  // 重载操作符*
const Matrix operator*(const double& s, const Matrix& rhs);    // 重载操作符*
const Matrix operator/(const Matrix& lhs, const Matrix& rhs);  // 重载操作符/
const bool isAugmented(const Matrix &A);						//检测函数
const double det(const Matrix& m);                             // 计算行列式
const double det(const Matrix& m, int start, int end);         // 计算子矩阵行列式
const Matrix abs(const Matrix& m);                             // 计算所有元素的绝对值
const double max(const Matrix& m);                             // 所有元素的最大值
const double max(const Matrix& m, int& row, int& col);          // 所有元素中的最大值及其下标
const double min(const Matrix& m);                             // 所有元素的最小值
const double min(const Matrix& m, int& row, int& col);          // 所有元素的最小值及其下标
const Matrix trans(const Matrix& m);                           // 返回转置矩阵
const Matrix submatrix(const Matrix& m, int rb, int re, int cb, int ce);  // 返回子矩阵
const Matrix inverse(const Matrix& m);                         // 计算逆矩阵
const Matrix displacement(const Matrix &m_input, const double &p);		//返回平移量为 p 的平移后矩阵
const Matrix LU(const Matrix& m);                              // 计算方阵的LU分解
const Matrix readMatrix(istream& in = std::cin);               // 从指定输入流读入矩阵
const Matrix readMatrix(string file);                          // 从文本文件读入矩阵
const Matrix loadMatrix(string file);                          // 从二进制文件读取矩阵
void  printMatrix(const Matrix& m, ostream& out = std::cout);  // 从指定输出流打印矩阵
void  printMatrix(const Matrix& m, string file);                // 将矩阵输出到文本文件
void  saveMatrix(const Matrix& m, string file);                 // 将矩阵保存为二进制文件

const Matrix& Matrix::operator+=(const Matrix& m)
{
	if (rows() != m.rows() || rows() != m.cols())
	{
		return *this;
	}

	int r = rows();
	int c = cols();

	for (int i = 0; i < r; ++i)
	{
		for (int j = 0; j < c; ++j)
		{
			array[i][j] += m[i][j];
		}
	}

	return *this;
}


const Matrix& Matrix::operator-=(const Matrix& m)
{
	if (rows() != m.rows() || cols() != m.cols())
	{
		return *this;
	}

	int r = rows();
	int c = cols();

	for (int i = 0; i < r; ++i)
	{
		for (int j = 0; j < c; ++j)
		{
			array[i][j] -= m[i][j];
		}
	}

	return *this;
}

const Matrix& Matrix::operator*=(const Matrix& m)
{
	if (cols() != m.rows() || !m.square())
	{
		return *this;
	}

	Matrix ret(rows(), cols());

	int r = rows();
	int c = cols();

	for (int i = 0; i < r; ++i)
	{
		for (int j = 0; j < c; ++j)
		{
			double sum = 0.0;
			for (int k = 0; k < c; ++k)
			{
				sum += array[i][k] * m[k][j];
			}
			ret[i][j] = sum;
		}
	}

	*this = ret;
	return *this;
}

const Matrix& Matrix::operator/=(const Matrix& m)
{
	Matrix tmp = inverse(m);
	return operator*=(tmp);
}


bool operator==(const Matrix& lhs, const Matrix& rhs)
{
	if (lhs.rows() != rhs.rows() || lhs.cols() != rhs.cols())
	{
		return false;
	}

	for (int i = 0; i < lhs.rows(); ++i)
	{
		for (int j = 0; j < lhs.cols(); ++j)
		{
			if (rhs[i][j] != rhs[i][j])
			{
				return false;
			}
		}
	}

	return true;
}

bool operator!=(const Matrix& lhs, const Matrix& rhs)
{
	return !(lhs == rhs);
}

const Matrix operator+(const Matrix& lhs, const Matrix& rhs)
{
	Matrix m;
	if (lhs.rows() != rhs.rows() || lhs.cols() != rhs.cols())
	{
		return m;
	}

	m = lhs;
	m += rhs;

	return m;
}

const Matrix operator-(const Matrix& lhs, const Matrix& rhs)
{
	Matrix m;
	if (lhs.rows() != rhs.rows() || lhs.cols() != rhs.cols())
	{
		return m;
	}

	m = lhs;
	m -= rhs;

	return m;
}

const Matrix operator*(const Matrix& lhs, const Matrix& rhs)
{
	Matrix m;
	if (lhs.cols() != rhs.rows())
	{
		return m;
	}

	m.resize(lhs.rows(), rhs.cols());

	int r = m.rows();
	int c = m.cols();
	int K = lhs.cols();

	for (int i = 0; i < r; ++i)
	{
		for (int j = 0; j < c; ++j)
		{
			double sum = 0.0;
			for (int k = 0; k < K; ++k)
			{
				sum += lhs[i][k] * rhs[k][j];
			}
			m[i][j] = sum;
		}
	}

	return m;
}

const Matrix operator*(const double& s, const Matrix& rhs)
{
	Matrix m;
	if (rhs.empty())
		return m;
	int r = rhs.rows();
	int c = rhs.cols();
	m.resize(r, c);
	for (int i = 0; i < r; i++)
	{
		for (int j = 0; j < c; j++)
		{
			m[i][j] = rhs[i][j] * s;
		}
	}
	return m;
}

const Matrix operator/(const Matrix& lhs, const Matrix& rhs)
{
	Matrix tmp = inverse(rhs);
	Matrix m;

	if (tmp.empty())
	{
		return m;
	}

	return m = lhs * tmp;
}

inline static double LxAbs(double d)
{
	return (d >= 0) ? (d) : (-d);
}

inline
static bool isSignRev(const vector<double>& v)
{
	int p = 0;
	int sum = 0;
	int n = (int)v.size();

	for (int i = 0; i < n; ++i)
	{
		p = (int)v[i];
		if (p >= 0)
		{
			sum += p + i;
		}
	}

	if (sum % 2 == 0) // 如果是偶数，说明不变号
	{
		return false;
	}
	return true;
}

const bool isAugmented(const Matrix &A)
{
	if (A.rows() + 1 != A.cols())
	{
		cout << "输入矩阵非增广矩阵" << endl;
		return false;
	}
	return true;
}

// 计算方阵行列式
const double det(const Matrix& m)
{
	double ret = 0.0;

	if (m.empty() || !m.square()) return ret;

	Matrix N = LU(m);

	if (N.empty()) return ret;

	ret = 1.0;
	for (int i = 0; i < N.cols(); ++i)
	{
		ret *= N[i][i];
	}

	if (isSignRev(N[N.rows() - 1]))
	{
		return -ret;
	}

	return ret;
}

// 计算矩阵指定子方阵的行列式
const double det(const Matrix& m, int start, int end)
{
	return det(submatrix(m, start, end, start, end));
}


// 计算矩阵转置
const Matrix trans(const Matrix& m)
{
	Matrix ret;
	if (m.empty()) return ret;

	int r = m.cols();
	int c = m.rows();

	ret.resize(r, c);
	for (int i = 0; i < r; ++i)
	{
		for (int j = 0; j < c; ++j)
		{
			ret[i][j] = m[j][i];
		}
	}

	return ret;
}

// 计算逆矩阵
const Matrix  inverse(const Matrix& m)
{
	Matrix ret;

	if (m.empty() || !m.square())
	{
		return ret;
	}

	int n = m.rows();

	ret.resize(n, n);
	Matrix A(m);

	for (int i = 0; i < n; ++i) ret[i][i] = 1.0;

	for (int j = 0; j < n; ++j)  //每一列
	{
		int p = j;
		double maxV = LxAbs(A[j][j]);
		for (int i = j + 1; i < n; ++i)  // 找到第j列中元素绝对值最大行
		{
			if (maxV < LxAbs(A[i][j]))
			{
				p = i;
				maxV = LxAbs(A[i][j]);
			}
		}

		if (maxV < 1e-20)
		{
			ret.resize(0, 0);
			return ret;
		}

		if (j != p)
		{
			A.swap_row(j, p);
			ret.swap_row(j, p);
		}

		double d = A[j][j];
		for (int i = j; i < n; ++i) A[j][i] /= d;
		for (int i = 0; i < n; ++i) ret[j][i] /= d;

		for (int i = 0; i < n; ++i)
		{
			if (i != j)
			{
				double q = A[i][j];
				for (int k = j; k < n; ++k)
				{
					A[i][k] -= q * A[j][k];
				}
				for (int k = 0; k < n; ++k)
				{
					ret[i][k] -= q * ret[j][k];
				}
			}
		}
	}

	return ret;
}

//返回平移量为 p 的平移后矩阵
const Matrix displacement(const Matrix &m_input, const double &p)
{
//	assert(m_input.square());
	int size = m_input.cols();
	Matrix m(m_input);
	for (int i = 0; i < size; i++)
	{
		m[2][i] -= p;
	}
	return m;
}

// 计算绝对值
const Matrix abs(const Matrix& m)
{
	Matrix ret;

	if (m.empty())
	{
		return ret;
	}

	int r = m.rows();
	int c = m.cols();
	ret.resize(r, c);

	for (int i = 0; i < r; ++i)
	{
		for (int j = 0; j < c; ++j)
		{
			double t = m[i][j];
			if (t < 0) ret[i][j] = -t;
			else ret[i][j] = t;
		}
	}

	return ret;
}

// 返回矩阵所有元素的最大值
const double max(const Matrix& m)
{
	if (m.empty()) return 0.;

	double ret = m[0][0];
	int r = m.rows();
	int c = m.cols();

	for (int i = 0; i < r; ++i)
	{
		for (int j = 0; j < c; ++j)
		{
			if (m[i][j] > ret) ret = m[i][j];
		}
	}
	return ret;
}

// 计算矩阵最大值，并返回该元素的引用
const double max(const Matrix& m, int& row, int& col)
{
	if (m.empty()) return 0.;

	double ret = m[0][0];
	row = 0;
	col = 0;

	int r = m.rows();
	int c = m.cols();

	for (int i = 0; i < r; ++i)
	{
		for (int j = 0; j < c; ++j)
		{
			if (m[i][j] > ret)
			{
				ret = m[i][j];
				row = i;
				col = j;
			}
		}
	}
	return ret;
}

// 计算矩阵所有元素最小值
const double min(const Matrix& m)
{
	if (m.empty()) return 0.;

	double ret = m[0][0];
	int r = m.rows();
	int c = m.cols();

	for (int i = 0; i < r; ++i)
	{
		for (int j = 0; j < c; ++j)
		{
			if (m[i][j] > ret) ret = m[i][j];
		}
	}

	return ret;
}

// 计算矩阵最小值，并返回该元素的引用
const double min(const Matrix& m, int& row, int& col)
{
	if (m.empty()) return 0.;

	double ret = m[0][0];
	row = 0;
	col = 0;
	int r = m.rows();
	int c = m.cols();

	for (int i = 0; i < r; ++i)
	{
		for (int j = 0; j < c; ++j)
		{
			if (m[i][j] > ret)
			{
				ret = m[i][j];
				row = i;
				col = j;
			}
		}
	}

	return ret;
}

// 取矩阵中指定位置的子矩阵,输入方式同数学用法，下标从 1 开始
const Matrix submatrix(const Matrix& m, int rb, int re, int cb, int ce)
{
	Matrix ret;
	if (m.empty()) return ret;

	if (rb < 1 || re > m.rows() || rb > re) return ret;
	if (cb < 1 || ce > m.cols() || cb > ce) return ret;

	ret.resize(re - rb + 1, ce - cb + 1);

	for (int i = rb; i <= re; ++i)
	{
		for (int j = cb; j <= ce; ++j)
		{
			ret[i - rb][j - cb] = m[i - 1][j - 1];
		}
	}

	return ret;
}


inline static
int max_idx(const Matrix& m, int k, int n)
{
	int p = k;
	for (int i = k + 1; i < n; ++i)
	{
		if (LxAbs(m[p][k]) < LxAbs(m[i][k]))
		{
			p = i;
		}
	}
	return p;
}

// 计算方阵 M 的 LU 分解
// 其中L为对角线元素全为1的下三角阵，U为对角元素依赖M的上三角阵
// 使得 M = LU
// 返回矩阵下三角部分存储L(对角元素除外)，上三角部分存储U(包括对角线元素)
const Matrix LU(const Matrix& m)
{
	Matrix ret;

	if (m.empty() || !m.square()) return ret;

	int n = m.rows();
	ret.resize(n + 1, n);

	for (int i = 0; i < n; ++i)
	{
		ret[n][i] = -1.0;
	}

	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			ret[i][j] = m[i][j];
		}
	}

	for (int k = 0; k < n - 1; ++k)
	{
		int p = max_idx(ret, k, n);
		if (p != k)              // 进行行交换
		{
			ret.swap_row(k, p);
			ret[n][k] = (double)p; // 记录交换信息
		}

		if (ret[k][k] == 0.0)
		{
			cout << "ERROR: " << endl;
			ret.resize(0, 0);
			return ret;
		}

		for (int i = k + 1; i < n; ++i)
		{
			ret[i][k] /= ret[k][k];
			for (int j = k + 1; j < n; ++j)
			{
				ret[i][j] -= ret[i][k] * ret[k][j];
			}
		}
	}

	return ret;
}

//---------------------------------------------------
//                      读取和打印
//---------------------------------------------------
// 从输入流读取矩阵
const Matrix readMatrix(istream& in)
{
	Matrix M;
	string str;
	double b;
	vector<double> v;

	while (getline(in, str))
	{
		for (string::size_type i = 0; i < str.size(); ++i)
		{
			if (str[i] == ',' || str[i] == ';')
			{
				str[i] = ' ';
			}
			else if (str[i] != '.' && (str[i] < '0' || str[i] > '9')
				&& str[i] != ' ' && str[i] != '\t' && str[i] != '-')
			{
				M.resize(0, 0);
				return M;
			}
		}

		istringstream sstream(str);
		v.resize(0);

		while (sstream >> b)
		{
			v.push_back(b);
		}
		if (v.size() == 0)
		{
			continue;
		}
		if (!M.push_back(v))
		{
			M.resize(0, 0);
			return M;
		}
	}

	return M;
}

// 从文本文件读入矩阵
const Matrix readMatrix(string file)
{
	ifstream fin(file.c_str());
	Matrix M;

	if (!fin)
	{
		cerr << "Error: open file " << file << " failed." << endl;
		return M;
	}

	M = readMatrix(fin);
	fin.close();

	return M;
}

// 将矩阵输出到指定输出流
void printMatrix(const Matrix& m, ostream& out)
{
	if (m.empty())
	{
		return;
	}

	int r = m.rows();
	int c = m.cols();

	int n = 0;              // 数据小数点前最大位数
	double maxV = max(abs(m));
	while (maxV >= 1.0)
	{
		maxV /= 10;
		++n;
	}
	if (n == 0) n = 1;    // 如果最大数绝对值小于1，这小数点前位数为1，为数字0
	int pre = 10;            // 小数点后数据位数
	int wid = n + pre + 7;  // 控制字符宽度=n+pre+符号位+小数点位

	out << std::showpoint;
	out << std::setiosflags(std::ios::fixed);
	out << std::setprecision(pre) << std::scientific;
	for (int i = 0; i < r; ++i)
	{
		for (int j = 0; j < c; ++j)
		{
			out << std::setw(wid) << m[i][j];
		}
		out << endl;
	}

	out << std::setprecision(6);
	out << std::noshowpoint;
}

// 将矩阵打印到指定文件
void printMatrix(const Matrix& m, string file)
{
	ofstream fout(file.c_str());
	if (!fout) return;

	printMatrix(m, fout);
	fout.close();
}

// 将矩阵数据存为二进制文件
void saveMatrix(const Matrix& m, string file)
{
	if (m.empty()) return;

	ofstream fout(file.c_str(), std::ios_base::out | std::ios::binary);
	if (!fout) return;

	int r = m.rows();
	int c = m.cols();
	char Flag[12] = "MATRIX_DATA";
	fout.write((char*)&Flag, sizeof(Flag));
	fout.write((char*)&r, sizeof(r));
	fout.write((char*)&c, sizeof(c));

	for (int i = 0; i < r; ++i)
	{
		for (int j = 0; j < c; ++j)
		{
			double t = m[i][j];
			fout.write((char*)&t, sizeof(t));
		}
	}

	fout.close();
}

// 从二进制文件load矩阵
const Matrix loadMatrix(string file)
{
	Matrix m;

	ifstream fin(file.c_str(), std::ios_base::in | std::ios::binary);
	if (!fin) return m;

	char Flag[12];
	fin.read((char*)&Flag, sizeof(Flag));

	string str(Flag);
	if (str != "MATRIX_DATA")
	{
		return m;
	}

	int r, c;
	fin.read((char*)&r, sizeof(r));
	fin.read((char*)&c, sizeof(c));

	if (r <= 0 || c <= 0) return m;

	m.resize(r, c);
	double t;

	for (int i = 0; i < r; ++i)
	{
		for (int j = 0; j < c; ++j)
		{
			fin.read((char*)&t, sizeof(t));
			m[i][j] = t;
		}
	}

	return m;
}

#endif