//Linear_Equatinos.h头文件，包含线性方程组的一般解法，
//包括顺序高斯消去法，选主元的高斯消去法，顺序doolittle分解法，选主元的doolittl分解法，
//带状矩阵的顺序高斯消去法等。
#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iomanip>
#include <assert.h>
#include "Matrix.h"
#include <algorithm>

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

namespace Linear_Equations
{
    //输入输出函数
	vector<double> readVector(istream& in);
	vector<double> readVector(string file);
    void printVector(const vector<double> &v, ostream& out = std::cout);
    void printVector(const vector<double> &v, string file);

    //过程中函数
    vector<double> Substitution_upper(const Matrix &A); //回带函数——要求输入增广矩阵，且系数阵必为上三角阵
    vector<double> Substitution_lower(const Matrix &A); //回带函数——要求输入增广矩阵，且系数阵必为下三角阵		

	vector<double> Substition_Doolittle(const Matrix &A, vector<double> right_end_vec);
	//算法函数
    Matrix LU_sequential_resolve(const Matrix &A);
	vector<double> solve_sequentioal_LU(const Matrix &m_input, const vector<double> &v_input);
}


/***********************************************/
/*                  输入输出                    */
/***********************************************/


vector<double> Linear_Equations::readVector(istream& in)
{
	Matrix m;
	m = readMatrix(cin);
	vector<double> v(m.rows());
	for (int i = 0; i < m.rows(); i++)
		v[i] = m[i][0];
	return v;
}


vector<double> Linear_Equations::readVector(string file)
{
	Matrix m;
	m = readMatrix(file);
	vector<double> v(m.rows());
	for (int i = 0; i < m.rows(); i++)
		v[i] = m[i][0];
	return v;
}


void Linear_Equations::printVector(const vector<double> &v, ostream &out)
{
    Matrix m;
    m.push_back(v);
    printMatrix(m, out);
}


void Linear_Equations::printVector(const vector<double> &v, string file)
{
    Matrix m;
    m.push_back(v);
    printMatrix(m, file);
}

/***********************************************/
/*                  过程中函数                  */
/***********************************************/

//回带过程：将系数阵为上三角阵的增广矩阵回带并返回解向量
vector<double> Linear_Equations::Substitution_upper(const Matrix &A)
{
    assert(isAugmented(A));
    int rows = A.rows();
    int cols = A.cols();
    vector<double> v(rows);

    for (int k = rows - 1; k >= 0; --k)
    {
        v[k] = A[k][cols - 1];
        for (int i = rows - 1; i > k; --i)
            v[k] = v[k] - A[k][i] * v[i];
        v[k] = v[k] / A[k][k];
    }
//    printVector(v);//检验解向量
    return v;
}

//回带过程：将系数阵为下三角阵的增广矩阵回带并返回解向量
vector<double> Linear_Equations::Substitution_lower(const Matrix &A)
{
    assert(isAugmented(A));
    int rows = A.rows();
    int cols = A.cols();
    vector<double> v(rows);

    for (int k = 0; k <= rows - 1; ++k)
    {
        v[k] = A[k][cols - 1];
        for (int i = 0; i < k; ++i)
            v[k] = v[k] - A[k][i] * v[i];
        v[k] = v[k] / A[k][k];
    }
//    printVector(v);//检验解向量
    return v;
}

/***********************************************/
/*                  算法实现                    */
/***********************************************/

//doolittle_LU分解法
Matrix Linear_Equations::LU_sequential_resolve(const Matrix &m_input)
{
	assert(m_input.square());
	Matrix m(m_input);
	int rows = m.rows();
	for (int k = 0; k < rows; k++)
	{
		for (int j = k; j < rows && k > 0; j++)
		{
				for (int s = 0; s <= k - 1; s++)
					m[k][j] -= m[k][s] * m[s][j];
		}
		for (int i = k + 1; i <= rows - 1 && k < rows - 1; i++)
		{
			for (int s = 0; s < k; s++)
				m[i][k] -= m[i][s] * m[s][k];
			m[i][k] /= m[k][k];
		}
	}
	return m;
}


vector<double> Linear_Equations::solve_sequentioal_LU(const Matrix &m_input, const vector<double> &v_input)
{
	Matrix m_seqLU;
	m_seqLU = LU_sequential_resolve(m_input);
	vector<double> sol_seqLU;
	return Substition_Doolittle(m_seqLU, v_input);
}

vector<double> Linear_Equations::Substition_Doolittle(const Matrix &m_input, vector<double> right_end_vec)
{
	Matrix m;
	m = Augmented_matrix_generating(m_input, right_end_vec);
	for (int i = 0; i < m_input.rows(); i++)//将LU矩阵的对角线化为1，以便进行下三角阵的消元过程
		m[i][i] = 1;
	m = Augmented_matrix_generating(m_input, Substitution_lower(m));//将消元结果传给新的矩阵以便生成上三角增广矩阵以便消元
	return Substitution_upper(m);
}