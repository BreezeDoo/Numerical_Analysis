//Linear_Equatinos.hͷ�ļ����������Է������һ��ⷨ��
//����˳���˹��ȥ����ѡ��Ԫ�ĸ�˹��ȥ����˳��doolittle�ֽⷨ��ѡ��Ԫ��doolittl�ֽⷨ��
//��״�����˳���˹��ȥ���ȡ�
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
    //�����������
	vector<double> readVector(istream& in);
	vector<double> readVector(string file);
    void printVector(const vector<double> &v, ostream& out = std::cout);
    void printVector(const vector<double> &v, string file);

    //�����к���
    vector<double> Substitution_upper(const Matrix &A); //�ش���������Ҫ���������������ϵ�����Ϊ��������
    vector<double> Substitution_lower(const Matrix &A); //�ش���������Ҫ���������������ϵ�����Ϊ��������		

	vector<double> Substition_Doolittle(const Matrix &A, vector<double> right_end_vec);
	//�㷨����
    Matrix LU_sequential_resolve(const Matrix &A);
	vector<double> solve_sequentioal_LU(const Matrix &m_input, const vector<double> &v_input);
}


/***********************************************/
/*                  �������                    */
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
/*                  �����к���                  */
/***********************************************/

//�ش����̣���ϵ����Ϊ����������������ش������ؽ�����
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
//    printVector(v);//���������
    return v;
}

//�ش����̣���ϵ����Ϊ����������������ش������ؽ�����
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
//    printVector(v);//���������
    return v;
}

/***********************************************/
/*                  �㷨ʵ��                    */
/***********************************************/

//doolittle_LU�ֽⷨ
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
	for (int i = 0; i < m_input.rows(); i++)//��LU����ĶԽ��߻�Ϊ1���Ա���������������Ԫ����
		m[i][i] = 1;
	m = Augmented_matrix_generating(m_input, Substitution_lower(m));//����Ԫ��������µľ����Ա�������������������Ա���Ԫ
	return Substitution_upper(m);
}