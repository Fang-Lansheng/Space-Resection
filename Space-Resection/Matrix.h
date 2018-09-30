#pragma once
#include "stdafx.h"
#include <iostream>
#include <math.h>
#include <assert.h>
#include <vector>

using namespace std;

class Matrix {

/// ����
public:
	int rows;				// �и�
	int cols;				// �п�
	int size;				// ������С
	vector<vector<double>> data;
	
/// ����
public:
	Matrix() {};											// Ĭ�Ϲ��캯��
	Matrix(int m, int n);									// ���캯��
	Matrix(int n);											// ���캯��������
	virtual ~Matrix();										// ��������
	double get(int row, int col) const;						// �õ������� row �� col �е�ֵ
	void set(int row, int col, double value);				// ������ row �� col �е�ֵ��Ϊ value

	bool operator=(const Matrix&);
	friend Matrix operator+(const Matrix&, const Matrix&);	// ��������� +
	friend Matrix operator-(const Matrix&, const Matrix&);	// ��������� -
	friend Matrix operator*(const Matrix&, const Matrix&);	// ��������� *
	friend Matrix operator*(double,		   const Matrix&);	// ��������� *
	friend Matrix operator*(const Matrix&, double);			// ��������� *
	friend Matrix operator/(const Matrix&, double);			// ��������� /

	double det();											// ���������ʽ
	Matrix t() const;										// �����ת��
	Matrix inv();											// �����������
	void print();											// ��ӡ����
};

Matrix::Matrix(int m, int n) {								// ���캯��

	assert(m > 0);
	assert(n > 0);

	rows = m;
	cols = n;
	size = m * n;

	vector<double> tempRow(n, 0);
	for (size_t i = 0; i < m; i++) {
		data.push_back(tempRow);
	}
}

Matrix::Matrix(int n) {										// ���캯��������

	assert(n > 0);

	rows = n;
	cols = n;
	size = n * n;

	vector<double> tempRow(n, 0);
	for (size_t i = 0; i < n; i++) {
		data.push_back(tempRow);
	}
}

Matrix::~Matrix() {											// ��������

	//delete[] p;
	//p = NULL;
}

double Matrix::get(int row, int col) const {				// get ��ȡ����Ԫ��ֵ

	assert(row >= 0 && row < this->rows);
	assert(col >= 0 && col < this->cols);

	return data[row][col];
}

void Matrix::set(int row, int col, double value) {			// set ���þ���Ԫ��ֵ

	assert(row >= 0 && row < rows);
	assert(col >= 0 && col < cols);

	data[row][col] = value;
}

bool Matrix::operator=(const Matrix& m) {

	assert(this->rows == m.rows);
	assert(this->cols == m.cols);

	for (int i = 0; i < this->rows; i++) {
		for (int j = 0; j < this->cols; j++) {
			// p[this->rows * i + j] = m.p[this->rows * i + j];
			data[i][j] = m.data[i][j];
		}
	}

	return true;
}

Matrix operator+(const Matrix& m1, const Matrix& m2) {		// �������

	assert(m1.rows == m2.rows);
	assert(m1.cols == m2.cols);

	Matrix ret(m1.rows, m1.cols);

	for (int i = 0; i < m1.rows; i++) {
		for (int j = 0; j < m1.cols; j++) {
			double val = m1.get(i, j) + m2.get(i, j);
			ret.set(i, j, val);
		}
	}
	return ret;
}

Matrix operator-(const Matrix& m1, const Matrix& m2) {		// �������

	assert(m1.rows == m2.rows);
	assert(m1.cols == m2.cols);

	Matrix ret(m1.rows, m1.cols);

	for (int i = 0; i < m1.rows; i++) {
		for (int j = 0; j < m1.cols; j++) {
			double val = m1.get(i, j) - m2.get(i, j);
			ret.set(i, j, val);
		}
	}
	return ret;
}

Matrix operator*(const Matrix& m1, const Matrix& m2) {		// �������
	
	assert(m1.size > 0);
	assert(m2.size > 0);
	assert(m1.cols == m2.rows);

	Matrix ret(m1.rows, m2.cols);
	double sum;

	for (int i = 0; i < m1.rows; i++) {
		for (int j = 0; j < m2.cols; j++) {
			sum = 0.0;
			for (int k = 0; k < m1.cols; k++) {				// m1.cols == m2.rows
				sum += m1.data[i][k] * m2.data[k][j];
			}
			ret.data[i][j] = sum;
		}
	}
	return ret;
}

vector<double> operator*(double n, vector<double> v) {

	vector<double> tempVec(v.size());
	for (size_t i = 0; i < (int)v.size(); i++) {
		tempVec[i] = v[i] * n;
	}

	return tempVec;
}

vector<double> operator/(vector<double> v, double n) {

	vector<double> tempVec(v.size());
	for (size_t i = 0; i < (int)v.size(); i++) {
		tempVec[i] = v[i] / n;
	}

	return tempVec;
}

Matrix operator*(double value, const Matrix& m1) {			// ���� �� ����

	Matrix ret(m1.rows, m1.cols);

	for (int i = 0; i < m1.size; i++) {
		// ret.p[i] = m1.p[i] * value;
		ret.data[i] = value * m1.data[i];
	}

	return ret;
}

Matrix operator*(const Matrix& m1, double value) {			// ���� �� ����

	Matrix ret(m1.rows, m1.cols);

	for (int i = 0; i < m1.size; i++) {
		// ret.p[i] = m1.p[i] * value;
		ret.data[i] = value * m1.data[i];
	}

	return ret;
}

Matrix operator/(const Matrix& m1, double value) {			// ���� �� ����

	Matrix ret(m1.rows, m1.cols);

	for (int i = 0; i < m1.size; i++) {
		// ret.p[i] = m1.p[i] / value;
		ret.data[i] = m1.data[i] / value;
	}

	return ret;
}

double dets(int n, vector<vector<double>> aa) {					// ���������ʽ

	if (n == 1)		
		return aa[0][0];	

	vector<vector<double>> bb((n - 1), vector<double>(n - 1));	// ���� n - 1 �׵Ĵ�������ʽ�� bb  
	int move = 0;												// �ж����Ƿ��ƶ�       	
	double sum = 0.0;											// sum Ϊ����ʽ��ֵ      	

	for (int a_row = 0; a_row < n; a_row++) {					// a �������Ѿ��� a(nn) ��ֵ�� b(n - 1)      	
		for (int b_row = 0; b_row < n - 1; b_row++) {			// �� aa ���һ�и�Ԫ�صĴ�������ʽ�浽bb      		
			move = a_row > b_row ? 0 : 1;						// bb ��С�� arow ���У�ͬ�и�ֵ�����ڵĴ�������ڵļ�һ      			
			for (int j = 0; j < n - 1; j++) {					// �� aa �ĵڶ��и�ֵ����n��      									
				bb[b_row][j] = aa[b_row + move][j + 1];
			}		
		}		
		int flag = (a_row % 2 == 0 ? 1 : -1);					// ��Ϊ����Ϊ 0������������ż��ʱ����������ʽΪ 1      		
		sum += flag * dets(n - 1, bb) * aa[a_row][0];			// aa ��һ�и�Ԫ�������������ʽ���ĺͼ�Ϊ����ʽ    	
	}	

	
	return sum;
}

double Matrix::det() {										// ���������ʽ
	assert(rows == cols);
	return dets(rows, data);
}

Matrix Matrix::t() const {									// ����ת��

	Matrix ret(cols, rows);
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			ret.data[j][i] = data[i][j];
		}
	}

	return ret;
}

Matrix Matrix::inv() {										// ��������
	
	double Det = det();
	assert(Det != 0.0);										// ����ʽ����Ϊ��

	Matrix tmp(rows);										// �൱������������߾���
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < rows; j++) {
			tmp.data[i][j] = data[i][j];
		}
	}
	Matrix ret(rows);										// ����󣨳�ʼ��Ϊ��λ����
	for (int i = 0; i < rows; i++) {
		ret.set(i, i, 1.0);
	}
	
	int r = 0;												// �к�
	double c = 1.0;											// ���г��ȱ任ʱ��ϵ�������ֵΪ 1

	// ���ȱ任��i ��ʾ�� i �У��������½�Ԫ�ػ�Ϊ 0 
	for (int i = 0; i < rows; i++) {
		if (tmp.data[i][i] == 0) {							// ����Խ����ϵ�Ԫ��Ϊ 0����ʱ�����б任	
			for (int j = i; j < rows; j++) {				// �� i ���Ժ�ĸ��н����жϣ��ҵ��� i ��Ԫ�ز�Ϊ 0 ���У������ i �н��н���
				if (j != i) {
					if (tmp.data[j][i] != 0) {
						r = j;								// ��ס���е��к�
						break;								// �˳�ѭ��
					}
				}
			}

			for (int k = 0; k < rows; k++) {				// �Ե� i �к͵� r �н��е���
				vector<vector<double>> t1(1, vector<double>(rows));
				t1[0][k] = tmp.data[i][k];
				tmp.data[i][k] = tmp.data[r][k];
				tmp.data[r][k] = t1[0][k];
				vector<vector<double>> t2(1, vector<double>(rows));
				t2[0][k] = ret.data[i][k];
				ret.data[i][k] = ret.data[r][k];
				ret.data[r][k] = t2[0][k];
			}
		}

		//for (int j = 0; j < rows; j++) {					// �������е������н��м��㣨�������еĸ��л�Ϊ 0��
		//	if (j != i) {									// ���в��������
		//		if (tmp.data[j][i] != 0) {					// ���� j �е� i ��Ϊ�㣬��������ʡ�Լ��㲽�裩
		//			c = tmp.data[j][i] / tmp.data[i][i];		// ����ϵ��
		//			for (int k = 0; k < rows; k++) {			// �� k �У��Ե� j �б�����
		//				tmp.data[j][k] -= c * tmp.data[i][k];	// �� j �� = �� j �� - c * �� i ��
		//				ret.data[j][k] -= c * ret.data[i][k];
		//			}
		//		}	
		//	}
		//}
		
		// ���Ž����������½ǡ���Ϊ 0 ����ʡȥ���������̣�Ҳ���������׳�����
		if (i > 0) {				// j -> �У� i -> ��
			for (int j = 0; j < i; j++) {							// �Ե� i �е� 0 - j �н��б������ֱ���ֵ��Ϊ 0
				if (tmp.data[i][j] != 0) {
					c = tmp.data[i][j] / tmp.data[j][j];
					for (int k = 0; k < rows; k++) {					// �����н��б���
						tmp.data[i][k] -= c * tmp.data[j][k];
						ret.data[i][k] -= c * tmp.data[j][k];
					}
				}
			}
		}
	}

	// ���Խ�Ԫ�ػ�Ϊ 1
	for (int i = 0; i < rows; i++) {
		c = 1 / tmp.data[i][i];
		for (int j = 0; j < rows; j++) {
			tmp.data[i][j] *= c;
			ret.data[i][j] *= c;
		}
	}

	// ת��Ϊ��λ����
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < rows; j++) {					// �Ը��������н��м���
			c = tmp.data[i][j];
			if (j != i) {									// �� [i][[i] �����м���
				for (int k = 0; k < rows; k++) {
					tmp.data[i][k] -= c * tmp.data[j][k];
					ret.data[i][k] -= c * ret.data[j][k];
				}
			}
		}
	}

	return ret;
}

void Matrix::print() {										// ��ӡ����

	cout << endl;
	for (size_t i = 0; i < (int)data.size(); i++) {
		for (size_t j = 0; j < (int)data[0].size(); j++) {
			// cout << data[i][j] << "\t";
			printf("%.6f\t", data[i][j]);
		}
		cout << endl;
	}
}