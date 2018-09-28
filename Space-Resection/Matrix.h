#pragma once
#include "stdafx.h"
#include <iostream>
#include <math.h>
#include <assert.h>

using namespace std;

class Matrix {
public:
	double *p;		// ָ��
	int rows;		// �и�
	int cols;		// �п�
	int size;		// ������С

public:
	Matrix();												// Ĭ�Ϲ��캯��
	Matrix(int rows, int cols, double value = 0.0);			// ���캯��
	Matrix(int n);											// ���캯��������
	virtual ~Matrix();										// ��������
	double get(int row, int col) const;						// �õ������� row �� col �е�ֵ
	void set(int row, int col, double value);				// ������ row �� col �е�ֵ��Ϊ value
	double at(int row, int col);							// ���ʾ���Ԫ��

	bool operator=(const Matrix&);
	friend Matrix operator+(const Matrix&, const Matrix&);	// ��������� +
	friend Matrix operator-(const Matrix&, const Matrix&);	// ��������� -
	friend Matrix operator*(const Matrix&, const Matrix&);	// ��������� *
	friend Matrix operator*(double,		   const Matrix&);
	friend Matrix operator*(const Matrix&, double);
	friend Matrix operator/(const Matrix&, double);			// ��������� /

	double det();											// ���������ʽ
	Matrix t() const;										// �����ת��
	Matrix adj();											// �����İ������
	Matrix inv();											// �����������
};

Matrix::Matrix(int rows, int cols, double value) {			// ���캯��	

	assert(rows > 0);
	assert(cols > 0);

	this->rows = rows;
	this->cols = cols;
	this->size = rows * cols;

	p = new double[rows * cols];
	
	for (int i = 0; i < rows * cols; i++) {
		p[i] = value;
	}
}

Matrix::Matrix(int n) {						// ���캯��������
	assert(n > 0);

	this->rows = n;
	this->cols = n;
	this->size = n * n;

	p = new double[n * n];

	for (int i = 0; i < n * n; i++) {
		p[i] = 0.0;
	}
}

Matrix::~Matrix() {											// ��������

	// delete[] p;
	// p = NULL;
}

double Matrix::get(int row, int col) const {				// get ��ȡ����Ԫ��ֵ

	assert(row >= 0 && row < this->rows);
	assert(col >= 0 && col < this->cols);

	return p[this->rows * row + col];
}

void Matrix::set(int row, int col, double value) {			// set ���þ���Ԫ��ֵ

	assert(row >= 0 && row < this->rows);
	assert(col >= 0 && col < this->cols);

	p[this->rows * row + col] = value;
}

double Matrix::at(int row, int col) {						// at ����ָ��
	return p[this->rows * row + col];
}

bool Matrix::operator=(const Matrix& m) {

	assert(this->rows == m.rows);
	assert(this->cols == m.cols);

	for (int i = 0; i < this->rows; i++) {
		for (int j = 0; j < this->cols; j++) {
			p[this->rows * i + j] = m.p[this->rows * i + j];
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
	
	assert(m1.cols == m2.rows);

	Matrix ret(m1.rows, m1.cols);

	for (int i = 0; i < m1.rows; i++) {
		for (int j = 0; j < m2.cols; j++) {
			for (int k = 0; k < m1.cols; k++) {	// m1.cols == m2.rows
				ret.p[i * m1.rows + j] = m1.p[i * m1.rows + k] * m2.p[k * m2.cols + j];
			}
		}
	}
	return ret;
}

Matrix operator*(double value, const Matrix& m1) {			// ���� �� ����

	Matrix ret(m1.rows, m1.cols);

	for (int i = 0; i < m1.size; i++) {
		ret.p[i] = m1.p[i] * value;
	}

	return ret;
}

Matrix operator*(const Matrix& m1, double value) {			// ���� �� ����

	Matrix ret(m1.rows, m1.cols);

	for (int i = 0; i < m1.size; i++) {
		ret.p[i] = m1.p[i] * value;
	}

	return ret;
}

Matrix operator/(const Matrix& m1, double value) {			// ���� �� ����

	Matrix ret(m1.rows, m1.cols);

	for (int i = 0; i < m1.size; i++) {
		ret.p[i] = m1.p[i] / value;
	}

	return ret;
}

double dets(int n, double *&aa) {							// ���������ʽ

	if (n == 1)		
		return aa[0];	

	double *bb = new double[(n - 1)*(n - 1)];				// ���� n - 1 �׵Ĵ�������ʽ�� bb        	
	int move = 0;											// �ж����Ƿ��ƶ�       	
	double sum = 0.0;										// sum Ϊ����ʽ��ֵ      	

	for (int a_row = 0; a_row < n; a_row++) {				// a �������Ѿ��� a(nn) ��ֵ�� b(n - 1)      	
		for (int b_row = 0; b_row < n - 1; b_row++) {		// �� aa ���һ�и�Ԫ�صĴ�������ʽ�浽bb      		
			move = a_row > b_row ? 0 : 1;					// bb ��С�� arow ���У�ͬ�и�ֵ�����ڵĴ�������ڵļ�һ      			
			for (int j = 0; j < n - 1; j++) {				// �� aa �ĵڶ��и�ֵ����n��      							
				bb[b_row * (n - 1) + j] = aa[(b_row + move) * n + j + 1];			
			}		
		}		
		int flag = (a_row % 2 == 0 ? 1 : -1);				// ��Ϊ����Ϊ 0������������ż��ʱ����������ʽΪ 1      		
		sum += flag * aa[a_row * n] * dets(n - 1, bb);		// aa ��һ�и�Ԫ�������������ʽ���ĺͼ�Ϊ����ʽ    	
	}	

	delete[] bb;	
	return sum;
}

double Matrix::det() {										// ���������ʽ
	assert(rows == cols);
	return dets(rows, p);
}

Matrix Matrix::t() const {									// ����ת��
	Matrix ret(cols, rows);

	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			ret.p[j * cols + i] = p[i * rows + j];
		}
	}
	return ret;
}

Matrix Matrix::adj() {										// �����İ������
	
	assert(rows == cols);

	const int n = rows;		// ����Ϊ����
	Matrix ret(n);			// �������
	double* bb = new double[(n - 1) * (n - 1)];				// ���� n - 1 �׵Ĵ�������ʽ

	int pi, pj, flag = 0;
	for (int ai = 0; ai < n; ai++) {
		for (int aj = 0; aj < n; aj++) {
			for (int bi = 0; bi < n - 1; bi++) {
				for (int bj = 0; bj < n - 1; bj++) {		// ai �еĴ�������ʽ��
					if (bi < ai)	
						pi = 0;								// С�� ai ���У�aa �� bb ͬ�и�ֵ
					else			
						pi = 1;								// ���ڵ��� ai ���У�ȡ aa ��� ai + 1 �и�ֵ�� bb ��ĵ� bi ��
					if (bj < aj)	
						pj = 0;								// С�� aj ���У�aa �� bb ͬ�и�ֵ
					else			
						pj = 1;								// ���ڵ��� aj ���У�ȡ aa ��� aj + 1 �и�ֵ�� bb ��ĵ� bj ��
					
					bb[bi * (n - 1) + bj] = p[(bi + pi) * n + bj + pj];
					
				}
			}

			flag = ((ai + aj) % 2 == 0) ? 1 : -1;
			ret.p[ai * n + aj] = flag * dets(n - 1, bb);
		}
	}

	delete[] bb;
	return ret;
}

Matrix Matrix::inv() {										// ��������
	
	double Det = det();
	assert(Det != 0.0);

	Matrix Adj = adj();										// �������
	Matrix ret(rows);										// �����

	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			ret.p[i * rows + j] = Adj.p[i * rows + j] / Det;
		}
	}

	return ret;
}
