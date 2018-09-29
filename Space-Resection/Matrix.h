#pragma once
#include "stdafx.h"
#include <iostream>
#include <math.h>
#include <assert.h>
#include <vector>

using namespace std;

class Matrix {
public:
	// double *p = NULL;		// 指针
	int rows;				// 行高
	int cols;				// 列宽
	int size;				// 容器大小
	vector<vector<double>> data;
	
public:
	Matrix() {};											// 默认构造函数
	// Matrix(int rows, int cols, double value);			// 构造函数
	Matrix(int m, int n);									// 构造函数
	Matrix(int n);											// 构造函数（方阵）
	virtual ~Matrix();										// 析构函数
	double get(int row, int col) const;						// 得到矩阵在 row 行 col 列的值
	void set(int row, int col, double value);				// 将矩阵 row 行 col 列的值赋为 value
	// double at(int row, int col);							// 访问矩阵元素

	bool operator=(const Matrix&);
	friend Matrix operator+(const Matrix&, const Matrix&);	// 重载运算符 +
	friend Matrix operator-(const Matrix&, const Matrix&);	// 重载运算符 -
	friend Matrix operator*(const Matrix&, const Matrix&);	// 重载运算符 *
	friend Matrix operator*(double,		   const Matrix&);
	friend Matrix operator*(const Matrix&, double);
	friend Matrix operator/(const Matrix&, double);			// 重载运算符 /

	double det();											// 求矩阵行列式
	Matrix t() const;										// 求矩阵转置
	Matrix adj();											// 求矩阵的伴随矩阵
	Matrix SwitchRow(int, int);
	Matrix SwitchCol(int, int);
	Matrix inv();											// 求矩阵的逆矩阵
	void print();											// 输出矩阵
};

//Matrix::Matrix(int rows, int cols, double value) {			// 构造函数	
//
//	assert(rows > 0);
//	assert(cols > 0);
//
//	this->rows = rows;
//	this->cols = cols;
//	this->size = rows * cols;
//
//	if (size > 0 && size <= 100) {
//		p = new double[size + 1 ];
//		
//		if (p != NULL) {
//			for (int i = 0; i < size; i++) {
//				p[i] = value;
//			}
//		}
//	}
//}

Matrix::Matrix(int m, int n) {								// 构造函数

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

Matrix::Matrix(int n) {										// 构造函数（方阵）

	assert(n > 0);

	rows = n;
	cols = n;
	size = n * n;

	vector<double> tempRow(n, 0);
	for (size_t i = 0; i < n; i++) {
		data.push_back(tempRow);
	}
}

Matrix::~Matrix() {											// 析构函数

	//delete[] p;
	//p = NULL;
}

double Matrix::get(int row, int col) const {				// get 获取矩阵元素值

	assert(row >= 0 && row < this->rows);
	assert(col >= 0 && col < this->cols);

	return data[row][col];
}

void Matrix::set(int row, int col, double value) {			// set 设置矩阵元素值

	assert(row >= 0 && row < this->rows);
	assert(col >= 0 && col < this->cols);

	data[row][col] = value;
}

//double Matrix::at(int row, int col) {						// at 返回指针
//	return p[this->rows * row + col];
//}

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

Matrix operator+(const Matrix& m1, const Matrix& m2) {		// 矩阵相加

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

Matrix operator-(const Matrix& m1, const Matrix& m2) {		// 矩阵相减

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

Matrix operator*(const Matrix& m1, const Matrix& m2) {		// 矩阵相乘
	
	assert(m1.size > 0);
	assert(m2.size > 0);
	assert(m1.cols == m2.rows);

	Matrix ret(m1.rows, m2.cols);

	for (int i = 0; i < m1.rows; i++) {
		for (int j = 0; j < m2.cols; j++) {
			for (int k = 0; k < m1.cols; k++) {	// m1.cols == m2.rows
				// ret.p[i * m1.rows + j] = m1.p[i * m1.rows + k] * m2.p[k * m2.cols + j];
				ret.data[i][j] = m1.data[i][k] + m2.data[k][j];
			}
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

Matrix operator*(double value, const Matrix& m1) {			// 常数 × 矩阵

	Matrix ret(m1.rows, m1.cols);

	for (int i = 0; i < m1.size; i++) {
		// ret.p[i] = m1.p[i] * value;
		ret.data[i] = value * m1.data[i];
	}

	return ret;
}

Matrix operator*(const Matrix& m1, double value) {			// 矩阵 × 常数

	Matrix ret(m1.rows, m1.cols);

	for (int i = 0; i < m1.size; i++) {
		// ret.p[i] = m1.p[i] * value;
		ret.data[i] = value * m1.data[i];
	}

	return ret;
}

Matrix operator/(const Matrix& m1, double value) {			// 矩阵 ÷ 常数

	Matrix ret(m1.rows, m1.cols);

	for (int i = 0; i < m1.size; i++) {
		// ret.p[i] = m1.p[i] / value;
		ret.data[i] = m1.data[i] / value;
	}

	return ret;
}

double dets(int n, vector<vector<double>> aa) {					// 求方阵的行列式

	if (n == 1)		
		return aa[0][0];	

	vector<vector<double>> bb((n - 1), vector<double>(n - 1));	// 创建 n - 1 阶的代数余子式阵 bb  
	int move = 0;												// 判断行是否移动       	
	double sum = 0.0;											// sum 为行列式的值      	

	for (int a_row = 0; a_row < n; a_row++) {					// a 的行数把矩阵 a(nn) 赋值到 b(n - 1)      	
		for (int b_row = 0; b_row < n - 1; b_row++) {			// 把 aa 阵第一列各元素的代数余子式存到bb      		
			move = a_row > b_row ? 0 : 1;						// bb 中小于 arow 的行，同行赋值，等于的错过，大于的加一      			
			for (int j = 0; j < n - 1; j++) {					// 从 aa 的第二列赋值到第n列      									
				bb[b_row][j] = aa[b_row + move][j + 1];
			}		
		}		
		int flag = (a_row % 2 == 0 ? 1 : -1);					// 因为列数为 0，所以行数是偶数时，代数余子式为 1      		
		sum += flag * dets(n - 1, bb) * aa[a_row][0];			// aa 第一列各元素与其代数余子式积的和即为行列式    	
	}	

	
	return sum;
}

double Matrix::det() {										// 求矩阵行列式
	assert(rows == cols);
	return dets(rows, data);
}

Matrix Matrix::t() const {									// 矩阵转置

	Matrix ret(cols, rows);
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			ret.data[j][i] = data[i][j];
		}
	}

	return ret;
}

Matrix Matrix::adj() {										// 求矩阵的伴随矩阵
	
	assert(rows == cols);

	const int n = rows;		// 矩阵为方阵
	Matrix ret(n);			// 伴随矩阵
	vector<vector<double>> bb((n - 1), vector<double>(n - 1));	// 创建 n - 1 阶的代数余子式

	int pi, pj, flag = 0;
	for (int ai = 0; ai < n; ai++) {
		for (int aj = 0; aj < n; aj++) {
			for (int bi = 0; bi < n - 1; bi++) {
				for (int bj = 0; bj < n - 1; bj++) {		// ai 行的代数余子式：
					if (bi < ai)	
						pi = 0;								// 小于 ai 的行，aa 与 bb 同行赋值
					else			
						pi = 1;								// 大于等于 ai 的行，取 aa 阵的 ai + 1 行赋值给 bb 阵的第 bi 行
					if (bj < aj)	
						pj = 0;								// 小于 aj 的列，aa 与 bb 同列赋值
					else			
						pj = 1;								// 大于等于 aj 的列，取 aa 阵的 aj + 1 列赋值给 bb 阵的第 bj 列
	
					bb[bi][bj] = data[bi + pi][bj + pj];				
				}
			}
			flag = ((ai + aj) % 2 == 0) ? 1 : -1;
			ret.data[ai][aj] = flag * dets(n - 1, bb);
		}
	}

	return ret;
}


Matrix Matrix::SwitchRow(int row1, int row2) {
	
	assert(row1 >= 0 && row1 < rows);
	assert(row2 >= 0 && row2 < rows);

	Matrix ret(rows, cols);
	ret.data = data;

	for (int i = 0; i < cols; i++) {
		ret.data[row2][i] = data[row1][i];
		ret.data[row1][i] = data[row2][i];
	}
	return ret;
}

Matrix Matrix::inv() {										// 矩阵求逆
	
	double Det = det();
	assert(Det != 0.0);

	Matrix Adj = adj();										// 伴随矩阵
	Matrix ret(rows);										// 逆矩阵

	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			// ret.p[i * rows + j] = Adj.p[i * rows + j] / Det;
			ret.data[i][j] = Adj.data[i][j] / Det;
		}
	}

	return ret;
}

void Matrix::print() {										// 打印矩阵

	cout << endl;
	for (size_t i = 0; i < (int)data.size(); i++)
	{
		for (size_t j = 0; j < (int)data[0].size(); j++)
		{
			cout << data[i][j] << "\t";
		}
		cout << endl;
	}
}

