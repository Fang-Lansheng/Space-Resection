#pragma once
#include "stdafx.h"
#include <iostream>
#include <math.h>
#include <assert.h>
#include <vector>

using namespace std;

class Matrix {

/// 变量
public:
	int rows;				// 行高
	int cols;				// 列宽
	int size;				// 容器大小
	vector<vector<double>> data;
	
/// 函数
public:
	Matrix() {};											// 默认构造函数
	Matrix(int m, int n);									// 构造函数
	Matrix(int n);											// 构造函数（方阵）
	virtual ~Matrix();										// 析构函数
	double get(int row, int col) const;						// 得到矩阵在 row 行 col 列的值
	void set(int row, int col, double value);				// 将矩阵 row 行 col 列的值赋为 value

	bool operator=(const Matrix&);
	friend Matrix operator+(const Matrix&, const Matrix&);	// 重载运算符 +
	friend Matrix operator-(const Matrix&, const Matrix&);	// 重载运算符 -
	friend Matrix operator*(const Matrix&, const Matrix&);	// 重载运算符 *
	friend Matrix operator*(double,		   const Matrix&);	// 重载运算符 *
	friend Matrix operator*(const Matrix&, double);			// 重载运算符 *
	friend Matrix operator/(const Matrix&, double);			// 重载运算符 /

	double det();											// 求矩阵行列式
	Matrix t() const;										// 求矩阵转置
	Matrix inv();											// 求矩阵的逆矩阵
	void print();											// 打印矩阵
};

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

Matrix Matrix::inv() {										// 矩阵求逆
	
	double Det = det();
	assert(Det != 0.0);										// 行列式不能为零

	Matrix tmp(rows);										// 相当于增广矩阵的左边矩阵
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < rows; j++) {
			tmp.data[i][j] = data[i][j];
		}
	}
	Matrix ret(rows);										// 逆矩阵（初始化为单位矩阵）
	for (int i = 0; i < rows; i++) {
		ret.set(i, i, 1.0);
	}
	
	int r = 0;												// 行号
	double c = 1.0;											// 进行初等变换时的系数，设初值为 1

	// 初等变换（i 表示第 i 行），将左下角元素化为 0 
	for (int i = 0; i < rows; i++) {
		if (tmp.data[i][i] == 0) {							// 如果对角线上的元素为 0，此时进行行变换	
			for (int j = i; j < rows; j++) {				// 对 i 行以后的各行进行判断，找到第 i 个元素不为 0 的行，并与第 i 行进行交换
				if (j != i) {
					if (tmp.data[j][i] != 0) {
						r = j;								// 记住该行的行号
						break;								// 退出循环
					}
				}
			}

			for (int k = 0; k < rows; k++) {				// 对第 i 行和第 r 行进行调换
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

		//for (int j = 0; j < rows; j++) {					// 对其他行的所有列进行计算（将其他行的该列化为 0）
		//	if (j != i) {									// 本行不参与计算
		//		if (tmp.data[j][i] != 0) {					// 若第 j 行第 i 列为零，则跳过（省略计算步骤）
		//			c = tmp.data[j][i] / tmp.data[i][i];		// 比例系数
		//			for (int k = 0; k < rows; k++) {			// 第 k 列（对第 j 行遍历）
		//				tmp.data[j][k] -= c * tmp.data[i][k];	// 第 j 行 = 第 j 行 - c * 第 i 行
		//				ret.data[j][k] -= c * ret.data[i][k];
		//			}
		//		}	
		//	}
		//}
		
		// 试着仅仅将“左下角”化为 0 ？（省去更多计算过程，也许结果更不易出错？）
		if (i > 0) {				// j -> 列； i -> 行
			for (int j = 0; j < i; j++) {							// 对第 i 行的 0 - j 列进行遍历，分别将其值化为 0
				if (tmp.data[i][j] != 0) {
					c = tmp.data[i][j] / tmp.data[j][j];
					for (int k = 0; k < rows; k++) {					// 对整行进行遍历
						tmp.data[i][k] -= c * tmp.data[j][k];
						ret.data[i][k] -= c * tmp.data[j][k];
					}
				}
			}
		}
	}

	// 将对角元素化为 1
	for (int i = 0; i < rows; i++) {
		c = 1 / tmp.data[i][i];
		for (int j = 0; j < rows; j++) {
			tmp.data[i][j] *= c;
			ret.data[i][j] *= c;
		}
	}

	// 转化为单位矩阵
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < rows; j++) {					// 对该行其他列进行计算
			c = tmp.data[i][j];
			if (j != i) {									// 对 [i][[i] 不进行计算
				for (int k = 0; k < rows; k++) {
					tmp.data[i][k] -= c * tmp.data[j][k];
					ret.data[i][k] -= c * ret.data[j][k];
				}
			}
		}
	}

	return ret;
}

void Matrix::print() {										// 打印矩阵

	cout << endl;
	for (size_t i = 0; i < (int)data.size(); i++) {
		for (size_t j = 0; j < (int)data[0].size(); j++) {
			// cout << data[i][j] << "\t";
			printf("%.6f\t", data[i][j]);
		}
		cout << endl;
	}
}