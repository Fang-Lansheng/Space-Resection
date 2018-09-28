// SpaceResection.cpp: 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "Matrix.h"
#include <iostream>
#include <fstream>
#include <math.h>

#define N 4					// 控制点数
#define PI 3.1415926
#define MAXITERATIONS 8		// 允许的最大迭代次数
#define PRECISION 3.0e-5	// 精度（角元素的改正数应 < 0.1' = 3.0e-5)

using namespace std;

struct EOEO {				// 外方位元素 elements of exterior orientation
	double Xs;				// 外方位元素线元素
	double Ys;
	double Zs;
	double phi;				// 外方位元素角元素	
	double omega;
	double kappa;
};

struct SrcData {			// 已知数据，共四对已知控制点的影像坐标及地面坐标
	double x;				// 影像 x 坐标
	double y;				// 影像 y 坐标
	double X;				// 地面 X 坐标
	double Y;				// 地面 Y 坐标
	double Z;				// 地面 Z 坐标
};

// 函数声明
void InitInterface();
bool CheckPrecison(Matrix &Correction);
void InitData(SrcData* sd, char* filename);
void Iterator(SrcData sd[N], double m, double f);
void OutputResult(EOEO* eoeo, Matrix &R, double* M, double m0);

/// --------------------------------------------------------------------------------
///	@ Function ---> 初始化坐标数据
///	@ Parameters -> sd：保存已知数据的结构体数组；filename：保存数据的文件名
///	@ Return -----> void
/// --------------------------------------------------------------------------------
void InitData(SrcData* sd, char* filename) {

	fstream datafile(filename, ios::in | ios::out);
	if (!datafile) {							// 检查文件是否存在
		cout << "Error: 打开文件失败，请检查文件是否存在！" << endl;
		exit(1);								// 非正常退出
	}

	memset(sd, 0, N * sizeof(SrcData));			// 为数组sd分配内存空间

	for (int i = 0; i < N; i++) {
		datafile >> sd[i].x >> sd[i].y >> sd[i].X >> sd[i].Y >> sd[i].Z;
	}
	datafile.close();							// 关闭文件

	cout << "\n> 数据读取完毕..." << endl;
	cout << "\t影像坐标(mm)\t\t\t" << "地面坐标(m)\t" << endl;
	for (int i = 0; i < N; i++) {
		cout << i << "\t" << sd[i].x << "\t" << sd[i].y << "\t\t" << sd[i].X << "\t\t" << sd[i].Y << "\t\t" << sd[i].Z << endl;
	}
}


/// --------------------------------------------------------------------------------
///	@ Function ---> 检查改正数是否已经达到了精度要求
///	@ Parameters -> corrections：保存改正数的数组
///	@ Return -----> 布尔值
/// --------------------------------------------------------------------------------
bool CheckPrecison(Matrix &Correction) {

	bool Boolean;
	Boolean = {
		fabs(Correction.at(0, 0)) < 0.000001 &&
		fabs(Correction.at(1, 0)) < 0.000001 &&
		fabs(Correction.at(2, 0)) < 0.000001 &&
		fabs(Correction.at(3, 0)) < PRECISION &&
		fabs(Correction.at(4, 0)) < PRECISION &&
		fabs(Correction.at(5, 0)) < PRECISION };
	return Boolean;
}

/// --------------------------------------------------------------------------------
///	@ Function ---> 迭代器，计算的主体部分
///	@ Parameters -> sd：保存原始数据的结构体数组；m：摄影比例尺的分母；f：摄影机主距
///	@ Return -----> void
/// --------------------------------------------------------------------------------
void Iterator(SrcData sd[N], double m, double f) {

	double phi, omega, kappa, Xs, Ys, Zs;	// 外方位元素
	Xs = Ys = Zs = 0.0;

	for (int i = 0; i < N; i++) {
		sd[i].x /= 1000;		// 单位换算 mm -> m
		sd[i].y /= 1000;
		Xs += sd[i].X;
		Ys += sd[i].Y;
	}

	// 确定未知数的初始值
	Xs /= N;
	Ys /= N;
	Zs = m * f / 1000.0;		// 主距 f 的单位为 mm，故除 1000
	phi = omega = kappa = 0.0;	// κ 的初始值也取 0

	cout << "\n> 外方位元素的初始值为：" << endl;
	cout << "Xs: " << Xs << "\t" << "Ys: " << Ys << "\t" << "Zs: " << Zs << "\t" << endl;
	cout << "φ: " << phi << "\t\t" << "ω: " << omega << "\t\t" << "κ: " << kappa << "\t\t" << endl;

	double x0(0), y0(0);	// 内方位元素，假设为 0
	double X0[N] = { 0.0 };
	double Y0[N] = { 0.0 };
	double Z0[N] = { 0.0 };

	Matrix R(3, 3);
	Matrix V(8, 1);
	Matrix A(2, 6);
	Matrix L(2, 1);
	Matrix Correction(6, 1);
	Matrix Buf_A(8, 6);
	Matrix Buf_L(8, 1);
	Matrix ATA(6, 6);
	Matrix ATL(6, 1);


	int iCount = 0;					// 迭代次数
	cout << "\n> 开始迭代运算..." << endl;
	while (!CheckPrecison(Correction)) {
		cout << "*************************************************" << endl;
		cout << endl << ">> 第 " << ++iCount << " 次迭代：" << endl;
		if (iCount == MAXITERATIONS) {
			cout << "迭代次数超限，可能不收敛" << endl;
			break; 
		}

		// 每次迭代前必须清空两个保存累加值的矩阵 ATA 和 ATL
		for (int i = 0; i < 6; i++) {
			ATL.set(i, 0, 0.0);
			for (int j = 0; j < 6; j++) {
				ATA.set(i, j, 0.0);
			}
		}


		// 计算旋转矩阵（9 个方向余弦）
		R.p[0 * R.rows + 0] = cos(phi) * cos(kappa) - sin(phi) * sin(omega) * sin(kappa);	// a1
		R.p[0 * R.rows + 1] = -cos(phi) * sin(kappa) - sin(phi) * sin(omega) * cos(kappa);	// a2
		R.p[0 * R.rows + 2] = -sin(phi) * cos(omega);										// a3
		R.p[1 * R.rows + 0] = cos(omega) * sin(kappa);										// b1
		R.p[1 * R.rows + 1] = cos(omega) * cos(kappa);										// b2
		R.p[1 * R.rows + 2] = -sin(omega);													// b3
		R.p[2 * R.rows + 0] = sin(phi) * cos(kappa) + cos(phi) * sin(omega) * sin(kappa);	// c1
		R.p[2 * R.rows + 1] = -sin(phi) * sin(kappa) + cos(phi) * sin(omega) * cos(kappa);	// c2
		R.p[2 * R.rows + 2] = cos(phi) * cos(omega);										// c3

		for (int i = 0; i < N; i++) {
			Z0[i] = R.p[0 * R.rows + 2] * (sd[i].X - Xs) + R.p[1 * R.rows + 2] * (sd[i].Y - Ys) + R.p[2 * R.rows + 2] * (sd[i].Z - Zs);	// Z0 = Z_bar = a3(X-Xs)+b3(Y-Ys)+c3(Z-Zs)
			X0[i] = x0 - f * (R.p[0 * R.rows + 0] * (sd[i].X - Xs) + R.p[1 * R.rows + 0] * (sd[i].Y - Ys) + R.p[2 * R.rows + 0] * (sd[i].Z - Zs)) / Z0[i];		// X0 为观测值
			Y0[i] = y0 - f * (R.p[0 * R.rows + 1] * (sd[i].X - Xs) + R.p[1 * R.rows + 1] * (sd[i].Y - Ys) + R.p[2 * R.rows + 1] * (sd[i].Z - Zs)) / Z0[i];		// Y0 为观测值

																																		// 计算误差方程式的系数
			A.p[0 * A.rows + 0] = (R.p[0 * R.rows + 0] * f + R.p[0 * R.rows + 2] * (sd[i].x - x0)) / Z0[i];								// a11
			A.p[0 * A.rows + 1] = (R.p[1 * R.rows + 0] * f + R.p[1 * R.rows + 2] * (sd[i].x - x0)) / Z0[i];								// a12
			A.p[0 * A.rows + 2] = (R.p[2 * R.rows + 0] * f + R.p[2 * R.rows + 2] * (sd[i].x - x0)) / Z0[i];								// a13
			A.p[0 * A.rows + 3] = (sd[i].y - y0) * sin(omega) - ((sd[i].x - x0) * ((sd[i].x - x0) * cos(kappa) - (sd[i].y - y0) * sin(kappa)) / f + f * cos(kappa)) * cos(omega);	// a14
			A.p[0 * A.rows + 4] = -f * sin(kappa) - (sd[i].x - x0) * ((sd[i].x - x0) * sin(kappa) + (sd[i].y - y0) * cos(kappa)) / f;	// a15
			A.p[0 * A.rows + 5] = sd[i].y - y0;																							// a16
			A.p[1 * A.rows + 0] = (R.p[0 * R.rows + 1] * f + R.p[0 * R.rows + 2] * (sd[i].y - y0)) / Z0[i];								// a21
			A.p[1 * A.rows + 1] = (R.p[1 * R.rows + 1] * f + R.p[0 * R.rows + 2] * (sd[i].y - y0)) / Z0[i];								// a22
			A.p[1 * A.rows + 2] = (R.p[2 * R.rows + 1] * f + R.p[0 * R.rows + 2] * (sd[i].y - y0)) / Z0[i];								// a23
			A.p[1 * A.rows + 3] = -(sd[i].x - x0) * sin(omega) - ((sd[i].y - y0) * ((sd[i].x - x0) * cos(kappa) - (sd[i].y - y0) * sin(kappa)) / f - f * sin(kappa)) * cos(omega);	// a24
			A.p[1 * A.rows + 4] = -f * cos(kappa) - (sd[i].y - y0) * ((sd[i].x - x0) * cos(kappa) - (sd[i].y - y0) * sin(kappa)) / f;	// a25
			A.p[1 * A.rows + 5] = -(sd[i].x - x0);																						// a26

			// 该循环保存A矩阵，最后评定精度用
			for (int j = 0; j < 2; j++) {
				for (int k = 0; k < 6; k++) {
					Buf_A.p[(2 * i + j) * Buf_A.rows + k] = Buf_A.p[j * A.rows + k];
				}
			}

			L.p[0 * L.rows + 0] = sd[i].x - X0[i];
			L.p[1 * L.rows + 0] = sd[i].y - Y0[i];


			// 保存 L 矩阵，最后评定精度用
			for (int j = 0; i < 2; i++) {
				Buf_L.p[(2 * i + j) * Buf_L.rows + 0] = L.p[j * L.rows + 0];
				// Buf_L.at<double>(2 * i + j, 0) = L.at<double>(j, 0);
			}

			// 所谓的逐步法化，即要在循环内部就将 ATA 计算出来累加，下面的 L 矩阵类似

			Matrix AT(6, 2);
			AT = A.t();

			// ATA = AT * A;
			for (int l = 0; l < AT.rows; l++) {
				for (int m = 0; m < A.cols; m++) {
					for (int n = 0; n < AT.cols; n++) {
						ATL.p[l * AT.rows + m] = AT.p[l * AT.rows + n] * A.p[n * A.cols + m];
					}
				}
			}
			// ATL = AT * L;
			for (int l = 0; l < AT.rows; l++) {
				for (int m = 0; m < L.cols; m++) {
					for (int n = 0; n < AT.cols; n++) {
						ATL.p[l * AT.rows + m] = AT.p[l * AT.rows + n] * L.p[n * L.cols + m];
					}
				}
			}
		}// for

		 // “逐步法化”的另一处不同，出循环即可直接计算 ATA 的逆乘 ATL
		Correction = ATA.inv() * ATL;

		// correction 即为改正数，每一次迭代时用未知数近似值与上次迭代计算的改正数之和作为新的近似值
		Xs		+= Correction.p[0 * Correction.rows + 0];
		Ys		+= Correction.p[1 * Correction.rows + 0];
		Zs		+= Correction.p[2 * Correction.rows + 0];
		phi		+= Correction.p[3 * Correction.rows + 0];
		omega	+= Correction.p[4 * Correction.rows + 0];
		kappa	+= Correction.p[5 * Correction.rows + 0];


		cout << "【改正数值】：" << endl;
		cout << "\tΔXs = " << Correction.p[0 * Correction.rows + 0] << "\tΔYs = " << Correction.p[1 * Correction.rows + 0] << "\tΔZs = " << Correction.p[2 * Correction.rows + 0] << endl;
		cout << "\tΔφ = " << Correction.p[3 * Correction.rows + 0] << "\tΔω = " << Correction.p[4 * Correction.rows + 0] << "\tΔκ = " << Correction.p[5 * Correction.rows + 0] << endl;

		cout << "\n【外方位元素值】：" << endl;
		cout << "\tXs = " << Xs << "\tYs = " << Ys << "\tZs = " << Zs << endl;
		cout << "\tφ = " << phi << "\tω = " << omega << "\tκ = " << kappa << endl;
	}// while

	EOEO eoeo;		// 外方位元素
	eoeo.Xs = Xs;
	eoeo.Ys = Ys;
	eoeo.Zs = Zs;
	eoeo.phi = phi;
	eoeo.omega = omega;
	eoeo.kappa = kappa;

	cout << ">>>正常退出迭代<<<" << endl << endl;

	// 精度评定
	double Q[6] = { 0.0 };
	for (int i = 0; i < 6; i++) {
		// Q[i] = ATA.at<double>(i, i);
		Q[i] = ATL.p[i * ATL.rows + i];
	}

	double m0(0);	// 单位权中误差
	double vv(0);	// [vv]，即平方和

	V = Buf_A * Correction;
	V = V - Buf_L;

	for (int i = 0; i < 8; i++) {
		// vv += pow(V.at<double>(i, 0), 2);
		vv += pow(V.p[i * V.rows + 0], 2);
	}
	m0 = sqrt(vv / (2 * N - 6));	// 中误差 m0

	double M[6] = { 0.0 };			// 保存六个值的中误差
	for (int i = 0; i < 6; i++) {
		M[i] = m0 * sqrt(Q[i]);
		if (i > 2) {
			M[i] = M[i] * 180 * 3600 / PI;	// 转换为角度制
		}
	}

	OutputResult(&eoeo, R, M, m0);		// 输出结果
	cout << endl << ">>>解算全部完成<<<" << endl << endl;
}

/// --------------------------------------------------------------------------------
///	@ Function ---> 初始化程序界面
///	@ Parameters -> None
///	@ Return -----> void
/// --------------------------------------------------------------------------------
void InitInterface() {
	cout << "单像空间后方交会程序（C++）" << endl;
	cout << "遥感信息工程学院 张济帆 2016302590060" << endl;
}

/// --------------------------------------------------------------------------------
///	@ Function ---> 输出解算结果
///	@ Parameters -> eoeo：指向最终解算结果的结构体数组；
///					R：旋转矩阵；
///					M：保存计算精度的数组；
///					m0：单位权中误差
///	@ Return -----> void
/// --------------------------------------------------------------------------------
void OutputResult(EOEO* eoeo, Matrix &R, double* M, double m0) {
	cout << "/////---------- 计算结果 ----------/////" << endl;
	cout << "六个外方位元素为：" << endl;
	printf("Xs = %.4f\n", eoeo->Xs);
	printf("Ys = %.4f\n", eoeo->Ys);
	printf("Zs = %.4f\n", eoeo->Zs);
	printf("φ = %.10f\n", eoeo->phi);
	printf("ω = %.10f\n", eoeo->omega);
	printf("κ = %.10f\n", eoeo->kappa);

	cout << endl << "旋转矩阵：" << endl;
	cout << "\t" << R.p[0 * R.rows + 0] << "\t" << R.p[0 * R.rows + 1] << "\t" << R.p[0 * R.rows + 2] << endl;
	cout << "\t" << R.p[1 * R.rows + 0] << "\t" << R.p[1 * R.rows + 1] << "\t" << R.p[1 * R.rows + 2] << endl;
	cout << "\t" << R.p[2 * R.rows + 0] << "\t" << R.p[2 * R.rows + 1] << "\t" << R.p[2 * R.rows + 2] << endl;

	cout << "单位权中误差：" << m0 << endl << endl;

	cout << "六个外方位元素的精度：" << endl;
	cout << "Xs:	" << M[0] << "m" << endl;
	cout << "Ys:	" << M[1] << "m" << endl;
	cout << "Zs:    " << M[2] << "m" << endl;
	cout << "φ:    " << M[3] << "'" << endl;
	cout << "ω:    " << M[4] << "'" << endl;
	cout << "κ:    " << M[5] << "'" << endl;
}



/////////////////////////////////////////////////////////////////////////////////////////
// 主函数
int main() {
	//InitInterface();

	//SrcData sd[N];
	//double m = 15000.0;		// 摄影比例尺 1/m = 1/15000
	//double f = 153.24;		// 主距 f = 153.24mm

	//char filename[50] = "data.txt";	// 储存已知坐标数据的文件名为 data.txt
	//InitData(sd, filename);
	//Iterator(sd, m, f);


	Matrix A(2, 2);
	A.set(1, 1, 2.0);
	Matrix B(2, 2);
	B.set(1, 1, 3.0);
	double val = B.get(1, 1);
	A = A + B;
	Matrix O(2, 2);
	O.set(0, 0, 2.0);
	O.p[0 * O.rows + 1] = 1.0;
	O.p[1 * O.rows + 0] = 1.0;
	O.p[1 * O.rows + 1] = 2.0;
	double det = O.det();
	cout << "\tO:" << endl;
	cout << "\t" << O.p[0 * O.rows + 0] << "\t" << O.p[0 * O.rows + 1] << endl;
	cout << "\t" << O.p[1 * O.rows + 0] << "\t" << O.p[1 * O.rows + 1] << endl;
	cout << "\tdet(O) = " << det << endl;
	Matrix Inv = O.inv();
	cout << "\tInv:" << endl;
	cout << "\t" << Inv.p[0 * Inv.rows + 0] << "\t" << Inv.p[0 * Inv.rows + 1] << endl;
	cout << "\t" << Inv.p[1 * Inv.rows + 0] << "\t" << Inv.p[1 * Inv.rows + 1] << endl;
	cout << "~~~~~~~" << A.at(1, 1) << endl;
	return 0;
}

