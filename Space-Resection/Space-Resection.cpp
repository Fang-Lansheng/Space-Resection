// SpaceResection.cpp: 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <cv.h>				// 利用 OpenCV 进行矩阵运算	
#include <highgui.h>

#define N 4					// 控制点数
#define PI 3.1415926
#define MAXITERATIONS 8		// 允许的最大迭代次数
#define PRECISION 3.0e-5	// 精度（角元素的改正数应 < 0.1' = 3.0e-5)

using namespace cv;
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
bool CheckPrecison(Mat &correction);
void InitData(SrcData* sd, char* filename);
void Iterator(SrcData sd[N], double m, double f);
void OutputResult(EOEO* eoeo, Mat &R, double* M, double m0);

/// --------------------------------------------------------------------------------
///	@ Function ---> 初始化坐标数据
///	@ Parameters -> sd：保存已知数据的结构体数组；filename：保存数据的文件名
///	@ Return -----> void
/// --------------------------------------------------------------------------------
void InitData(SrcData* sd, char* filename) {

	//cout << "开始读取坐标数据" << endl;			// 读取坐标数据
	//FILE *fp;
	//if ((fp = fopen(filename, "r")) == NULL) {
	//	cout << "Error: 打开文件失败，请检查文件是否存在！" << endl;
	//	exit(1);
	//}
	//while (!feof(fp)) {
	//	for (int i = 0; i < 4; i++) {
	//		fscanf_s(fp, "%lf %lf %lf %lf %lf", &sd[i].x, &sd[i].y, &sd[i].X, &sd[i].Y, &sd[i].Z);
	//	}
	//}
	//fclose(fp);

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
bool CheckPrecison(Mat &correction) {

	bool Boolean;
	Boolean = {
		fabs(correction.at<double>(0, 0)) < 0.000001 &&
		fabs(correction.at<double>(1, 0)) < 0.000001 &&
		fabs(correction.at<double>(2, 0)) < 0.000001 &&
		fabs(correction.at<double>(3, 0)) < PRECISION &&
		fabs(correction.at<double>(4, 0)) < PRECISION &&
		fabs(correction.at<double>(5, 0)) < PRECISION };
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

	Mat R = Mat::zeros(Size(3, 3), CV_64FC1);	// 旋转矩阵 R

	Mat V = Mat::zeros(Size(8, 1), CV_64FC1);
	Mat A = Mat::zeros(Size(2, 6), CV_64FC1);			// 系数矩阵，即误差方程式 V = AX - L 中的 A 矩阵
	Mat correction = Mat::ones(Size(6, 1), CV_64FC1);	// 改正数矩阵，即误差方程式 V = AX - L 中的 X 矩阵
	Mat L = Mat::zeros(Size(2, 1), CV_64FC1);

	Mat Buf_A = Mat::zeros(Size(8, 6), CV_64FC1);	// 存储 8 × 6 的 A 矩阵
	Mat Buf_L = Mat::zeros(Size(8, 1), CV_64FC1);	// 存储 8 × 1 的 L 矩阵

	Mat ATA = Mat::zeros(Size(6, 6), CV_64FC1);
	Mat ATL = Mat::zeros(Size(6, 1), CV_64FC1);
	uchar* pATA = ATA.data;
	uchar* pATL = ATL.data;


	int iCount = 0;					// 迭代次数
	cout << "\n> 开始迭代运算..." << endl;
	while (!CheckPrecison(correction)) {
		cout << "*************************************************" << endl;
		cout << endl << ">> 第 " << ++iCount << " 次迭代：" << endl;
		if (iCount == MAXITERATIONS) {
			cout << "迭代次数超限，可能不收敛" << endl;
			break; 
		}

		// 每次迭代前必须清空两个保存累加值的矩阵 ATA 和 ATL
		for (int i = 0; i < 6; i++) {
			pATL[i * 6 + 0] = 0;
			for (int j = 0; j < 6; j++) {
				pATA[i * 6 + j] = 0;
			}
		}
		//for (int i = 0; i < 6; i++) {
		//	ATL.at<double>(i, 0) = 0;
		//	for (int j = 0; j < 6; j++) {
		//		ATA.at<double>(i, j) = 0;
		//	}
		//}


		// 计算旋转矩阵（9 个方向余弦）
		R.at<double>(0, 0) = cos(phi) * cos(kappa) - sin(phi) * sin(omega) * sin(kappa);	// a1
		R.at<double>(0, 1) = -cos(phi) * sin(kappa) - sin(phi) * sin(omega) * cos(kappa);	// a2
		R.at<double>(0, 2) = -sin(phi) * cos(omega);										// a3
		R.at<double>(1, 0) = cos(omega) * sin(kappa);										// b1
		R.at<double>(1, 1) = cos(omega) * cos(kappa);										// b2
		R.at<double>(1, 2) = -sin(omega);													// b3
		R.at<double>(2, 0) = sin(phi) * cos(kappa) + cos(phi) * sin(omega) * sin(kappa);	// c1
		R.at<double>(2, 1) = -sin(phi) * sin(kappa) + cos(phi) * sin(omega) * cos(kappa);	// c2
		R.at<double>(2, 2) = cos(phi) * cos(omega);											// c3

		for (int i = 0; i < N; i++) {
			Z0[i] = R.at<double>(0, 2) * (sd[i].X - Xs) + R.at<double>(1, 2) * (sd[i].Y - Ys) + R.at<double>(2, 2) * (sd[i].Z - Zs);	// Z0 = Z_bar = a3(X-Xs)+b3(Y-Ys)+c3(Z-Zs)
			X0[i] = x0 - f * (R.at<double>(0, 0) * (sd[i].X - Xs) + R.at<double>(1, 0) * (sd[i].Y - Ys) + R.at<double>(2, 0) * (sd[i].Z - Zs)) / Z0[i];		// X0 为观测值
			Y0[i] = y0 - f * (R.at<double>(0, 1) * (sd[i].X - Xs) + R.at<double>(1, 1) * (sd[i].Y - Ys) + R.at<double>(2, 1) * (sd[i].Z - Zs)) / Z0[i];		// Y0 为观测值

																																							// 计算误差方程式的系数
			A.at<double>(0, 0) = (R.at<double>(0, 0) * f + R.at<double>(0, 2) * (sd[i].x - x0)) / Z0[i];								// a11
			A.at<double>(0, 1) = (R.at<double>(1, 0) * f + R.at<double>(1, 2) * (sd[i].x - x0)) / Z0[i];								// a12
			A.at<double>(0, 2) = (R.at<double>(2, 0) * f + R.at<double>(2, 2) * (sd[i].x - x0)) / Z0[i];								// a13
			A.at<double>(0, 3) = (sd[i].y - y0) * sin(omega) - ((sd[i].x - x0) * ((sd[i].x - x0) * cos(kappa)
				- (sd[i].y - y0) * sin(kappa)) / f + f * cos(kappa)) * cos(omega);										// a14
			A.at<double>(0, 4) = -f * sin(kappa) - (sd[i].x - x0) * ((sd[i].x - x0) * sin(kappa) + (sd[i].y - y0) * cos(kappa)) / f;	// a15
			A.at<double>(0, 5) = sd[i].y - y0;																							// a16
			A.at<double>(1, 0) = (R.at<double>(0, 1) * f + R.at<double>(0, 2) * (sd[i].y - y0)) / Z0[i];								// a21
			A.at<double>(1, 1) = (R.at<double>(1, 1) * f + R.at<double>(1, 2) * (sd[i].y - y0)) / Z0[i];								// a22
			A.at<double>(1, 2) = (R.at<double>(2, 1) * f + R.at<double>(2, 2) * (sd[i].y - y0)) / Z0[i];								// a23
			A.at<double>(1, 3) = -(sd[i].x - x0) * sin(omega) - ((sd[i].y - y0) * ((sd[i].x - x0) * cos(kappa)
				- (sd[i].y - y0) * sin(kappa)) / f - f * sin(kappa)) * cos(omega);										// a24
			A.at<double>(1, 4) = -f * cos(kappa) - (sd[i].y - y0) * ((sd[i].x - x0) * cos(kappa) - (sd[i].y - y0) * sin(kappa)) / f;	// a25
			A.at<double>(1, 5) = -(sd[i].x - x0);																						// a26

																																		// 该循环保存A矩阵，最后评定精度用
			for (int j = 0; j < 2; j++) {
				for (int k = 0; k < 6; k++) {
					Buf_A.at<double>(2 * i + j, k) = A.at<double>(j, k);
				}
			}

			L.at<double>(0, 0) = sd[i].x - X0[i];
			L.at<double>(1, 0) = sd[i].y - Y0[i];


			// 保存 L 矩阵，最后评定精度用
			for (int j = 0; i < 2; i++) {
				Buf_L.at<double>(2 * i + j, 0) = L.at<double>(j, 0);
			}

			// 所谓的逐步法化，即要在循环内部就将 ATA 计算出来累加，下面的 L 矩阵类似
			ATA = A.t() * A;
			ATL = A.t() * L;
		}// for

		 // “逐步法化”的另一处不同，出循环即可直接计算 ATA 的逆乘 ATL
		correction = ATA.inv() * ATL;

		// correction 即为改正数，每一次迭代时用未知数近似值与上次迭代计算的改正数之和作为新的近似值
		Xs += correction.at<double>(0, 0);
		Ys += correction.at<double>(1, 0);
		Zs += correction.at<double>(2, 0);
		phi += correction.at<double>(3, 0);
		omega += correction.at<double>(4, 0);
		kappa += correction.at<double>(5, 0);

		cout << "【改正数值】：" << endl;
		cout << "\tΔXs = " << correction.at<double>(0, 0) << "\tΔYs = " << correction.at<double>(1, 0) << "\tΔZs = " << correction.at<double>(2, 0) << endl;
		cout << "\tΔφ = " << correction.at<double>(3, 0) << "\tΔω = " << correction.at<double>(4, 0) << "\tΔκ = " << correction.at<double>(5, 0) << endl;

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
		Q[i] = ATA.at<double>(i, i);
	}

	double m0(0);	// 单位权中误差
	double vv(0);	// [vv]，即平方和

	V = Buf_A * correction;
	V = V - Buf_L;

	for (int i = 0; i < 8; i++) {
		vv += pow(V.at<double>(i, 0), 2);
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
void OutputResult(EOEO* eoeo, Mat &R, double* M, double m0) {
	cout << "/////---------- 计算结果 ----------/////" << endl;
	cout << "六个外方位元素为：" << endl;
	printf("Xs = %.4f\n", eoeo->Xs);
	printf("Ys = %.4f\n", eoeo->Ys);
	printf("Zs = %.4f\n", eoeo->Zs);
	printf("φ = %.10f\n", eoeo->phi);
	printf("ω = %.10f\n", eoeo->omega);
	printf("κ = %.10f\n", eoeo->kappa);

	cout << endl << "旋转矩阵：" << endl;
	cout << "\t" << R.at<double>(0, 0) << "\t" << R.at<double>(0, 1) << "\t" << R.at<double>(0, 2) << endl;
	cout << "\t" << R.at<double>(1, 0) << "\t" << R.at<double>(1, 1) << "\t" << R.at<double>(1, 2) << endl;
	cout << "\t" << R.at<double>(2, 0) << "\t" << R.at<double>(2, 1) << "\t" << R.at<double>(2, 2) << endl;

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
	InitInterface();

	SrcData sd[N];
	double m = 15000.0;		// 摄影比例尺 1/m = 1/15000
	double f = 153.24;		// 主距 f = 153.24mm

	char filename[50] = "data.txt";	// 储存已知坐标数据的文件名为 data.txt
	InitData(sd, filename);
	Iterator(sd, m, f);

	return 0;
}

