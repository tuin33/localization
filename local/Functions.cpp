#include "Functions.h"
#include <iostream>
#include <map>
#include <algorithm>
#include <mat.h>
using namespace std;

//排序规则结构体定义
struct cmpByValue
{
	bool operator()(const pair<int, double>& leftPair, const pair<int, double>& rightPair)
	{
		return leftPair.second < rightPair.second;
	}
};
//获取向量中最大值 的一个函数
int mymax(const vector<int>& v)
{
	int m = v.front();
	for (auto x : v)
		m = max(m, x);
	return m;
}
//正态分布生成函数，生成期望为mu，标准差为sigma的正态分布
double gaussrand(const double& mu, const double& sigma)
{
	static double V1, V2, S;
	static int phase = 0;
	double X;
	if (phase == 0)
	{
		do
		{
			double U1 = (double)rand() / RAND_MAX;
			double U2 = (double)rand() / RAND_MAX;
			V1 = 2 * U1 - 1;
			V2 = 2 * U2 - 1;
			S = V1 * V1 + V2 * V2;
		} while (S >= 1 || S == 0);
		X = V1 * sqrt(-2 * log(S) / S);
	}
	else
		X = V2 * sqrt(-2 * log(S) / S);
	phase = 1 - phase;
	X = X * sigma + mu;
	return X;
}
//计算距离，v1与v2是两个大小相同的向量
double distance(const vector<double>& v1, const vector<double>& v2)
{
	double sum = 0;
	for (int i = 0; i < v1.size(); ++i)
		sum += (v1[i] - v2[i]) * (v1[i] - v2[i]);
	return sqrt(sum);
}
//实现k-临近算法
vector<int> knnsearch(vector<vector<double>> testData, vector<double> targetData, int k)
{
	map<int, double> mapIndexDistance;
	for (int i = 0; i < testData.size(); ++i)
		mapIndexDistance[i] = distance(testData[i], targetData);
	vector<pair<int, double>> pairIndexDistance(mapIndexDistance.begin(), mapIndexDistance.end());
	sort(pairIndexDistance.begin(), pairIndexDistance.end(), cmpByValue());
	vector<int> vec(k);
	for (int i = 0; i < k; ++i)
		vec[i] = pairIndexDistance[i].first;
	return vec;
}
//判断点是否落入多边形内
vector<bool> inpolygon(const vector<vector<double>>& pointData, const vector<int>& optionalTagFlag, const vector<vector<double>>& testData,const int& numOfTag)
{
	vector<double> xPoint(numOfTag);
	vector<double> yPoint(numOfTag);
	for (int k = 0; k < numOfTag; ++k)
	{
		xPoint[k] = pointData[optionalTagFlag[k]][0];
		yPoint[k] = pointData[optionalTagFlag[k]][1];
	}
	vector<double> xLine(4);
	vector<double> yLine(4);
	for (int k = 0; k < 4; ++k)
	{
		xLine[k] = testData[k][2];
		yLine[k] = testData[k][3];
	}

	//开始k-近邻算法
	vector<bool> result(numOfTag);
	for (int i = 0; i < numOfTag; ++i)
	{
		bool flag = false;	   //判断结果（true；点落在多边形内；false:点未落在多边形内）
		int k = testData.size() - 1; //是多边形的最后一个顶点
		for (int j = 0; j < testData.size(); ++j)
		{
			//判断点是否在线段的两侧
			if ((yLine[j] < yPoint[i] && yLine[k] >= yPoint[i]) || (yLine[k] < yPoint[i] && yLine[j] >= yPoint[i]))
			{
				//根据两点式方程计算出过点P且平行于X轴的直线与线段的交点，两点式方程：x = x1 +  (y - y1) * (x2 - x1) / (y2 - y1);
				if (xLine[j] + (yPoint[i] - yLine[j]) * (xLine[k] - xLine[j]) / (yLine[k] - yLine[j]) < xPoint[i])
					flag = !flag;
			}
			//进行下一线段判断
			k = j;
		}
		result[i] = flag;
	}
	return result;
}
//找到里程计的时间在视觉时间序列中的前后两个最近的位姿点
int findpso(const vector<double>& visionTimeDiff, const vector<double>& odometerTimeDiff, const int& index, bool& flag)
{
	int point = 0;
	for (auto it = visionTimeDiff.begin(); it != visionTimeDiff.end(); ++it)
	{
		if ((*it) >= odometerTimeDiff[index] && !flag)
		{
			point = it - visionTimeDiff.begin();
			flag = true;
		}
		if (flag)
			break;
	}
	return point;
}
int findpre(const vector<double>& visionTimeDiff, const vector<double>& odometerTimeDiff, const int& index, bool& flag)
{
	int point = 0;
	for (auto it = visionTimeDiff.end() - 1; it >= visionTimeDiff.begin(); --it)
	{
		if (((*it) <= odometerTimeDiff[index]) && !flag)
		{
			point = it - visionTimeDiff.begin();
			flag = true;
		}
		if (flag)
			break;
	}
	return point;
}
//计算机器人所在大致范围（用于初始化粒子）
vector<vector<double>> calPF_Scope(vector<vector<double>>& referenceTag, vector<int>& randData, vector<int>& colLeft, vector<int>& colRight)
{
	vector<vector<double>> PF_Scope(2, vector<double>(2));
	double leftMax = referenceTag[colLeft[0]][0];
	for (auto x : colLeft)
		leftMax = max(leftMax, referenceTag[x][0]);
	double rightMax = referenceTag[colRight[0]][0];
	for (auto x : colRight)
		rightMax = max(rightMax, referenceTag[x][0]);
	PF_Scope[0][0] = max(leftMax, rightMax) - 160 + randData[0];

	double leftMin = referenceTag[colLeft[0]][0];
	for (auto x : colLeft)
		leftMin = min(leftMin, referenceTag[x][0]);
	double rightMin = referenceTag[colRight[0]][0];
	for (auto x : colRight)
		rightMin = min(rightMin, referenceTag[x][0]);
	PF_Scope[0][1] = min(leftMin, rightMin) + 160 + randData[2];

	leftMax = referenceTag[colLeft[0]][1];
	for (auto x : colLeft)
		leftMax = max(leftMax, referenceTag[x][1]);
	rightMax = referenceTag[colRight[0]][1];
	for (auto x : colRight)
		rightMax = max(rightMax, referenceTag[x][1]);
	PF_Scope[1][0] = max(leftMax, rightMax) - 160 + randData[1];

	leftMin = referenceTag[colLeft[0]][1];
	for (auto x : colLeft)
		leftMin = min(leftMin, referenceTag[x][1]);
	rightMin = referenceTag[colRight[0]][1];
	for (auto x : colRight)
		rightMin = min(rightMin, referenceTag[x][1]);
	PF_Scope[1][1] = min(leftMin, rightMin) + 160 + randData[3];

	return PF_Scope;
}
//读取.mat文件
//void readmat(const char* matPath, const char** fieldName, double*** q,int len)
//{
//	//视觉数据
//	//前期工作
//	//定义文件路径
//	//const char* p = R"(E:\MATLAB workspace\hgdw\0319_3.mat)";
//	//定义域名
//	//const char* name[] = { "AGV_relative", "AtennaL_POS", "AtennaR_POS", "AGV11", "TAG_POS" };
//	//打开视觉数据文件
//	MATFile* pM = matOpen(matPath, "r");
//	//获取名为"my_vision_total"的变量，它在matlab中是一个结构体
//	mxArray* pS = matGetVariable(pM, "my_vision_total");
//	//MyVisionTotal myVisionTotal;
//	//结构体首元素的地址
//	//double*** q = &(myVisionTotal.AGV_Relative);
//
//	//循环存储
//	for (int i = 0; i < len; ++i)
//	{
//		//打开对应的域，这里的PF存放的是各field的指针
//		mxArray* pF = mxGetField(pS, 0, fieldName[i]);
//		//把域中的数据存入一个临时变量
//		double* temp = (double*)mxGetData(pF);
//		//获取数据的行数
//		auto row = mxGetM(pF);
//		//获取数据的列数
//		auto col = mxGetN(pF);
//		//新建指针数组，这还是个二级指针，这个数组内存放的是每一行的首地址
//		*(q + i) = new double* [row];
//		for (int j = 0; j < row; ++j)
//		{
//			//新建 double 类型的数组
//			*(*(q + i) + j) = new double[col];
//			for (int k = 0; k < col; ++k)
//				//把临时变量的值赋给结构体的变量
//				*(*(*(q + i) + j) + k) = temp[k * row + j];
//		}
//	}
//	//至此，数据已按照.mat文件中的排布全部导入进了C++中的myVisionTotal
//}