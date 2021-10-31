#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <map>
#include <numeric>
#include <algorithm>
#include <malloc.h>
#include <mat.h>
#include <windows.h>
// #include <unistd.h>
using namespace std;
struct Data
{
	vector<int> no;			//1
	vector<double> phase;		//2
	vector<double> myThree;		//3
	vector<double> myFour;		//4
	vector<double> Xcoordinate; 	//5
	vector<double> Ycoordinate; 	//6
	vector<double> Zcoordinate; 	//7
	vector<double> readerTime;	//8
	vector<double> windowsTime; 	//9
	vector<double> myTen;		//10
	vector<int> RSSI;		//11 接收的信号强度

struct TextData
{
	vector<int> no;
	vector<string> epcData;
};
struct MyDataResult
{
	Data data;
	TextData textData;
};
struct MyVisionTotal
{
	double** AGV_Relative;
	double** AntennaL_POS;
	double** AntennaR_POS;
	double** AGV11;
	double** TAG_POS;
};
//排序规则结构体定义
struct cmpByValue
{
	bool operator()(const pair<int, double>& leftPair, const pair<int, double>& rightPair)
	{
		return leftPair.second < rightPair.second;
	}
};
//数据精度测试阈值
const double ep = 0.0001;
//PI
const double PI = 3.1415926538;
//这个没有用到
const double Resolution = 0.0015;
//光速
const int vLight = 300000000;
//单位 raidan/s,标准差，应该是观测噪声的标准差
const double PhaseGauss = 0.1; 
//文件路径，可以把需要的文件放在工程文件夹中，简化路径的定义
const string MyDataLeftAntenna = R"(E:\MATLAB workspace\hgdw\48\myData - leftAntenna.txt)";
const string MyDataRightAntenna = R"(E:\MATLAB workspace\hgdw\48\myData - rightAntenna.txt)";
const string OdometerVisionStartPoint = R"(E:\MATLAB workspace\hgdw\48\odometer_vision_start_point.txt)";
//粒子滤波循环次数，虽然并没有循环这么多次
const int Times = 400;
//获取向量中最大值 的一个函数
int mymax(const vector<int>&);
//正态分布生成函数，生成期望为mu，标准差为sigma的正态分布
double gaussrand(const double& mu, const double& sigma);
//计算距离，v1与v2是两个大小相同的向量
double distance(const vector<double>& v1, const vector<double>& v2);
//实现k-临近算法
vector<int> knnsearch(vector<vector<double>> testData, vector<double> targetData, int k);
//判断点是否落入多边形内
vector<bool> inpolygon(const vector<double>& pointX, const vector<double>& pointY, const vector<double>& lineX, const vector<double>& lineY);
int main()
{
	//系统数据设定
	double antennaHLeft = 0;
	double antennaHRight = 0;
	double antennaHLeftError = antennaHLeft + 0;
	double antennaHRightError = antennaHRight - 0;
	double antennaAlphaLeft = 0;
	double antennaAlphaRight = 0;
	double antennaAlphaLeftError = antennaAlphaLeft + 0 / 180 * PI;
	double antennaAlphaRightError = antennaAlphaRight + 0 / 180 * PI;
	vector<double> waveLengthVar;
	vector<double> frequencyVar;
	for (double l = 920.625; l <= 924.375; l = l + 0.25)
	{
		double f = vLight / l * 0.0001;
		waveLengthVar.push_back(f);
		frequencyVar.push_back(l);
	}

	//系统数据导入
	//RFID和机器人数据
	fstream fin;

	fin.open(MyDataLeftAntenna, ios::in);
	MyDataResult myDataResultLeft;
	while (1)
	{
		int n;
		string s;
		double d;
		//到达文件结尾需再次读入，eof才会变为true
		fin >> n;
		if (fin.eof())
			break;
		myDataResultLeft.textData.no.push_back(n);
		fin >> s;
		myDataResultLeft.textData.epcData.push_back(s);
		fin >> n;
		myDataResultLeft.data.no.push_back(n);
		fin >> d;
		myDataResultLeft.data.phase.push_back(d);
		fin >> d;
		myDataResultLeft.data.myThree.push_back(d);
		fin >> d;
		myDataResultLeft.data.myFour.push_back(d);
		fin >> d;
		myDataResultLeft.data.Xcoordinate.push_back(d);
		fin >> d;
		myDataResultLeft.data.Ycoordinate.push_back(d);
		fin >> d;
		myDataResultLeft.data.Zcoordinate.push_back(d);
		fin >> d;
		myDataResultLeft.data.readerTime.push_back(d);
		fin >> d;
		myDataResultLeft.data.windowsTime.push_back(d);
		fin >> d;
		myDataResultLeft.data.myTen.push_back(d);
		fin >> n;
		myDataResultLeft.data.RSSI.push_back(n);
	}
	fin.close();

	fin.open(MyDataRightAntenna, ios::in);
	MyDataResult myDataResultRight;
	while (1)
	{
		int n;
		string s;
		double d;
		fin >> n;
		if (fin.eof())
			break;
		myDataResultRight.textData.no.push_back(n);
		fin >> s;
		myDataResultRight.textData.epcData.push_back(s);
		fin >> n;
		myDataResultRight.data.no.push_back(n);
		fin >> d;
		myDataResultRight.data.phase.push_back(d);
		fin >> d;
		myDataResultRight.data.myThree.push_back(d);
		fin >> d;
		myDataResultRight.data.myFour.push_back(d);
		fin >> d;
		myDataResultRight.data.Xcoordinate.push_back(d);
		fin >> d;
		myDataResultRight.data.Ycoordinate.push_back(d);
		fin >> d;
		myDataResultRight.data.Zcoordinate.push_back(d);
		fin >> d;
		myDataResultRight.data.readerTime.push_back(d);
		fin >> d;
		myDataResultRight.data.windowsTime.push_back(d);
		fin >> d;
		myDataResultRight.data.myTen.push_back(d);
		fin >> n;
		myDataResultRight.data.RSSI.push_back(n);
	}
	fin.close();

	fin.open(OdometerVisionStartPoint, ios::in);
	//先定义vector再尾插
	vector<int> odometerVisionStartPoint;
	{
		int OVSP;
		while (fin >> OVSP)
			odometerVisionStartPoint.push_back(OVSP);
	}
	fin.close();//皓宇哥哥细节啊

	//以下两段代码用于建立所有的EPC，不难看出共有36个
	vector<int> referenceTag_EPC_No;
	for (int i = 10; i <= 45; ++i)
		referenceTag_EPC_No.push_back(i);

	vector<string> referenceTag_EPC;
	//另一种遍历方式罢了
	for (auto x : referenceTag_EPC_No)
		referenceTag_EPC.push_back("1037-9654-FFFF-FFFF-FFFF-00" + to_string(x));

	int visionStartPoint = odometerVisionStartPoint[0];
	int odometerRightStartPoint = odometerVisionStartPoint[1];
	int odometerLeftStartPoint = odometerVisionStartPoint[2];

	//视觉数据
	//前期工作
	//定义文件路径
	const char* p = R"(E:\MATLAB workspace\hgdw\0319_3.mat)";
	//定义域名
	const char* name[] = { "AGV_relative", "AtennaL_POS", "AtennaR_POS", "AGV11", "TAG_POS" };
	//打开视觉数据文件
	MATFile* pM = matOpen(p, "r");
	//获取名为"my_vision_total"的变量，它在matlab中是一个结构体
	mxArray* pS = matGetVariable(pM, "my_vision_total");
	MyVisionTotal myVisionTotal;
	//结构体首元素的地址
	double*** q = &(myVisionTotal.AGV_Relative);

	//循环存储
	for (int i = 0; i < sizeof(myVisionTotal) / sizeof(double**); ++i)
	{
		//打开对应的域，这里的PF存放的是各field的指针
		mxArray* pF = mxGetField(pS, 0, name[i]);
		//把域中的数据存入一个临时变量
		double* temp = (double*)mxGetData(pF);
		//获取数据的行数
		auto row = mxGetM(pF);
		//获取数据的列数
		auto col = mxGetN(pF);
		//新建指针数组，这还是个二级指针，这个数组内存放的是每一行的首地址
		*(q + i) = new double* [row];
		for (int j = 0; j < row; ++j)
		{
			//新建 double 类型的数组
			*(*(q + i) + j) = new double[col];
			for (int k = 0; k < col; ++k)
				//把临时变量的值赋给结构体的变量
				*(*(*(q + i) + j) + k) = temp[k * row + j];
		}
	}
	//至此，数据已按照.mat文件中的排布全部导入进了C++中的myVisionTotal

	/*数据检测*/
	// for (int i = 0; i < sizeof(myVisionTotal) / sizeof(double **); i++)
	// {
	// 	mxArray *pF = mxGetField(pS, 0, name[i]);
	// 	auto row = mxGetM(pF);
	// 	auto col = mxGetN(pF);
	// 	for (int j = 0; j < row; j++)
	// 	{
	// 		for (int k = 0; k < col; k++)
	// 			cout << *(*(*(q + i) + j) + k) << ' ';
	// 		cout << endl;
	// 	}
	// }
	/**************/

	//_msize()函数用于获取new出来内存的大小
	//逗号前后分别用来获取行数和列数，注意设定列数时vector<double>不能少
	vector<vector<double>> myVisionCurrent_AGV(_msize(myVisionTotal.AGV11) / sizeof(myVisionTotal.AGV11[0]), vector<double>(_msize(myVisionTotal.AGV11[0]) / sizeof(myVisionTotal.AGV11[0][0])));
	for (int i = 0; i < _msize(myVisionTotal.AGV11) / sizeof(myVisionTotal.AGV11[0]); ++i)
		for (int j = 0; j < _msize(myVisionTotal.AGV11[0]) / sizeof(myVisionTotal.AGV11[0][0]); ++j)
			myVisionCurrent_AGV[i][j] = myVisionTotal.AGV11[i][j];

	//数组转置了一次
	vector<vector<double>> referenceTag(_msize(myVisionTotal.TAG_POS[0]) / sizeof(myVisionTotal.TAG_POS[0][0]), vector<double>(_msize(myVisionTotal.TAG_POS) / sizeof(myVisionTotal.TAG_POS[0])));
	for (int i = 0; i < _msize(myVisionTotal.TAG_POS[0]) / sizeof(myVisionTotal.TAG_POS[0][0]); ++i)
		for (int j = 0; j < _msize(myVisionTotal.TAG_POS) / sizeof(myVisionTotal.TAG_POS[0]); ++j)
			referenceTag[i][j] = myVisionTotal.TAG_POS[j][i] * 100;

	//这里用的是转置前的矩阵，于是求标签数的话就求列数就好（因为这时还没有转置）
	int referenceTagNum = _msize(myVisionTotal.TAG_POS[0]) / sizeof(myVisionTotal.TAG_POS[0][0]);

	for (int i = 0; i < 3; i++)
		for (int j = 0; j < myVisionCurrent_AGV[0].size(); ++j)
			myVisionCurrent_AGV[i][j] *= 100;

	double myVisionRightAntennaHigh = myVisionTotal.AntennaR_POS[1][0] * 100;
	double myVisionLeftAntennaHigh = myVisionTotal.AntennaL_POS[1][0] * 100;

	//1 2 3 4 -> 3 2 1 4
	myVisionCurrent_AGV[0].swap(myVisionCurrent_AGV[2]);
	//3 2 1 4 -> 3 4 1 2
	myVisionCurrent_AGV[1].swap(myVisionCurrent_AGV[3]);
	//3 4 1 2 -> 3 1 4 2
	myVisionCurrent_AGV[1].swap(myVisionCurrent_AGV[2]);

	//又是这种遍历方法，若不加引用则无法更改
	for (auto& vec : referenceTag)
	{
		//注意每一层循环中，vec都是referenceTag中的某一行。故循环共需要执行行数次
		auto temp = vec[2];
		vec[2] = vec[1];
		vec[1] = vec[0];
		vec[0] = temp;
	}

	//计算天线在移动机器人坐标系下的坐标
	//不难发现relative的[3][0]即第四个数是角度
	vector<double> antennaR_Robot(2);
	antennaR_Robot[0] = (myVisionTotal.AntennaR_POS[2][0] - myVisionTotal.AGV_Relative[2][0]) * cos(-myVisionTotal.AGV_Relative[3][0]) - (myVisionTotal.AntennaR_POS[0][0] - myVisionTotal.AGV_Relative[0][0]) * sin(-myVisionTotal.AGV_Relative[3][0]); //机器人下的x坐标
	antennaR_Robot[1] = (myVisionTotal.AntennaR_POS[2][0] - myVisionTotal.AGV_Relative[2][0]) * sin(-myVisionTotal.AGV_Relative[3][0]) + (myVisionTotal.AntennaR_POS[0][0] - myVisionTotal.AGV_Relative[0][0]) * cos(-myVisionTotal.AGV_Relative[3][0]); //机器人下的y坐标

	vector<double> antennaL_Robot(2);
	antennaL_Robot[0] = (myVisionTotal.AntennaL_POS[2][0] - myVisionTotal.AGV_Relative[2][0]) * cos(-myVisionTotal.AGV_Relative[3][0]) - (myVisionTotal.AntennaL_POS[0][0] - myVisionTotal.AGV_Relative[0][0]) * sin(-myVisionTotal.AGV_Relative[3][0]); //机器人下的x坐标
	antennaL_Robot[1] = (myVisionTotal.AntennaL_POS[2][0] - myVisionTotal.AGV_Relative[2][0]) * sin(-myVisionTotal.AGV_Relative[3][0]) + (myVisionTotal.AntennaL_POS[0][0] - myVisionTotal.AGV_Relative[0][0]) * cos(-myVisionTotal.AGV_Relative[3][0]); //机器人下的y坐标

	//天线高度和角度
	//勾股定理
	//下述角度与高度为何？
	antennaHRight = sqrt(antennaR_Robot[0] * antennaR_Robot[0] + antennaR_Robot[1] * antennaR_Robot[1]) * 100;
	antennaHRightError = antennaHRight - 0;
	antennaAlphaRight = atan(abs(antennaR_Robot[0] / antennaR_Robot[1])) + PI / 2; //45 / 180 * PI
	antennaAlphaRightError = antennaAlphaRight + 0 / 180 * PI;
	antennaHLeft = sqrt(antennaL_Robot[0] * antennaL_Robot[0] + antennaL_Robot[1] * antennaL_Robot[1]) * 100;
	antennaHLeftError = antennaHLeft + 0;
	antennaAlphaLeft = atan(abs(antennaL_Robot[0] / antennaL_Robot[1])) + PI / 2; //45 / 180 * PI
	antennaAlphaLeftError = antennaAlphaLeft + 0 / 180 * PI;

	//绘制 - 未完成

	//读取，并对测试数据分类
	//提取左右天线读取标签的epc、该epc标签的读取次数、标签的数量
	for (auto& x : myDataResultLeft.data.phase)
		//将相位转化为弧度
		x = x / 180 * PI;

	//textData中存放的是读取的次数和每一次读取到的epc，下面一段代码只是一个提取操作
	vector<string> epcDataLeft(myDataResultLeft.textData.epcData.size());
	for (int i = 0; i < epcDataLeft.size(); ++i)
		epcDataLeft[i] = myDataResultLeft.textData.epcData[i];

	//用来统计每个标签被读取的次数
	vector<int> referenceTagCountLeft(referenceTagNum);
	vector<vector<int>> referenceTagNumLeft;
	for (int i = 0; i < myDataResultLeft.data.no.size(); ++i)
		for (int j = 0; j < referenceTagNum; ++j)
			if (referenceTag_EPC[j] == epcDataLeft[i])
			{
				++referenceTagCountLeft[j];
				if (referenceTagNumLeft.size() < referenceTagCountLeft[j])
				{
					vector<int> v(referenceTagNum);
					referenceTagNumLeft.push_back(v);
				}
				// -1是为了在vector规定的内存里
				referenceTagNumLeft[referenceTagCountLeft[j] - 1][j] = i;
				break;
			}

	for (auto& x : myDataResultRight.data.phase)
		x = x / 180 * PI;

	vector<string> epcDataRight(myDataResultRight.textData.epcData.size());
	for (int i = 0; i < epcDataRight.size(); ++i)
		epcDataRight[i] = myDataResultRight.textData.epcData[i];

	vector<int> referenceTagCountRight(referenceTagNum);
	vector<vector<int>> referenceTagNumRight;
	for (int i = 0; i < myDataResultRight.data.no.size(); ++i)
		for (int j = 0; j < referenceTagNum; ++j)
			if (referenceTag_EPC[j] == epcDataRight[i])
			{
				referenceTagCountRight[j] = referenceTagCountRight[j] + 1;
				if (referenceTagNumRight.size() < referenceTagCountRight[j])
				{
					vector<int> v(referenceTagNum);
					referenceTagNumRight.push_back(v);
				}
				// -1是为了在vector规定的内存里
				referenceTagNumRight[referenceTagCountRight[j] - 1][j] = i;
				break;
			}

	//时间
	vector<double> readerTimeLeft = myDataResultLeft.data.readerTime;
	vector<double> windowsTimeLeft = myDataResultLeft.data.windowsTime;
	vector<double> readerTimeRight = myDataResultRight.data.readerTime;
	vector<double> windowsTimeRight = myDataResultRight.data.windowsTime;

	//标签读取速度
	double tagReaderRate = myDataResultLeft.data.no.size() / (myDataResultLeft.data.readerTime.back() - myDataResultLeft.data.readerTime.front()) * 1000000;

	//相位预处理
	double indexLb = 0.5; //0.5;%1/2 3/4 2/3
	double indexUb = 1.5; //3/2 5/4 4/3

	vector<vector<double>> phaseMeasurementAdjustLeft(mymax(referenceTagCountLeft), vector<double>(referenceTagNum));
	for (int j = 0; j < referenceTagNum; ++j)
		if (referenceTagCountLeft[j] > 1)
		{
			phaseMeasurementAdjustLeft[0][j] = myDataResultLeft.data.phase[referenceTagNumLeft[0][j]];
			for (int i = 1; i < referenceTagCountLeft[j]; ++i)
			{
				auto dif = myDataResultLeft.data.phase[referenceTagNumLeft[i][j]] - phaseMeasurementAdjustLeft[i - 1][j];
				if ((dif > PI * indexLb) && (dif < PI * indexUb))
					phaseMeasurementAdjustLeft[i][j] = myDataResultLeft.data.phase[referenceTagNumLeft[i][j]] - PI;
				else if ((-dif > PI * indexLb) && (-dif < PI * indexUb))
					phaseMeasurementAdjustLeft[i][j] = myDataResultLeft.data.phase[referenceTagNumLeft[i][j]] + PI;
				else
					phaseMeasurementAdjustLeft[i][j] = myDataResultLeft.data.phase[referenceTagNumLeft[i][j]];
			}
		}

	if (accumulate(referenceTagCountLeft.begin(), referenceTagCountLeft.end(), 0) > 0)
		for (auto& vec : phaseMeasurementAdjustLeft)
			for (auto& x : vec)
				x = 2 * PI - x;

	vector<vector<double>> phasePsoLeft(mymax(referenceTagCountLeft), vector<double>(referenceTagNum));
	for (int j = 0; j < referenceTagNum; ++j)
		if (referenceTagCountLeft[j] > 1)
		{
			phasePsoLeft[0][j] = phaseMeasurementAdjustLeft[0][j];
			for (int i = 1; i < referenceTagCountLeft[j]; i++)
			{
				int k = round((phaseMeasurementAdjustLeft[i][j] - phasePsoLeft[i - 1][j]) / PI);
				phasePsoLeft[i][j] = phaseMeasurementAdjustLeft[i][j] - PI * k;
			}
		}

	vector<vector<double>> phaseMeasurementAdjustRight(mymax(referenceTagCountRight), vector<double>(referenceTagNum));
	for (int j = 0; j < referenceTagNum; ++j)
		if (referenceTagCountRight[j] > 1)
		{
			phaseMeasurementAdjustRight[0][j] = myDataResultRight.data.phase[referenceTagNumRight[0][j]];
			for (int i = 1; i < referenceTagCountRight[j]; ++i)
			{
				auto dif = myDataResultRight.data.phase[referenceTagNumRight[i][j]] - phaseMeasurementAdjustRight[i - 1][j];
				if ((dif > PI * indexLb) && (dif < PI * indexUb))
					phaseMeasurementAdjustRight[i][j] = myDataResultRight.data.phase[referenceTagNumRight[i][j]] - PI;
				else if ((-dif > PI * indexLb) && (-dif < PI * indexUb))
					phaseMeasurementAdjustRight[i][j] = myDataResultRight.data.phase[referenceTagNumRight[i][j]] + PI;
				else
					phaseMeasurementAdjustRight[i][j] = myDataResultRight.data.phase[referenceTagNumRight[i][j]];
			}
		}

	for (auto& vec : phaseMeasurementAdjustRight)
		for (auto& x : vec)
			x = 2 * PI - x;

	vector<vector<double>> phasePsoRight(mymax(referenceTagCountRight), vector<double>(referenceTagNum));
	for (int j = 0; j < referenceTagNum; ++j)
		if (referenceTagCountRight[j] > 1)
		{
			phasePsoRight[0][j] = phaseMeasurementAdjustRight[0][j];
			for (int i = 1; i < referenceTagCountRight[j]; ++i)
			{
				auto k = round((phaseMeasurementAdjustRight[i][j] - phasePsoRight[i - 1][j]) / PI);
				phasePsoRight[i][j] = phaseMeasurementAdjustRight[i][j] - PI * k;
			}
		}

	// 右天线：利用时间差信息，通过插值法得到相位测量时刻真实的移动机器人和天线的位姿/位置
	vector<vector<double>> trackMobileRobotRight;
	trackMobileRobotRight.push_back(myDataResultRight.data.Xcoordinate);
	trackMobileRobotRight.push_back(myDataResultRight.data.Ycoordinate);
	trackMobileRobotRight.push_back(myDataResultRight.data.Zcoordinate);

	//将里程计得到的机器人航迹分发
	vector<vector<vector<double>>> trackMobileRobotRightTag(referenceTagNum, vector<vector<double>>(mymax(referenceTagCountRight), vector<double>(3)));
	for (int j = 0; j < referenceTagNum; ++j)
		if (referenceTagCountRight[j] > 1)
			for (int i = 0; i < referenceTagCountRight[j]; ++i)
				for (int k = 0; k < 3; ++k)
					trackMobileRobotRightTag[j][i][k] = trackMobileRobotRight[k][referenceTagNumRight[i][j]];

	//计算出视觉的时间差
	vector<double> visionTimeDiff(myVisionCurrent_AGV[0].size());
	for (int i = 0; i < myVisionCurrent_AGV[0].size() - visionStartPoint + 1; ++i)
		visionTimeDiff[i + visionStartPoint - 1] = 1.0 / 120 * (i + 1);

	//去除NaN
	for (int i = 0; i < myVisionCurrent_AGV[0].size();)
		if (isnan(myVisionCurrent_AGV[0][i]))
		{
			for (auto& vec : myVisionCurrent_AGV)
				vec.erase(vec.begin() + i, vec.begin() + i + 1);
			visionTimeDiff.erase(visionTimeDiff.begin() + i, visionTimeDiff.begin() + i + 1);
		}
		else
			++i;

	//纠正视觉所获得的移动机器人朝向角数据
	myVisionCurrent_AGV[2][0] = myVisionCurrent_AGV[2][0];
	for (int i = 1; i < myVisionCurrent_AGV[2].size(); ++i)
	{
		auto k = round((myVisionCurrent_AGV[2][i] - myVisionCurrent_AGV[2][i - 1]) / PI);
		myVisionCurrent_AGV[2][i] = myVisionCurrent_AGV[2][i] - PI * k;
	}

	//计算出里程计的时间差
	vector<double> odometerTimeDiff(readerTimeRight.size() - odometerRightStartPoint + 1);
	for (int i = 0; i < odometerTimeDiff.size(); ++i)
		odometerTimeDiff[i] = (readerTimeRight[odometerRightStartPoint - 1 + i] - readerTimeRight[odometerRightStartPoint - 1]) / 1000000;

	//通过插值法获得移动机器人位姿点
	vector<vector<double>> mobileRobotPoseVisionRight(odometerTimeDiff.size() + odometerRightStartPoint - 1, vector<double>(4));
	for (int i = 0; i < odometerRightStartPoint - 1; ++i)
		for (int j = 0; j < myVisionCurrent_AGV.size(); ++j)
			mobileRobotPoseVisionRight[i][j] = myVisionCurrent_AGV[j][0];
	for (int j = 0; j < odometerTimeDiff.size(); ++j)
	{
		int pointPre = 0;
		int pointPso = 0;
		bool preFlag = false;
		bool psoFlag = false;
		//先找到 里程计的时间 在 视觉时间序列中 的 前后两个最近 的位姿点
		for (auto it = visionTimeDiff.begin(); it != visionTimeDiff.end(); ++it)
		{
			if ((*it) >= odometerTimeDiff[j] && !psoFlag)
			{
				pointPso = it - visionTimeDiff.begin();
				psoFlag = true;
			}
			if (psoFlag)
				break;
		}
		for (auto it = visionTimeDiff.end() - 1; it >= visionTimeDiff.begin(); --it)
		{
			if (((*it) <= odometerTimeDiff[j]) && !preFlag)
			{
				pointPre = it - visionTimeDiff.begin();
				preFlag = true;
			}
			if (preFlag)
				break;
		}
		if (psoFlag)
			for (int i = 0; i < 4; ++i)
				mobileRobotPoseVisionRight[j + odometerRightStartPoint - 1][i] = (odometerTimeDiff[j] - visionTimeDiff[pointPso]) / (visionTimeDiff[pointPre] - visionTimeDiff[pointPso]) * myVisionCurrent_AGV[i][pointPre] + (odometerTimeDiff[j] - visionTimeDiff[pointPre]) / (visionTimeDiff[pointPso] - visionTimeDiff[pointPre]) * myVisionCurrent_AGV[i][pointPso];
	}

	//右天线：由上述插值得到的机器人实际路径，得到理论的相位序列
	vector<vector<vector<double>>> trackMobileRobotRightAntennaReal(referenceTagNum, vector<vector<double>>(mymax(referenceTagCountRight), vector<double>(2)));
	vector<vector<double>> distanceTheory(mymax(referenceTagCountRight), vector<double>(referenceTagNum));
	vector<vector<double>> phaseTagTheoryRight(mymax(referenceTagCountRight), vector<double>(referenceTagNum));
	for (int j = 0; j < referenceTagNum; ++j)
		if (referenceTagCountRight[j] > 1)
			for (int i = 0; i < referenceTagCountRight[j]; ++i)
			{
				trackMobileRobotRightAntennaReal[j][i][0] = mobileRobotPoseVisionRight[referenceTagNumRight[i][j]][0] + antennaHRightError * cos(-antennaAlphaRightError + mobileRobotPoseVisionRight[referenceTagNumRight[i][j]][2]);
				trackMobileRobotRightAntennaReal[j][i][1] = mobileRobotPoseVisionRight[referenceTagNumRight[i][j]][1] + antennaHRightError * sin(-antennaAlphaRightError + mobileRobotPoseVisionRight[referenceTagNumRight[i][j]][2]);
				distanceTheory[i][j] = sqrt(pow((trackMobileRobotRightAntennaReal[j][i][0] - referenceTag[j][0]), 2) + pow((trackMobileRobotRightAntennaReal[j][i][1] - referenceTag[j][1]), 2) + pow((myVisionRightAntennaHigh - referenceTag[j][2]), 2)) * 2;
				phaseTagTheoryRight[i][j] = distanceTheory[i][j] * 2 * PI / waveLengthVar[0];
			}

	//右天线：信号强度分析
	vector<vector<int>> RSSI_TagRight(mymax(referenceTagCountRight), vector<int>(referenceTagNum));
	for (int j = 0; j < referenceTagNum; ++j)
		if (referenceTagCountRight[j] > 1)
			for (int i = 0; i < referenceTagCountRight[j]; ++i)
				RSSI_TagRight[i][j] = myDataResultRight.data.RSSI[referenceTagNumRight[i][j]];

	//左天线：利用时间差信息，通过插值法得到相位测鲁时刻真实的移动机器人和天线的位姿/位置
	vector<vector<double>> trackMobileRobotLeft;
	trackMobileRobotLeft.push_back(myDataResultLeft.data.Xcoordinate);
	trackMobileRobotLeft.push_back(myDataResultLeft.data.Ycoordinate);
	trackMobileRobotLeft.push_back(myDataResultLeft.data.Zcoordinate);

	//里程计得到的机器人航迹分发
	vector<vector<vector<double>>> trackMobileRobotLeftTag(referenceTagNum, vector<vector<double>>(mymax(referenceTagCountLeft), vector<double>(3)));
	for (int j = 0; j < referenceTagNum; ++j)
		if (referenceTagCountLeft[j] > 1)
			for (int i = 0; i < referenceTagCountLeft[j]; ++i)
				for (int k = 0; k < 3; ++k)
					trackMobileRobotLeftTag[j][i][k] = trackMobileRobotLeft[k][referenceTagNumLeft[i][j]];

	//计算出里程计的时间差
	odometerTimeDiff.clear();
	odometerTimeDiff.resize(readerTimeLeft.size() - odometerLeftStartPoint + 1);
	for (int i = 0; i < odometerTimeDiff.size(); ++i)
		odometerTimeDiff[i] = (readerTimeLeft[odometerLeftStartPoint - 1 + i] - readerTimeLeft[odometerLeftStartPoint - 1]) / 1000000;

	//通过插值法获得移动机器人位姿点
	vector<vector<double>> mobileRobotPoseVisionLeft(odometerTimeDiff.size() + odometerLeftStartPoint - 1, vector<double>(4));
	for (int i = 0; i < odometerLeftStartPoint - 1; ++i)
		for (int j = 0; j < myVisionCurrent_AGV.size(); ++j)
			mobileRobotPoseVisionLeft[i][j] = myVisionCurrent_AGV[j][0];
	for (int j = 0; j < odometerTimeDiff.size(); ++j)
	{
		int pointPre = 0;
		int pointPso = 0;
		bool preFlag = false;
		bool psoFlag = false;
		//先找到 里程计的时间 在 视觉时间序列中 的 前后两个最近 的位姿点
		for (auto it = visionTimeDiff.begin(); it != visionTimeDiff.end(); ++it)
		{
			if ((*it) >= odometerTimeDiff[j] && !psoFlag)
			{
				pointPso = it - visionTimeDiff.begin();
				psoFlag = true;
			}
			if (psoFlag)
				break;
		}
		for (auto it = visionTimeDiff.end() - 1; it != visionTimeDiff.begin(); --it)
		{
			if ((*it <= odometerTimeDiff[j]) && !preFlag)
			{
				pointPre = it - visionTimeDiff.begin();
				preFlag = true;
			}
			if (preFlag)
				break;
		}
		if (psoFlag)
			for (int i = 0; i < 4; ++i)
				mobileRobotPoseVisionLeft[j + odometerLeftStartPoint - 1][i] = (odometerTimeDiff[j] - visionTimeDiff[pointPso]) / (visionTimeDiff[pointPre] - visionTimeDiff[pointPso]) * myVisionCurrent_AGV[i][pointPre] + (odometerTimeDiff[j] - visionTimeDiff[pointPre]) / (visionTimeDiff[pointPso] - visionTimeDiff[pointPre]) * myVisionCurrent_AGV[i][pointPso];
	}

	//左天线：由上述插值得到的机器人实际路径，得到理论的相位序列
	vector<vector<vector<double>>> trackMobileRobotLeftAntennaReal(referenceTagNum, vector<vector<double>>(mymax(referenceTagCountLeft), vector<double>(2)));

	distanceTheory.clear();
	distanceTheory.resize(mymax(referenceTagCountLeft));
	for (auto& vec : distanceTheory)
		vec.resize(referenceTagNum);

	vector<vector<double>> phaseTagTheoryLeft(mymax(referenceTagCountLeft), vector<double>(referenceTagNum));
	for (int j = 0; j < referenceTagNum; ++j)
		if (referenceTagCountLeft[j] > 1)
			for (int i = 0; i < referenceTagCountLeft[j]; ++i)
			{
				trackMobileRobotLeftAntennaReal[j][i][0] = mobileRobotPoseVisionLeft[referenceTagNumLeft[i][j]][0] - antennaHLeftError * cos(PI - antennaAlphaLeftError - mobileRobotPoseVisionLeft[referenceTagNumLeft[i][j]][2]);
				trackMobileRobotLeftAntennaReal[j][i][1] = mobileRobotPoseVisionLeft[referenceTagNumLeft[i][j]][1] + antennaHLeftError * sin(PI - antennaAlphaLeftError - mobileRobotPoseVisionLeft[referenceTagNumLeft[i][j]][2]);
				distanceTheory[i][j] = sqrt(pow((trackMobileRobotLeftAntennaReal[j][i][0] - referenceTag[j][0]), 2) + pow((trackMobileRobotLeftAntennaReal[j][i][1] - referenceTag[j][1]), 2) + pow((myVisionLeftAntennaHigh - referenceTag[j][2]), 2)) * 2;
				phaseTagTheoryLeft[i][j] = distanceTheory[i][j] * 2 * PI / waveLengthVar[0];
			}

	//左天线：信号强度分析
	vector<vector<int>> RSSI_TagLeft(mymax(referenceTagCountLeft), vector<int>(referenceTagNum));
	for (int j = 0; j < referenceTagNum; ++j)
		if (referenceTagCountLeft[j] > 1)
			for (int i = 0; i < referenceTagCountLeft[j]; ++i)
				RSSI_TagLeft[i][j] = myDataResultLeft.data.RSSI[referenceTagNumLeft[i][j]];

	//执行粒子滤波定位算法
	int PF_Count = 500;

	//粒子状态（移动机器人位姿），第一、二、三行为机器人的x y th值
	vector<vector<vector<double>>> PF_ParticleRobot(500, vector<vector<double>>(3, vector<double>(PF_Count)));
	//粒子状态（变换矩阵），分别对应两个平移、一个旋转
	vector<vector<vector<double>>> PF_ParticleTransformation(500, vector<vector<double>>(3, vector<double>(PF_Count)));
	//粒子状态（左天线位置），第一、二行为左天线的x y值
	vector<vector<vector<double>>> PF_ParticleAntennaLeft(20000, vector<vector<double>>(2, vector<double>(PF_Count)));
	//粒子状态（右天线位置），第一、二行为右天线的x y值
	vector<vector<vector<double>>> PF_ParticleAntennaRight(20000, vector<vector<double>>(2, vector<double>(PF_Count)));

	//观测：保存左天线两个时刻之间的相位差值,以行排列
	vector<double> PF_ObserveGradientLeft(referenceTag.size());
	//观测：保存右天线两个时刻之间的相位差值,以行排列
	vector<double> PF_ObserveGradientRight(referenceTag.size());

	//预测：第1行是左天线前一个时刻点的相位（解缠相位），第2行是左天线当前时刻点的相位（解缠相位），第3为左天线两个时刻之间的相位差值
	vector<vector<vector<double>>> PF_PredictionLeft(referenceTag.size(), vector<vector<double>>(3, vector<double>(PF_Count)));
	//预测：第1行是右天线前一个时刻点的相位（解缠相位），第2行是右天线当前时刻点的相位（解缠相位），第3为右天线两个时刻之间的相位差值
	vector<vector<vector<double>>> PF_PredictionRight(referenceTag.size(), vector<vector<double>>(3, vector<double>(PF_Count)));

	//单位cm
	int gradientLen = 2;
	//单位cm
	int distanceFarThreshold = 200;
	//时间限制，单位s
	int gradientTimeLen = 10;
	//过程噪声，标准差，单位cm, rad
	vector<vector<double>> PF_Q = { {3, 0}, {0, 0.1} };
	//测量噪声，方差，单位rad
	double PF_R = (PhaseGauss + 0.4) * (PhaseGauss + 0.4);
	//粒子与观测值之间的差距,每一行对应一个参考标签
	vector<vector<double>> PF_DistanceLeft(referenceTag.size(), vector<double>(PF_Count));
	//粒子与观测值之间的差距,每一行对应一个参考标签
	vector<vector<double>> PF_DistanceRight(referenceTag.size(), vector<double>(PF_Count));
	//粒子与观测值之间的权重因子，第三行用于综合评价
	vector<vector<double>> PF_W_Left(referenceTag.size(), vector<double>(PF_Count));
	//粒子与观测值之间的权重因子，第三行用于综合评价
	vector<vector<double>> PF_W_Right(referenceTag.size(), vector<double>(PF_Count));
	//第一行：本次观测得权重；第二行：本次最终的权重
	vector<vector<vector<double>>> PF_W(Times, vector<vector<double>>(2, vector<double>(PF_Count)));
	//协方差矩阵
	vector<vector<vector<double>>> PF_CenterVar(200, vector<vector<double>>(3, vector<double>(3)));
	double neffRatioThreshold = 0.6;

	//标签可读性的相关判定值
	//不可读的半径，一定大于此值
	int tagUnreadableRadius = 20;
	//可读的半径，一定小于此值
	int tagReadableRadius = 110;

	//移动机器人本体的遮挡区域
	vector<vector<double>> leftShadowUnreadablePointCor = { {120, 20}, {-120, 20}, {-120, -120}, {120, -120} };

	vector<double> leftShadowUnreadablePointH(leftShadowUnreadablePointCor.size());
	for (int i = 0; i < leftShadowUnreadablePointH.size(); ++i)
		leftShadowUnreadablePointH[i] = sqrt(pow(leftShadowUnreadablePointCor[i][0], 2) + pow(leftShadowUnreadablePointCor[i][1], 2));

	vector<double> leftShadowUnreadablePointAlpha(leftShadowUnreadablePointCor.size());
	for (int i = 0; i < leftShadowUnreadablePointAlpha.size(); ++i)
		leftShadowUnreadablePointAlpha[i] = atan2(leftShadowUnreadablePointCor[i][1], leftShadowUnreadablePointCor[i][0]);

	vector<vector<double>> leftShadowReadablePointCor = { {120, -20}, {-120, -20}, {-120, -120}, {120, -120} };

	vector<double> leftShadowReadablePointH(leftShadowReadablePointCor.size());
	for (int i = 0; i < leftShadowReadablePointH.size(); ++i)
		leftShadowReadablePointH[i] = sqrt(pow(leftShadowReadablePointCor[i][0], 2) + pow(leftShadowReadablePointCor[i][1], 2));

	vector<double> leftShadowReadablePointAlpha(leftShadowReadablePointCor.size());
	for (int i = 0; i < leftShadowReadablePointAlpha.size(); ++i)
		leftShadowReadablePointAlpha[i] = atan2(leftShadowReadablePointCor[i][1], leftShadowReadablePointCor[i][0]);

	vector<vector<double>> rightShadowUnreadablePointCor = { {120, -20}, {-120, -20}, {-120, 120}, {120, 120} };

	vector<double> rightShadowUnreadablePointH(rightShadowUnreadablePointCor.size());
	for (int i = 0; i < rightShadowUnreadablePointH.size(); ++i)
		rightShadowUnreadablePointH[i] = sqrt(pow(rightShadowUnreadablePointCor[i][0], 2) + pow(rightShadowUnreadablePointCor[i][1], 2));

	vector<double> rightShadowUnreadablePointAlpha(rightShadowUnreadablePointCor.size());
	for (int i = 0; i < rightShadowUnreadablePointAlpha.size(); ++i)
		rightShadowUnreadablePointAlpha[i] = atan2(rightShadowUnreadablePointCor[i][1], rightShadowUnreadablePointCor[i][0]);

	vector<vector<double>> rightShadowReadablePointCor = { {120, 20}, {-120, 20}, {-120, 120}, {120, 120} };

	vector<double> rightShadowReadablePointH(rightShadowReadablePointCor.size());
	for (int i = 0; i < rightShadowReadablePointH.size(); ++i)
		rightShadowReadablePointH[i] = sqrt(pow(rightShadowReadablePointCor[i][0], 2) + pow(rightShadowReadablePointCor[i][1], 2));

	vector<double> rightShadowReadablePointAlpha(rightShadowReadablePointCor.size());
	for (int i = 0; i < rightShadowReadablePointAlpha.size(); ++i)
		rightShadowReadablePointAlpha[i] = atan2(rightShadowReadablePointCor[i][1], rightShadowReadablePointCor[i][0]);

	//在粒子初始化中，PF_tag_readable_count
	int PF_TagReadableCountInit = 5;
	//在粒子滤波更新中，标签可读性需要判断的周围的标签数量
	int PF_TagReadableCount = 5;
	//在粒子滤波中，数据段长度的下界
	int distanceIntervalLp = 1;

	//为了降低重采样的频率
	double recoveryAlphaSlow = 0.001; //0.05
	double recoveryAlphaFast = 0.2;	  //0.2
	int PF_W_Slow = 0;
	int PF_W_Fast = 0;

	//右天线的标号，为方便索引故减一
	int numberFlag = 0;
	//左天线的标号，为方便索引故减一
	int numberFlagVice = 0;

	vector<vector<double>> PF_ObserveNativeLeft = phasePsoLeft;
	vector<vector<double>> PF_ObserveNativeRight = phasePsoRight;

	vector<double> robotXtAssume(Times);
	vector<double> robotYtAssume(Times);
	vector<double> robotThtAssume(Times);
	vector<double> robotXt(Times);
	vector<double> robotYt(Times);
	vector<double> robotTht(Times);

	vector<vector<int>> numberFlagRight(2, vector<int>(referenceTagNum));
	vector<vector<int>> numberFlagLeft(2, vector<int>(referenceTagNum));
	vector<vector<double>> antennaLeft(2, vector<double>(Times));
	vector<vector<double>> antennaRight(2, vector<double>(Times));
	vector<vector<double>> robotEstimationError(5, vector<double>(Times));
	vector<int> readableTagRightNum(Times);
	vector<int> readableTagLeftNum(Times);

	vector<vector<double>> PF_CenterMean(4, vector<double>(Times));
	vector<vector<double>> robotPositionAssume(3, vector<double>(Times));
	vector<vector<int>> PF_W_RSSI(Times, vector<int>(PF_Count));
	vector<vector<vector<double>>> PF_ParticleAntennaRightPso;
	vector<vector<vector<double>>> PF_ParticleAntennaRightPre;
	vector<vector<vector<double>>> PF_ParticleAntennaLeftPso;
	vector<vector<vector<double>>> PF_ParticleAntennaLeftPre;
	vector<bool> phaseFlag(Times);
	vector<double> neffRatioPre(Times);
	vector<int> PF_ReSample(Times);
	vector<double> time(Times);

	/*验证数据导入：随机数*/
	vector<double> randNum(Times);
	vector<vector<vector<double>>> ran(Times, vector<vector<double>>(3, vector<double>(PF_Count)));
	MATFile* pM0 = matOpen("E:/MATLAB workspace/hgdw/test400.mat", "r");
	mxArray* ppV1 = matGetVariable(pM0, "ran");
	double* temp1 = (double*)mxGetData(ppV1);
	mxArray* ppV2 = matGetVariable(pM0, "randnum");
	double* temp2 = (double*)mxGetData(ppV2);
	for (int i = 0; i < Times; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			for (int k = 0; k < 500; ++k)
			{
				ran[i][j][k] = temp1[i * 500 * 3 + k * 3 + j];
			}
		}
		randNum[i] = temp2[i];
	}
	/*******检测******/

	for (int i = 0; i < Times; ++i)
	{
		auto startTime = clock();
		if (i == 0)
		{
			robotXtAssume[i] = trackMobileRobotRight[0][numberFlag];
			robotYtAssume[i] = trackMobileRobotRight[1][numberFlag];
			robotThtAssume[i] = trackMobileRobotRight[2][numberFlag];

			robotPositionAssume[0][i] = robotXtAssume[i];
			robotPositionAssume[1][i] = robotYtAssume[i];
			robotPositionAssume[2][i] = 1;

			//找到静止条件下，所能够读到的标签，并依据这些标签的坐标，直接对粒子进行初始化
			//右天线
			epcDataRight.clear();
			epcDataRight.resize(odometerRightStartPoint);
			for (int j = 0; j < odometerRightStartPoint; ++j)
				epcDataRight[j] = myDataResultRight.textData.epcData[j];

			//用完就扔
			vector<double> RSSI_DataRightCurrent(odometerRightStartPoint);
			for (int j = 0; j < odometerRightStartPoint; ++j)
				RSSI_DataRightCurrent[j] = myDataResultRight.data.RSSI[j];

			vector<double> referenceTagCountRightCurrent(referenceTagNum);
			for (int j = 0; j < odometerRightStartPoint; ++j)
				for (int k = 0; k < referenceTagNum; ++k)
					if ((referenceTag_EPC[k] == epcDataRight[j]) && (RSSI_DataRightCurrent[j] >= -53))
					{
						//数量加1
						++referenceTagCountRightCurrent[k];
						break;
					}

			//统计referenceTagCountRightCurrent中非零元素个数
			int count = 0;
			for (auto x : referenceTagCountRightCurrent)
				if (x)
					++count;
			readableTagRightNum[i] = count;

			//列数就是标签的位置
			vector<int> colRight;
			for (auto it = referenceTagCountRightCurrent.begin(); it != referenceTagCountRightCurrent.end(); ++it)
				if (*it)
					colRight.push_back(it - referenceTagCountRightCurrent.begin());

			//左天线
			epcDataLeft.clear();
			epcDataLeft.resize(odometerLeftStartPoint);
			for (int j = 0; j < odometerLeftStartPoint; ++j)
				epcDataLeft[j] = myDataResultLeft.textData.epcData[j];

			//用完就扔
			vector<double> RSSI_DataLeftCurrent(odometerLeftStartPoint);
			for (int j = 0; j < odometerLeftStartPoint; ++j)
				RSSI_DataLeftCurrent[j] = myDataResultLeft.data.RSSI[j];

			vector<double> referenceTagCountLeftCurrent(referenceTagNum);
			for (int j = 0; j < odometerLeftStartPoint; ++j)
				for (int k = 0; k < referenceTagNum; ++k)
					if ((referenceTag_EPC[k] == epcDataLeft[j]) && (RSSI_DataLeftCurrent[j] >= -53))
					{
						++referenceTagCountLeftCurrent[k];
						break;
					}

			//列数就是标签的位置
			vector<int> colLeft;
			for (auto it = referenceTagCountLeftCurrent.begin(); it != referenceTagCountLeftCurrent.end(); ++it)
				if (*it)
					colLeft.push_back(it - referenceTagCountLeftCurrent.begin());

			//稍微有点繁琐了
			vector<int> randData(4);
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

			//依据该范围，生成粒子
			//for (int j = 0; j < PF_Count; ++j)
			//{
			//	//粒子x坐标
			//	PF_ParticleRobot[i][0][j] = PF_Scope[0][0] + rand() / double(RAND_MAX) * (PF_Scope[0][1] - PF_Scope[0][0]);
			//	//粒子y坐标
			//	PF_ParticleRobot[i][1][j] = PF_Scope[1][0] + rand() / double(RAND_MAX) * (PF_Scope[1][1] - PF_Scope[1][0]);
			//	//粒子方向
			//	PF_ParticleRobot[i][2][j] = (-PI) + rand() / double(RAND_MAX) * (PI - (-PI));
			//}

			/*验证数据导入*/
			for (int j = 0; j < PF_Count; ++j)
			{
				//粒子x坐标
				PF_ParticleRobot[i][0][j] = PF_Scope[0][0] + (PF_Scope[0][1] - PF_Scope[0][0]) * ran[i][0][j];
				//粒子y坐标
				PF_ParticleRobot[i][1][j] = PF_Scope[1][0] + (PF_Scope[1][1] - PF_Scope[1][0]) * ran[i][1][j];
				//粒子方向
				PF_ParticleRobot[i][2][j] = (-PI) + (PI - (-PI)) * ran[i][2][j];
			}
			/********验证********/

			//机器人位姿所对应的天线的位置
			for (int j = 0; j < PF_Count; ++j)
			{
				PF_ParticleAntennaLeft[i][0][j] = PF_ParticleRobot[i][0][j] - antennaHLeftError * cos(PI - antennaAlphaLeftError - PF_ParticleRobot[i][2][j]);
				PF_ParticleAntennaLeft[i][1][j] = PF_ParticleRobot[i][1][j] + antennaHLeftError * sin(PI - antennaAlphaLeftError - PF_ParticleRobot[i][2][j]);
			}
			for (int j = 0; j < PF_Count; ++j)
			{
				PF_ParticleAntennaRight[i][0][j] = PF_ParticleRobot[i][0][j] + antennaHRightError * cos(-antennaAlphaRightError + PF_ParticleRobot[i][2][j]);
				PF_ParticleAntennaRight[i][1][j] = PF_ParticleRobot[i][1][j] + antennaHRightError * sin(-antennaAlphaRightError + PF_ParticleRobot[i][2][j]);
			}

			for (int j = 0; j < PF_Count; ++j)
				PF_W[i][1][j] = 1 / double(PF_Count);

			//x y坐标的均值及其协方差
			for (int j = 0; j < 2; ++j)
				for (int k = 0; k < PF_Count; ++k)
					PF_CenterMean[j][i] += PF_ParticleRobot[i][j][k] * PF_W[i][1][k];
			
			double weightSqr = 0;
			for (int j = 0; j < PF_Count; ++j)
				weightSqr += PF_W[i][1][j] * PF_W[i][1][j];
			double factor = 0;
			if (abs(weightSqr - 1.0) < sqrt(DBL_EPSILON))
				factor = 1.0;
			else
				factor = 1 / (1 - weightSqr);

			vector<vector<double>> meanDiff(2, vector<double>(PF_Count));
			for (int j = 0; j < 2; ++j)
				for (int k = 0; k < PF_Count; ++k)
					meanDiff[j][k] = PF_ParticleRobot[i][j][k] - PF_CenterMean[j][i];

			for (int j = 0; j < 2; ++j)
				for (int k = 0; k < 2; ++k)
				{
					for (int m = 0; m < PF_Count; ++m)
						PF_CenterVar[i][j][k] += meanDiff[j][m] * PF_W[i][1][m] * meanDiff[k][m];
					PF_CenterVar[i][j][k] *= factor;
				}

			//角度均值及协方差：
			double sinsum = 0;
			for (int j = 0; j < PF_Count; ++j)
				sinsum += PF_W[i][1][j] * sin(PF_ParticleRobot[i][2][j]);

			double cossum = 0;
			for (int j = 0; j < PF_Count; ++j)
				cossum += PF_W[i][1][j] * cos(PF_ParticleRobot[i][2][j]);

			double resultantLength = sqrt(pow(sinsum, 2) + pow(cossum, 2));
			PF_CenterMean[2][i] = atan2(sinsum, cossum);
			PF_CenterVar[i][2][2] = -2 * log(resultantLength);
			for (int j = 0; j < PF_Count; ++j)
				PF_CenterMean[3][i] += PF_ParticleRobot[i][2][j] * PF_W[i][1][j];

			//定位误差
			robotXt[i] = mobileRobotPoseVisionRight[numberFlag][0];
			robotYt[i] = mobileRobotPoseVisionRight[numberFlag][1];
			robotTht[i] = mobileRobotPoseVisionRight[numberFlag][2];

			antennaLeft[0][i] = robotXt[i] - antennaHLeftError * cos(PI - antennaAlphaLeftError - robotTht[i]);
			antennaLeft[1][i] = robotYt[i] + antennaHLeftError * sin(PI - antennaAlphaLeftError - robotTht[i]);
			antennaRight[0][i] = robotXt[i] + antennaHRightError * cos(-antennaAlphaRightError + robotTht[i]);
			antennaRight[1][i] = robotYt[i] + antennaHRightError * sin(-antennaAlphaRightError + robotTht[i]);

			robotEstimationError[0][i] = PF_CenterMean[0][i] - robotXt[i];
			robotEstimationError[1][i] = PF_CenterMean[1][i] - robotYt[i];
			robotEstimationError[2][i] = sqrt(pow(robotEstimationError[0][i], 2) + pow(robotEstimationError[1][i], 2));

			robotEstimationError[3][i] = PF_CenterMean[2][i] - robotTht[i];

			robotEstimationError[4][i] = PF_CenterMean[3][i] - robotTht[i];

			numberFlagRight[0].assign(referenceTagNum, -1);
			numberFlagRight[1].assign(referenceTagNum, -1);

			//绘图
			// 	 figure
			// h1 = plot(reference_tag(:,1), reference_tag(:,2), 'k.', 'markersize',25);  %系统状态位置
			// hold on
			// h2 = plot(robot_xt(1:i),robot_yt(1:i),'r.', 'markersize',50);%机器人的位置
			// h3 = plot(PF_particle_robot(1, :, i), PF_particle_robot(2, :, i),'.','markersize',10);
			// h4 = quiver(robot_xt(1:i),robot_yt(1:i), cos(robot_tht(i)), sin(robot_tht(i)),15, 'color', 'r', 'linewidth', 3);
			// h5 = quiver(PF_particle_robot(1, :, i), PF_particle_robot(2, :, i), cos(PF_particle_robot(3, :, i)), sin(PF_particle_robot(3, :, i)),0.3, 'color', 'c', 'linewidth', 3);
			// legend('参考标签阵列', '移动机器人真实位姿','移动机器人可能位姿');
			// set(gca,'child',[h2 h4 h3 h5 h1])
			// xlabel('x (cm)');ylabel('y (cm)');
			// legend('参考标签阵列', '移动机器人真实位姿','移动机器人可能位姿');

			// figure
			// h1 = plot(reference_tag(:,1), reference_tag(:,2), 'k.', 'marker', 'pentagram', 'markersize', 10);  %系统状态位置
			// hold on
			// h2 = plot(robot_xt(1:1),robot_yt(1:1),'r.', 'markersize',35);%机器人的位置
			// h3 = plot(PF_particle_robot(1, :, 1), PF_particle_robot(2, :, 1),'.','markersize',15);
			// h4 = quiver(robot_xt(1:1),robot_yt(1:1), cos(robot_tht(1)), sin(robot_tht(1)),3, 'color', 'r', 'linewidth', 3,'MaxHeadSize',5);
			// h5 = quiver(PF_particle_robot(1, :, 1), PF_particle_robot(2, :, 1), cos(PF_particle_robot(3, :, 1)), sin(PF_particle_robot(3, :, 1)),0.3, 'color', 'c', 'linewidth', 3);
			// legend('Reference tag', 'Ture pose of mobile robot','Distribution of particle');

			// set(gca,'child',[h2 h4 h3 h5 h1]);
			// xlabel('X (cm)');
			// ylabel('Y (cm)');
			// axis([-170 -110 -55 -00]);
		}
		else
		{
			//cout << i << endl;
			int numberFlagPre = numberFlag;

			while (1)
			{
				++numberFlag;
				if (numberFlag < myDataResultRight.data.readerTime.size())
				{
					double timeDiff = myDataResultRight.data.readerTime[numberFlag] - myDataResultRight.data.readerTime[numberFlagPre];
					if (timeDiff >= 200 * 1000)
						break;
				}
				else
					//跳出当前循环
					break;
			}
			//跳出大循环
			if (numberFlag >= myDataResultRight.data.readerTime.size())
				break;

			robotXtAssume[i] = trackMobileRobotRight[0][numberFlag];
			robotYtAssume[i] = trackMobileRobotRight[1][numberFlag];
			robotThtAssume[i] = trackMobileRobotRight[2][numberFlag];

			robotPositionAssume[0][i] = robotXtAssume[i];
			robotPositionAssume[1][i] = robotYtAssume[i];
			robotPositionAssume[2][i] = 1;

			//状态转移：预测
			double robotXtDiff = robotXtAssume[i] - robotXtAssume[i - 1];
			double robotYtDiff = robotYtAssume[i] - robotYtAssume[i - 1];
			double robotThtDiff = robotThtAssume[i] - robotThtAssume[i - 1];

			vector<double> positionDifference1(2);
			positionDifference1[0] = robotXtDiff * cos(-robotThtAssume[i - 1]) - robotYtDiff * sin(-robotThtAssume[i - 1]);
			positionDifference1[1] = robotXtDiff * sin(-robotThtAssume[i - 1]) + robotYtDiff * cos(-robotThtAssume[i - 1]);

			vector<vector<double>> positionDifference2(2, vector<double>(PF_Count));
			//拆开来写比较讲究。。。
			for (int j = 0; j < PF_Count; ++j)
			{
				positionDifference2[0][j] = positionDifference1[0] * cos(PF_ParticleRobot[i - 1][2][j]) - positionDifference1[1] * sin(PF_ParticleRobot[i - 1][2][j]);
				positionDifference2[1][j] = positionDifference1[0] * sin(PF_ParticleRobot[i - 1][2][j]) + positionDifference1[1] * cos(PF_ParticleRobot[i - 1][2][j]);
			}

			//for (int j = 0; j < PF_Count; ++j)
			//{
			//	PF_ParticleRobot[i][0][j] = PF_ParticleRobot[i - 1][0][j] + positionDifference2[0][j] + gaussrand(0, PF_Q[0][0]);
			//	PF_ParticleRobot[i][1][j] = PF_ParticleRobot[i - 1][1][j] + positionDifference2[1][j] + gaussrand(0, PF_Q[0][0]);
			//	PF_ParticleRobot[i][2][j] = PF_ParticleRobot[i - 1][2][j] + robotThtDiff + gaussrand(0, PF_Q[1][1]);
			//}

			/*验证数据导入*/
			for (int j = 0; j < PF_Count; ++j)
			{
				//粒子x坐标
				PF_ParticleRobot[i][0][j] = PF_ParticleRobot[i - 1][0][j] + positionDifference2[0][j] + ran[i][0][j];
				//粒子y坐标
				PF_ParticleRobot[i][1][j] = PF_ParticleRobot[i - 1][1][j] + positionDifference2[1][j] + ran[i][1][j];
				//粒子方向
				PF_ParticleRobot[i][2][j] = PF_ParticleRobot[i - 1][2][j] + robotThtDiff + ran[i][2][j];
			}
			/********验证********/

			//右天线
			epcDataRight.clear();
			epcDataRight.resize(numberFlag - numberFlagPre);
			for (int j = 0; j < (numberFlag - numberFlagPre); ++j)
				epcDataRight[j] = myDataResultRight.textData.epcData[numberFlagPre + 1 + j];

			//用完就扔
			vector<double> RSSI_DataRightCurrent(numberFlag - numberFlagPre);
			for (int j = 0; j < (numberFlag - numberFlagPre); ++j)
				RSSI_DataRightCurrent[j] = myDataResultRight.data.RSSI[numberFlagPre + 1 + j];

			vector<double> referenceTagCountRightCurrent(referenceTagNum);
			for (int j = 0; j < (numberFlag - numberFlagPre); ++j)
				for (int k = 0; k < referenceTagNum; ++k)
					if ((referenceTag_EPC[k] == epcDataRight[j]) && (RSSI_DataRightCurrent[j] >= -52))
					{
						//数量加1
						++referenceTagCountRightCurrent[k];
						break;
					}

			//统计referenceTagCountRightCurrent中非零元素个数
			int count = 0;
			for (auto x : referenceTagCountRightCurrent)
				if (x)
					++count;
			readableTagRightNum[i] = count;

			//左天线
			int numberFlagVicePre = numberFlagVice;

			//寻找满足条件的索引并附带上0，然后寻找最大值。
			int m;
			for (auto it = myDataResultLeft.data.readerTime.begin(); it != myDataResultLeft.data.readerTime.end(); ++it)
				if ((*it) < myDataResultRight.data.readerTime[numberFlag])
					m = it - myDataResultLeft.data.readerTime.begin();
			numberFlagVice = (m > 0) ? m : 0;

			//左天线
			epcDataLeft.clear();
			epcDataLeft.resize(numberFlagVice - numberFlagVicePre);
			for (int j = 0; j < (numberFlagVice - numberFlagVicePre); ++j)
				epcDataLeft[j] = myDataResultLeft.textData.epcData[numberFlagVicePre + 1 + j];

			//用完就扔
			vector<double> RSSI_DataLeftCurrent(numberFlagVice - numberFlagVicePre);
			for (int j = 0; j < (numberFlagVice - numberFlagVicePre); ++j)
				RSSI_DataLeftCurrent[j] = myDataResultLeft.data.RSSI[numberFlagVicePre + 1 + j];

			vector<double> referenceTagCountLeftCurrent(referenceTagNum);
			for (int j = 0; j < (numberFlagVice - numberFlagVicePre); ++j)
				for (int k = 0; k < referenceTagNum; ++k)
					if ((referenceTag_EPC[k] == epcDataLeft[j]) && (RSSI_DataLeftCurrent[j] >= -52))
					{
						//数量加1
						++referenceTagCountLeftCurrent[k];
						break;
					}

			//统计referenceTagCountRightCurrent中非零元素个数
			count = 0;
			for (auto x : referenceTagCountLeftCurrent)
				if (x)
					++count;
			readableTagLeftNum[i] = count;

			//利用新标签的可读性，对粒子进行初步筛选--不符合可读性的，权重变为0
			//寻找指定范围内的最大值以及索引
			//最大值
			int RSSI_LeftMaxIndex = myDataResultLeft.data.RSSI[numberFlagVicePre + 1];
			//最大值在指定范围内的索引
			int I_Left = 0;
			for (int j = numberFlagVicePre + 1; j < numberFlagVice; ++j)
				if (myDataResultLeft.data.RSSI[j] > RSSI_LeftMaxIndex)
				{
					RSSI_LeftMaxIndex = myDataResultLeft.data.RSSI[j];
					I_Left = j - numberFlagVicePre - 1;
				}

			//寻找相等并记录行与列的信息
			int rowLeft = 0;
			int colLeft = 0;
			for (int j = 0; j < referenceTagNumLeft.size(); ++j)
				for (int k = 0; k < referenceTagNumLeft[0].size(); ++k)
					if (referenceTagNumLeft[j][k] == (I_Left + numberFlagVicePre + 1))
					{
						rowLeft = j;
						colLeft = k;
					}

			//实现k-临近算法 简化为一个函数
			//找到离RSSI最大的点的最近的9个参考标签
			vector<int> optionalTagLeftFlag = knnsearch(referenceTag, referenceTag[colLeft], PF_TagReadableCount);

			//找出右天线中探寻的若干个点
			//寻找指定范围内的最大值以及索引
			//最大值
			int RSSI_RightMaxIndex = myDataResultRight.data.RSSI[numberFlagPre + 1];
			//最大值在指定范围内的索引
			int I_Right = 0;
			for (int j = numberFlagPre + 1; j < numberFlag; ++j)
				if (myDataResultRight.data.RSSI[j] > RSSI_RightMaxIndex)
				{
					RSSI_RightMaxIndex = myDataResultRight.data.RSSI[j];
					I_Right = j - numberFlagPre - 1;
				}

			//寻找相等并记录行与列的信息
			int rowRight = 0;
			int colRight = 0;
			for (int j = 0; j < referenceTagNumRight.size(); ++j)
				for (int k = 0; k < referenceTagNumRight[0].size(); ++k)
					if (referenceTagNumRight[j][k] == (I_Right + numberFlagPre + 1))
					{
						rowRight = j;
						colRight = k;
					}

			//找到离RSSI最大的点的最近的9个参考标签
			vector<int> optionalTagRightFlag = knnsearch(referenceTag, referenceTag[colRight], PF_TagReadableCount);

			//所探寻的9个标签的可读性
			vector<bool> leftTagReadFlag(PF_TagReadableCount);
			for (int j = 0; j < PF_TagReadableCount; ++j)
				leftTagReadFlag[j] = (referenceTagCountLeft[optionalTagLeftFlag[j]] > 1) && (referenceTagCountLeftCurrent[optionalTagLeftFlag[j]] > 0);
			vector<bool> rightTagReadFlag(PF_TagReadableCount);
			for (int j = 0; j < PF_TagReadableCount; ++j)
				rightTagReadFlag[j] = (referenceTagCountRight[optionalTagRightFlag[j]] > 1) && (referenceTagCountRightCurrent[optionalTagRightFlag[j]] > 0);

			//计算当前位姿下，左右天线和阴影控制点位置
			vector<vector<double>> PF_ParticleAntennaLeftCurrent(2, vector<double>(PF_Count));
			vector<vector<double>> PF_ParticleAntennaRightCurrent(2, vector<double>(PF_Count));
			for (int j = 0; j < PF_Count; ++j)
			{
				PF_ParticleAntennaLeftCurrent[0][j] = PF_ParticleRobot[i][0][j] - antennaHLeftError * cos(PI - antennaAlphaLeftError - PF_ParticleRobot[i][2][j]);
				PF_ParticleAntennaLeftCurrent[1][j] = PF_ParticleRobot[i][1][j] + antennaHLeftError * sin(PI - antennaAlphaLeftError - PF_ParticleRobot[i][2][j]);
				PF_ParticleAntennaRightCurrent[0][j] = PF_ParticleRobot[i][0][j] + antennaHRightError * cos(-antennaAlphaRightError + PF_ParticleRobot[i][2][j]);
				PF_ParticleAntennaRightCurrent[1][j] = PF_ParticleRobot[i][1][j] + antennaHRightError * sin(-antennaAlphaRightError + PF_ParticleRobot[i][2][j]);
			}

			vector<vector<double>> PF_ParticleShadowLeftCurrent(4, vector<double>(4));
			vector<vector<double>> PF_ParticleShadowRightCurrent(4, vector<double>(4));
			for (int j = 0; j < PF_Count; ++j)
			{
				//cout << j << ' ';
				for (int k = 0; k < 4; ++k)
				{
					PF_ParticleShadowLeftCurrent[k][0] = PF_ParticleRobot[i][0][j] + leftShadowUnreadablePointH[k] * cos(leftShadowUnreadablePointAlpha[k] + PF_ParticleRobot[i][2][j]);
					PF_ParticleShadowLeftCurrent[k][1] = PF_ParticleRobot[i][1][j] + leftShadowUnreadablePointH[k] * sin(leftShadowUnreadablePointAlpha[k] + PF_ParticleRobot[i][2][j]);
					PF_ParticleShadowLeftCurrent[k][2] = PF_ParticleRobot[i][0][j] + leftShadowReadablePointH[k] * cos(leftShadowReadablePointAlpha[k] + PF_ParticleRobot[i][2][j]);
					PF_ParticleShadowLeftCurrent[k][3] = PF_ParticleRobot[i][1][j] + leftShadowReadablePointH[k] * sin(leftShadowReadablePointAlpha[k] + PF_ParticleRobot[i][2][j]);

					PF_ParticleShadowRightCurrent[k][0] = PF_ParticleRobot[i][0][j] + rightShadowUnreadablePointH[k] * cos(rightShadowUnreadablePointAlpha[k] + PF_ParticleRobot[i][2][j]);
					PF_ParticleShadowRightCurrent[k][1] = PF_ParticleRobot[i][1][j] + rightShadowUnreadablePointH[k] * sin(rightShadowUnreadablePointAlpha[k] + PF_ParticleRobot[i][2][j]);
					PF_ParticleShadowRightCurrent[k][2] = PF_ParticleRobot[i][0][j] + rightShadowReadablePointH[k] * cos(rightShadowReadablePointAlpha[k] + PF_ParticleRobot[i][2][j]);
					PF_ParticleShadowRightCurrent[k][3] = PF_ParticleRobot[i][1][j] + rightShadowReadablePointH[k] * sin(rightShadowReadablePointAlpha[k] + PF_ParticleRobot[i][2][j]);
				}

				//逐个标签确定是否当前的移动机器人位姿是否符合要求
				vector<bool> optionalFlag(PF_TagReadableCount);

				//参考标签与天线的距离
				vector<double> leftDistanceThresholdFlag(PF_TagReadableCount);
				vector<double> rightDistanceThresholdFlag(PF_TagReadableCount);
				for (int k = 0; k < PF_TagReadableCount; ++k)
				{
					leftDistanceThresholdFlag[k] = sqrt(pow((referenceTag[optionalTagLeftFlag[k]][0] - PF_ParticleAntennaLeftCurrent[0][j]), 2) + pow((referenceTag[optionalTagLeftFlag[k]][1] - PF_ParticleAntennaLeftCurrent[1][j]), 2));
					rightDistanceThresholdFlag[k] = sqrt(pow((referenceTag[optionalTagRightFlag[k]][0] - PF_ParticleAntennaRightCurrent[0][j]), 2) + pow((referenceTag[optionalTagRightFlag[k]][1] - PF_ParticleAntennaRightCurrent[1][j]), 2));
				}

				//参考标签与阴影区的关系
				//参考标签是否在阴影区内
				vector<double> xPointLeft(PF_TagReadableCount);
				vector<double> yPointLeft(PF_TagReadableCount);
				for (int k = 0; k < PF_TagReadableCount; ++k)
				{
					xPointLeft[k] = referenceTag[optionalTagLeftFlag[k]][0];
					yPointLeft[k] = referenceTag[optionalTagLeftFlag[k]][1];
				}
				vector<double> xLineLeft(4);
				vector<double> yLineLeft(4);
				for (int k = 0; k < 4; ++k)
				{
					xLineLeft[k] = PF_ParticleShadowLeftCurrent[k][2];
					yLineLeft[k] = PF_ParticleShadowLeftCurrent[k][3];
				}
				vector<bool> leftInpolygonReadableFlag = inpolygon(xPointLeft, yPointLeft, xLineLeft, yLineLeft);

				//参考标签是否在阴影区内
				vector<double> xPointRight(PF_TagReadableCount);
				vector<double> yPointRight(PF_TagReadableCount);
				for (int k = 0; k < PF_TagReadableCount; ++k)
				{
					xPointRight[k] = referenceTag[optionalTagRightFlag[k]][0];
					yPointRight[k] = referenceTag[optionalTagRightFlag[k]][1];
				}
				vector<double> xLineRight(4);
				vector<double> yLineRight(4);
				for (int k = 0; k < 4; ++k)
				{
					xLineRight[k] = PF_ParticleShadowRightCurrent[k][2];
					yLineRight[k] = PF_ParticleShadowRightCurrent[k][3];
				}
				vector<bool> rightInpolygonReadableFlag = inpolygon(xPointRight, yPointRight, xLineRight, yLineRight);

				for (int k = 0; k < PF_TagReadableCount; ++k)
					//左右天线都未读到
					if ((leftTagReadFlag[k] == 0) && (rightTagReadFlag[k] == 0))
					{
						//距离大于阈值，或者 在阴影区内|| left_inpolygon_unreadable_flag(k) == 1
						if (leftDistanceThresholdFlag[k] > tagUnreadableRadius)
							//距离大于阈值，或者 在阴影区内|| right_inpolygon_unreadable_flag(k) == 1
							if (rightDistanceThresholdFlag[k] > tagUnreadableRadius)
								optionalFlag[k] = 1;
					}
				//左天线读到，右天线未读到
					else if ((leftTagReadFlag[k] == 1) && (rightTagReadFlag[k] == 0))
					{
						//距离小于阈值，并且 不在阴影区内
						if ((leftDistanceThresholdFlag[k] < tagReadableRadius) && (leftInpolygonReadableFlag[k] == 0))
							//距离大于阈值，或者 在阴影区内|| right_inpolygon_unreadable_flag(k) == 1
							if (rightDistanceThresholdFlag[k] > tagUnreadableRadius)
								optionalFlag[k] = 1;
					}
				//左天线未读到，右天线读到
					else if ((leftTagReadFlag[k] == 0) && (rightTagReadFlag[k] == 1))
					{
						//距离大于阈值，或者 在阴影区内|| left_inpolygon_unreadable_flag(k) == 1
						if (leftDistanceThresholdFlag[k] > tagUnreadableRadius)
							//距离小于阈值，并且 不在阴影区内
							if ((rightDistanceThresholdFlag[k] < tagReadableRadius) && (rightInpolygonReadableFlag[k] == 0))
								optionalFlag[k] = 1;
					}
				//左右天线都读到
					else if ((leftTagReadFlag[k] == 1) && (rightTagReadFlag[k] == 1))
					{
						//距离小于阈值，并且 不再阴影区内
						if ((leftDistanceThresholdFlag[k] < tagReadableRadius) && (leftInpolygonReadableFlag[k] == 0))
							//距离小于阈值，并且 不再阴影区内
							if ((rightDistanceThresholdFlag[k] < tagReadableRadius) && (rightInpolygonReadableFlag[k] == 0))
								optionalFlag[k] = 1;
					}

				//所有标签的可读性都符合常理
				int count = 0;
				for (auto flag : optionalFlag)
					if (flag)
						++count;
				if (count == PF_TagReadableCount)
					PF_W_RSSI[i][j] = 1;
				else
					PF_W_RSSI[i][j] = 0;
				
			}

			//右天线
			auto numberFlagRightPre = numberFlagRight;
			for (int j = 0; j < referenceTagNum; ++j)
			{
				if (referenceTagCountRightCurrent[j] >= 1)
				{
					//从后向前找第一个小于numberFlag的值的索引
					int numberFlagOption = 0;
					bool f = false;
					for (int k = referenceTagCountRight[j] - 1; k >= 0; --k)
					{
						if ((referenceTagNumRight[k][j] < numberFlag) && !f)
						{
							numberFlagOption = k;
							f = true;
						}
						if (f)
							break;
					}
					if (!isnan(double(numberFlagOption)))
						//因为是C++里的索引，所以不是==1而是==0
						if (numberFlagOption == 0)
						{
							//同理，预判为可能索引的存在，故赋值为0
							numberFlagRight[1][j] = 0;
							numberFlagRight[0][j] = 0;
						}
						else
						{
							numberFlagRight[1][j] = numberFlagOption;
							numberFlagRight[0][j] = numberFlagRight[1][j];
							double distanceInterval = sqrt(pow((trackMobileRobotRight[0][referenceTagNumRight[numberFlagRight[1][j]][j]] - trackMobileRobotRight[0][numberFlag]), 2) +
								pow((trackMobileRobotRight[1][referenceTagNumRight[numberFlagRight[1][j]][j]] - trackMobileRobotRight[1][numberFlag]), 2));
							if (distanceInterval < distanceFarThreshold)
								while (1)
								{
									distanceInterval = sqrt(pow((trackMobileRobotRight[0][referenceTagNumRight[numberFlagRight[1][j]][j]] - trackMobileRobotRight[0][referenceTagNumRight[numberFlagRight[0][j]][j]]), 2) +
										pow((trackMobileRobotRight[1][referenceTagNumRight[numberFlagRight[1][j]][j]] - trackMobileRobotRight[1][referenceTagNumRight[numberFlagRight[0][j]][j]]), 2));
									if (distanceInterval > gradientLen)
										break;
									else
										--numberFlagRight[0][j];
									if (numberFlagRight[0][j] < 1)
										break;
									double rightPhaseDis = sqrt(pow((trackMobileRobotRight[0][referenceTagNumRight[numberFlagRight[0][j]][j]] - trackMobileRobotRight[0][referenceTagNumRight[numberFlagRight[0][j] + 1][j]]), 2) +
										pow((trackMobileRobotRight[1][referenceTagNumRight[numberFlagRight[0][j]][j]] - trackMobileRobotRight[1][referenceTagNumRight[numberFlagRight[0][j] + 1][j]]), 2));
									if (rightPhaseDis > 7)
										if (distanceInterval > distanceIntervalLp)
										{
											++numberFlagRight[0][j];
											break;
										}
										else
										{
											numberFlagRight[1][j] = 0;
											numberFlagRight[0][j] = 0;
											break;
										}
								}
							else
							{
								numberFlagRight[1][j] = 0;
								numberFlagRight[0][j] = 0;
							}
						}
					else
					{
						numberFlagRight[1][j] = 0;
						numberFlagRight[0][j] = 0;
					}
				}
				else
				{
					numberFlagRight[1][j] = 0;
					numberFlagRight[0][j] = 0;
				}
			}

			//左天线：找主天线下，到每个标签对应的数据段--左天线。注：超过一定时间的数据不用，该逻辑还未实现
			for (int j = 0; j < referenceTagNum; ++j)
			{
				if ((referenceTagCountLeftCurrent[j] >= 1) && (numberFlagVice > 1))
				{
					//从后向前找第一个小于numberFlag的值的索引
					int numberFlagOption = 0;
					bool f = false;
					for (int k = referenceTagCountLeft[j] - 1; k >= 0; --k)
					{
						if ((referenceTagNumLeft[k][j] < numberFlagVice) && !f)
						{
							numberFlagOption = k;
							f = true;
						}
						if (f)
							break;
					}
					if (!isnan(double(numberFlagOption)))
						//因为是C++里的索引，所以不是==1而是==0
						if (numberFlagOption == 0)
						{
							//同理，预判为可能索引的存在，故赋值为0
							numberFlagLeft[1][j] = 0;
							numberFlagLeft[0][j] = 0;
						}
						else
						{
							numberFlagLeft[1][j] = numberFlagOption;
							numberFlagLeft[0][j] = numberFlagLeft[1][j];
							//所以为什么是left和right相减呢
							double distanceInterval = sqrt(pow((trackMobileRobotLeft[0][referenceTagNumLeft[numberFlagLeft[1][j]][j]] - trackMobileRobotRight[0][numberFlag]), 2) +
								pow((trackMobileRobotLeft[1][referenceTagNumLeft[numberFlagLeft[1][j]][j]] - trackMobileRobotRight[1][numberFlag]), 2));
							if (distanceInterval < distanceFarThreshold)
								while (1)
								{
									distanceInterval = sqrt(pow((trackMobileRobotLeft[0][referenceTagNumLeft[numberFlagLeft[1][j]][j]] - trackMobileRobotLeft[0][referenceTagNumLeft[numberFlagLeft[0][j]][j]]), 2) +
										pow((trackMobileRobotLeft[1][referenceTagNumLeft[numberFlagLeft[1][j]][j]] - trackMobileRobotLeft[1][referenceTagNumLeft[numberFlagLeft[0][j]][j]]), 2));
									double leftTimeInterval = (readerTimeLeft[referenceTagNumLeft[numberFlagLeft[1][j]][j]] - readerTimeLeft[referenceTagNumLeft[numberFlagLeft[0][j]][j]]) / 1000000;
									//距离间隔不能大于...，时间间隔不能大于...
									if ((distanceInterval > gradientLen) || (leftTimeInterval > gradientTimeLen))
										break;
									else
										--numberFlagLeft[0][j];
									if (numberFlagLeft[0][j] < 1)
										break;
									double leftPhaseDis = sqrt(pow((trackMobileRobotLeft[0][referenceTagNumLeft[numberFlagLeft[0][j]][j]] - trackMobileRobotLeft[0][referenceTagNumLeft[numberFlagLeft[0][j] + 1][j]]), 2) +
										pow((trackMobileRobotLeft[1][referenceTagNumLeft[numberFlagLeft[0][j]][j]] - trackMobileRobotLeft[1][referenceTagNumLeft[numberFlagLeft[0][j] + 1][j]]), 2));
									if (leftPhaseDis > 7)
										if (distanceInterval > distanceIntervalLp)
										{
											++numberFlagLeft[0][j];
											break;
										}
										else
										{
											numberFlagLeft[1][j] = 0;
											numberFlagLeft[0][j] = 0;
											break;
										}
								}
							else
							{
								numberFlagLeft[1][j] = 0;
								numberFlagLeft[0][j] = 0;
							}
						}
					else
					{
						numberFlagLeft[1][j] = 0;
						numberFlagLeft[0][j] = 0;
					}
				}
				else
				{
					numberFlagLeft[1][j] = 0;
					numberFlagLeft[0][j] = 0;
				}
			}

			//右天线
			vector<double> RSSI_MeanTagRight(referenceTagNum, -200);
			for (int j = 0; j < referenceTagNum; ++j)
				if (numberFlagRight[1][j] != numberFlagRight[0][j])
				{
					double sum = 0;
					for (int k = numberFlagRight[0][j]; k <= numberFlagRight[1][j]; ++k)
						sum += RSSI_TagRight[k][j];
					RSSI_MeanTagRight[j] = sum / (numberFlagRight[1][j] - numberFlagRight[0][j] + 1);
				}

			//左天线
			vector<double> RSSI_MeanTagLeft(referenceTagNum, -200);
			for (int j = 0; j < referenceTagNum; ++j)
				if (numberFlagLeft[1][j] != numberFlagLeft[0][j])
				{
					double sum = 0;
					for (int k = numberFlagLeft[0][j]; k <= numberFlagLeft[1][j]; ++k)
						sum += RSSI_TagLeft[k][j];
					RSSI_MeanTagLeft[j] = sum / (numberFlagLeft[1][j] - numberFlagLeft[0][j] + 1);
				}

			//升序排序
			vector<int> I_RSSI_Right(referenceTagNum);
			vector<double> RMTR(RSSI_MeanTagRight.begin(), RSSI_MeanTagRight.end());
			for (int j = 0; j < referenceTagNum; ++j)
				I_RSSI_Right[j] = j;
			for (int j = 0; j < referenceTagNum; ++j)
			{
				bool f = false;
				for (int k = referenceTagNum - 1; k > j; --k)
					if (RMTR[k] < RMTR[k - 1])
					{
						int temp = I_RSSI_Right[k];
						I_RSSI_Right[k] = I_RSSI_Right[k - 1];
						I_RSSI_Right[k - 1] = temp;
						double tem = RMTR[k];
						RMTR[k] = RMTR[k - 1];
						RMTR[k - 1] = tem;
						f = true;
					}
				if (!f)
					break;
			}

			//右天线：计算每个标签的前后两个时刻点，移动机器人的位姿，及其对应的天线位置, 并计算其权值
			vector<double> PF_W_RightImprovement(PF_Count, -1);

			//选出最合适的标签
			if (RSSI_MeanTagRight[I_RSSI_Right.back()] > -53)
			{
				int j = I_RSSI_Right.back();
				if (numberFlagRight[1][j] != numberFlagRight[0][j])
				{
					//后一个时刻，粒子状态--右天线
					double robotXtAssumePso = trackMobileRobotRight[0][referenceTagNumRight[numberFlagRight[1][j]][j]];
					double robotYtAssumePso = trackMobileRobotRight[1][referenceTagNumRight[numberFlagRight[1][j]][j]];
					double robotThtAssumePso = trackMobileRobotRight[2][referenceTagNumRight[numberFlagRight[1][j]][j]];

					robotXtDiff = robotXtAssumePso - robotXtAssume[i];
					robotYtDiff = robotYtAssumePso - robotYtAssume[i];
					robotThtDiff = robotThtAssumePso - robotThtAssume[i];

					positionDifference1[0] = robotXtDiff * cos(-robotThtAssume[i]) - robotYtDiff * sin(-robotThtAssume[i]);
					positionDifference1[1] = robotXtDiff * sin(-robotThtAssume[i]) + robotYtDiff * cos(-robotThtAssume[i]);

					for (int k = 0; k < PF_Count; ++k)
					{
						positionDifference2[0][k] = positionDifference1[0] * cos(PF_ParticleRobot[i][2][k]) - positionDifference1[1] * sin(PF_ParticleRobot[i][2][k]);
						positionDifference2[1][k] = positionDifference1[0] * sin(PF_ParticleRobot[i][2][k]) + positionDifference1[1] * cos(PF_ParticleRobot[i][2][k]);
					}

					vector<vector<double>> PF_ParticleRobotRightPso(3, vector<double>(PF_Count));
					for (int k = 0; k < PF_Count; ++k)
					{
						PF_ParticleRobotRightPso[0][k] = PF_ParticleRobot[i][0][k] + positionDifference2[0][k];
						PF_ParticleRobotRightPso[1][k] = PF_ParticleRobot[i][1][k] + positionDifference2[1][k];
						PF_ParticleRobotRightPso[2][k] = PF_ParticleRobot[i][2][k] + robotThtDiff;
					}

					if (PF_ParticleAntennaRightPso.size() < (j + 1))
						PF_ParticleAntennaRightPso.resize(j + 1);
					vector<vector<double>> temp(2, vector<double>(PF_Count));
					for (int k = 0; k < PF_Count; ++k)
					{
						//x坐标
						temp[0][k] = PF_ParticleRobotRightPso[0][k] + antennaHRightError * cos(-antennaAlphaRightError + PF_ParticleRobotRightPso[2][k]);
						//y坐标
						temp[1][k] = PF_ParticleRobotRightPso[1][k] + antennaHRightError * sin(-antennaAlphaRightError + PF_ParticleRobotRightPso[2][k]);
					}
					PF_ParticleAntennaRightPso[j] = temp;

					//前一个时刻，粒子状态--右天线
					double robotXtAssumePre = trackMobileRobotRight[0][referenceTagNumRight[numberFlagRight[0][j]][j]];
					double robotYtAssumePre = trackMobileRobotRight[1][referenceTagNumRight[numberFlagRight[0][j]][j]];
					double robotThtAssumePre = trackMobileRobotRight[2][referenceTagNumRight[numberFlagRight[0][j]][j]];

					robotXtDiff = robotXtAssumePre - robotXtAssume[i];
					robotYtDiff = robotYtAssumePre - robotYtAssume[i];
					robotThtDiff = robotThtAssumePre - robotThtAssume[i];

					positionDifference1[0] = robotXtDiff * cos(-robotThtAssume[i]) - robotYtDiff * sin(-robotThtAssume[i]);
					positionDifference1[1] = robotXtDiff * sin(-robotThtAssume[i]) + robotYtDiff * cos(-robotThtAssume[i]);

					for (int k = 0; k < PF_Count; ++k)
					{
						positionDifference2[0][k] = positionDifference1[0] * cos(PF_ParticleRobot[i][2][k]) - positionDifference1[1] * sin(PF_ParticleRobot[i][2][k]);
						positionDifference2[1][k] = positionDifference1[0] * sin(PF_ParticleRobot[i][2][k]) + positionDifference1[1] * cos(PF_ParticleRobot[i][2][k]);
					}

					vector<vector<double>> PF_ParticleRobotRightPre(3, vector<double>(PF_Count));
					for (int k = 0; k < PF_Count; ++k)
					{
						PF_ParticleRobotRightPre[0][k] = PF_ParticleRobot[i][0][k] + positionDifference2[0][k];
						PF_ParticleRobotRightPre[1][k] = PF_ParticleRobot[i][1][k] + positionDifference2[1][k];
						PF_ParticleRobotRightPre[2][k] = PF_ParticleRobot[i][2][k] + robotThtDiff;
					}

					if (PF_ParticleAntennaRightPre.size() < (j + 1))
						PF_ParticleAntennaRightPre.resize(j + 1);
					for (int k = 0; k < PF_Count; ++k)
					{
						//x坐标
						temp[0][k] = PF_ParticleRobotRightPre[0][k] + antennaHRightError * cos(-antennaAlphaRightError + PF_ParticleRobotRightPre[2][k]);
						//y坐标
						temp[1][k] = PF_ParticleRobotRightPre[1][k] + antennaHRightError * sin(-antennaAlphaRightError + PF_ParticleRobotRightPre[2][k]);
					}
					PF_ParticleAntennaRightPre[j] = temp;

					//观测相位梯度
					PF_ObserveGradientRight[j] = PF_ObserveNativeRight[numberFlagRight[1][j]][j] - PF_ObserveNativeRight[numberFlagRight[0][j]][j];

					//理论相位梯度
					vector<double> PF_DistanceRightAntennaPre(PF_Count);
					vector<double> PF_DistanceRightAntennaPso(PF_Count);
					for (int k = 0; k < PF_Count; ++k)
					{
						PF_DistanceRightAntennaPre[k] = sqrt(pow((referenceTag[j][0] - PF_ParticleAntennaRightPre[j][0][k]), 2) + pow((referenceTag[j][1] - PF_ParticleAntennaRightPre[j][1][k]), 2) + pow((referenceTag[j][2] - myVisionRightAntennaHigh), 2)) * 2;
						PF_DistanceRightAntennaPso[k] = sqrt(pow((referenceTag[j][0] - PF_ParticleAntennaRightPso[j][0][k]), 2) + pow((referenceTag[j][1] - PF_ParticleAntennaRightPso[j][1][k]), 2) + pow((referenceTag[j][2] - myVisionRightAntennaHigh), 2)) * 2;
						PF_PredictionRight[j][2][k] = (PF_DistanceRightAntennaPso[k] - PF_DistanceRightAntennaPre[k]) * 2 * PI / waveLengthVar[0];
					}

					//粒子权重评估：用相位差
					for (int k = 0; k < PF_Count; ++k)
					{
						PF_DistanceRight[j][k] = PF_PredictionRight[j][2][k] - PF_ObserveGradientRight[j];
						//求权重
						PF_W_Right[j][k] = (1 / sqrt(2 * PF_R) / sqrt(2 * PI)) * exp(-pow(PF_DistanceRight[j][k], 2) / 2 / (2 * PF_R));
						PF_W_RightImprovement[k] = PF_W_Right[j][k];
					}
				}
				else
					for (auto& x : PF_W_Right[j])
						x = -1;
			}

			//升序排序
			vector<int> I_RSSI_Left(referenceTagNum);
			vector<double> RMTL(RSSI_MeanTagLeft.begin(), RSSI_MeanTagLeft.end());
			for (int j = 0; j < referenceTagNum; ++j)
				I_RSSI_Left[j] = j;
			for (int j = 0; j < referenceTagNum; ++j)
			{
				bool f = false;
				for (int k = referenceTagNum - 1; k > j; --k)
					if (RMTL[k] < RMTL[k - 1])
					{
						int temp = I_RSSI_Left[k];
						I_RSSI_Left[k] = I_RSSI_Left[k - 1];
						I_RSSI_Left[k - 1] = temp;
						double tem = RMTL[k];
						RMTL[k] = RMTL[k - 1];
						RMTL[k - 1] = tem;
						f = true;
					}
				if (!f)
					break;
			}

			//左天线：计算每个标签的前后两个时刻点，移动机器人的位姿，及其对应的天线位置, 并计算其权值
			vector<double> PF_W_LeftImprovement(PF_Count, -1);

			if (RSSI_MeanTagLeft[I_RSSI_Left.back()] > -53)
			{
				int j = I_RSSI_Left.back();
				if (numberFlagLeft[1][j] != numberFlagLeft[0][j])
				{
					//后一个时刻，粒子状态--左天线
					double robotXtAssumePso = trackMobileRobotLeft[0][referenceTagNumLeft[numberFlagLeft[1][j]][j]];
					double robotYtAssumePso = trackMobileRobotLeft[1][referenceTagNumLeft[numberFlagLeft[1][j]][j]];
					double robotThtAssumePso = trackMobileRobotLeft[2][referenceTagNumLeft[numberFlagLeft[1][j]][j]];

					robotXtDiff = robotXtAssumePso - robotXtAssume[i];
					robotYtDiff = robotYtAssumePso - robotYtAssume[i];
					robotThtDiff = robotThtAssumePso - robotThtAssume[i];

					positionDifference1[0] = robotXtDiff * cos(-robotThtAssume[i]) - robotYtDiff * sin(-robotThtAssume[i]);
					positionDifference1[1] = robotXtDiff * sin(-robotThtAssume[i]) + robotYtDiff * cos(-robotThtAssume[i]);

					for (int k = 0; k < PF_Count; ++k)
					{
						positionDifference2[0][k] = positionDifference1[0] * cos(PF_ParticleRobot[i][2][k]) - positionDifference1[1] * sin(PF_ParticleRobot[i][2][k]);
						positionDifference2[1][k] = positionDifference1[0] * sin(PF_ParticleRobot[i][2][k]) + positionDifference1[1] * cos(PF_ParticleRobot[i][2][k]);
					}

					vector<vector<double>> PF_ParticleRobotLeftPso(3, vector<double>(PF_Count));
					for (int k = 0; k < PF_Count; ++k)
					{
						PF_ParticleRobotLeftPso[0][k] = PF_ParticleRobot[i][0][k] + positionDifference2[0][k];
						PF_ParticleRobotLeftPso[1][k] = PF_ParticleRobot[i][1][k] + positionDifference2[1][k];
						PF_ParticleRobotLeftPso[2][k] = PF_ParticleRobot[i][2][k] + robotThtDiff;
					}

					if (PF_ParticleAntennaLeftPso.size() < (j + 1))
						PF_ParticleAntennaLeftPso.resize(j + 1);
					vector<vector<double>> temp(2, vector<double>(PF_Count));
					for (int k = 0; k < PF_Count; ++k)
					{
						//x坐标
						temp[0][k] = PF_ParticleRobotLeftPso[0][k] - antennaHLeftError * cos(PI - antennaAlphaLeftError - PF_ParticleRobotLeftPso[2][k]);
						//y坐标
						temp[1][k] = PF_ParticleRobotLeftPso[1][k] + antennaHLeftError * sin(PI - antennaAlphaLeftError - PF_ParticleRobotLeftPso[2][k]);
					}
					PF_ParticleAntennaLeftPso[j] = temp;

					//前一个时刻，粒子状态--左天线
					double robotXtAssumePre = trackMobileRobotLeft[0][referenceTagNumLeft[numberFlagLeft[0][j]][j]];
					double robotYtAssumePre = trackMobileRobotLeft[1][referenceTagNumLeft[numberFlagLeft[0][j]][j]];
					double robotThtAssumePre = trackMobileRobotLeft[2][referenceTagNumLeft[numberFlagLeft[0][j]][j]];

					robotXtDiff = robotXtAssumePre - robotXtAssume[i];
					robotYtDiff = robotYtAssumePre - robotYtAssume[i];
					robotThtDiff = robotThtAssumePre - robotThtAssume[i];

					positionDifference1[0] = robotXtDiff * cos(-robotThtAssume[i]) - robotYtDiff * sin(-robotThtAssume[i]);
					positionDifference1[1] = robotXtDiff * sin(-robotThtAssume[i]) + robotYtDiff * cos(-robotThtAssume[i]);

					for (int k = 0; k < PF_Count; ++k)
					{
						positionDifference2[0][k] = positionDifference1[0] * cos(PF_ParticleRobot[i][2][k]) - positionDifference1[1] * sin(PF_ParticleRobot[i][2][k]);
						positionDifference2[1][k] = positionDifference1[0] * sin(PF_ParticleRobot[i][2][k]) + positionDifference1[1] * cos(PF_ParticleRobot[i][2][k]);
					}

					vector<vector<double>> PF_ParticleRobotLeftPre(3, vector<double>(PF_Count));
					for (int k = 0; k < PF_Count; ++k)
					{
						PF_ParticleRobotLeftPre[0][k] = PF_ParticleRobot[i][0][k] + positionDifference2[0][k];
						PF_ParticleRobotLeftPre[1][k] = PF_ParticleRobot[i][1][k] + positionDifference2[1][k];
						PF_ParticleRobotLeftPre[2][k] = PF_ParticleRobot[i][2][k] + robotThtDiff;
					}

					if (PF_ParticleAntennaLeftPre.size() < (j + 1))
						PF_ParticleAntennaLeftPre.resize(j + 1);
					for (int k = 0; k < PF_Count; ++k)
					{
						//x坐标
						temp[0][k] = PF_ParticleRobotLeftPre[0][k] - antennaHLeftError * cos(PI - antennaAlphaLeftError - PF_ParticleRobotLeftPre[2][k]);
						//y坐标
						temp[1][k] = PF_ParticleRobotLeftPre[1][k] + antennaHLeftError * sin(PI - antennaAlphaLeftError - PF_ParticleRobotLeftPre[2][k]);
					}
					PF_ParticleAntennaLeftPre[j] = temp;

					//观测相位梯度
					PF_ObserveGradientLeft[j] = PF_ObserveNativeLeft[numberFlagLeft[1][j]][j] - PF_ObserveNativeLeft[numberFlagLeft[0][j]][j];

					//理论相位梯度
					vector<double> PF_DistanceLeftAntennaPre(PF_Count);
					vector<double> PF_DistanceLeftAntennaPso(PF_Count);
					for (int k = 0; k < PF_Count; ++k)
					{
						PF_DistanceLeftAntennaPre[k] = sqrt(pow((referenceTag[j][0] - PF_ParticleAntennaLeftPre[j][0][k]), 2) + pow((referenceTag[j][1] - PF_ParticleAntennaLeftPre[j][1][k]), 2) + pow((referenceTag[j][2] - myVisionLeftAntennaHigh), 2)) * 2;
						PF_DistanceLeftAntennaPso[k] = sqrt(pow((referenceTag[j][0] - PF_ParticleAntennaLeftPso[j][0][k]), 2) + pow((referenceTag[j][1] - PF_ParticleAntennaLeftPso[j][1][k]), 2) + pow((referenceTag[j][2] - myVisionLeftAntennaHigh), 2)) * 2;
						PF_PredictionLeft[j][2][k] = (PF_DistanceLeftAntennaPso[k] - PF_DistanceLeftAntennaPre[k]) * 2 * PI / waveLengthVar[0];
					}

					//粒子权重评估：用相位差
					for (int k = 0; k < PF_Count; ++k)
					{
						PF_DistanceLeft[j][k] = PF_PredictionLeft[j][2][k] - PF_ObserveGradientLeft[j];
						//求权重
						PF_W_Left[j][k] = (1 / sqrt(2 * PF_R) / sqrt(2 * PI)) * exp(-pow(PF_DistanceLeft[j][k], 2) / 2 / (2 * PF_R));
						PF_W_LeftImprovement[k] = PF_W_Left[j][k];
					}
				}
				else
					for (auto& x : PF_W_Left[j])
						x = -1;
			}

			//下面的是不是太过繁琐了？写成一个函数？
			//粒子权重评估
			vector<vector<double>> PF_W_PrePe;
			int n = 0;
			count = 0;
			for (auto x : PF_W_RightImprovement)
				if (x == -1)
					++count;
			if (count != PF_Count)
			{
				++n;
				PF_W_PrePe.resize(n);
				PF_W_PrePe[n - 1] = PF_W_RightImprovement;
			}
			count = 0;
			for (auto x : PF_W_LeftImprovement)
				if (x == -1)
					++count;
			if (count != PF_Count)
			{
				++n;
				PF_W_PrePe.resize(n);
				PF_W_PrePe[n - 1] = PF_W_LeftImprovement;
			}

			vector<double> PF_W_Pre(PF_Count, 1);
			if (PF_W_PrePe.size() == 0)
			{
				phaseFlag[i] = false;
			}
			else if (PF_W_PrePe.size() == 1)
			{
				double sum = 1;
				for (int j = 0; j < PF_Count; ++j)
					sum *= PF_W_PrePe[0][j];
				for (int j = 0; j < PF_Count; ++j)
					PF_W_Pre[j] = sum;
				phaseFlag[i] = true;
			}
			else
			{
				for (int j = 0; j < PF_Count; ++j)
					PF_W_Pre[j] = PF_W_PrePe[0][j] * PF_W_PrePe[1][j];
				phaseFlag[i] = true;
			}

			for (int j = 0; j < PF_Count; ++j)
			{
				PF_W[i][0][j] = PF_W_Pre[j] + 1e-99;
				PF_W[i][1][j] = PF_W[i][0][j] * PF_W[i - 1][1][j] * PF_W_RSSI[i][j];
			}

			double sum = accumulate(PF_W[i][1].begin(), PF_W[i][1].end(), 0.0);
			for (auto& x : PF_W[i][1])
				x /= sum;

			//x，y坐标的均值及其协方差
			//所有粒子的中心位置--x，y坐标
			for (int j = 0; j < 2; ++j)
				for (int k = 0; k < PF_Count; ++k)
					PF_CenterMean[j][i] += PF_ParticleRobot[i][j][k] * PF_W[i][1][k];
			double weightSqr = 0;
			for (int j = 0; j < PF_Count; ++j)
				weightSqr += PF_W[i][1][j] * PF_W[i][1][j];

			double factor = 0;
			if (abs(weightSqr - 1.0) < sqrt(DBL_EPSILON))
				factor = 1.0;
			else
				factor = 1 / (1 - weightSqr);

			vector<vector<double>> meanDiff(2, vector<double>(PF_Count));
			for (int j = 0; j < 2; ++j)
				for (int k = 0; k < PF_Count; ++k)
					meanDiff[j][k] = PF_ParticleRobot[i][j][k] - PF_CenterMean[j][i];

			for (int j = 0; j < 2; ++j)
				for (int k = 0; k < 2; ++k)
				{
					for (int m = 0; m < PF_Count; ++m)
						PF_CenterVar[i][j][k] += meanDiff[j][m] * PF_W[i][1][m] * meanDiff[k][m];
					PF_CenterVar[i][j][k] *= factor;
				}

			//角度均值及协方差：
			double sinsum = 0;
			for (int j = 0; j < PF_Count; ++j)
				sinsum += PF_W[i][1][j] * sin(PF_ParticleRobot[i][2][j]);
			double cossum = 0;
			for (int j = 0; j < PF_Count; ++j)
				cossum += PF_W[i][1][j] * cos(PF_ParticleRobot[i][2][j]);
			double resultantLength = sqrt(pow(sinsum, 2) + pow(cossum, 2));
			PF_CenterMean[2][i] = atan2(sinsum, cossum);
			PF_CenterVar[i][2][2] = -2 * log(resultantLength);
			for (int j = 0; j < PF_Count; ++j)
				PF_CenterMean[3][i] += PF_ParticleRobot[i][2][j] * PF_W[i][1][j];

			//判断当前是否启动重采样
			double squaredSum = weightSqr;

			//有效粒子的数量
			double neff = 1 / squaredSum;
			neffRatioPre[i] = neff / PF_Count;
			//cout << neffRatioPre[i] << ' ';
			//记得还原！！！
			if (neffRatioPre[i] < neffRatioThreshold)
			{
				//执行系统重采样
				int numNewParticles = PF_Count;
				//double randNum = rand() / double(RAND_MAX);
				vector<double> randSamples(numNewParticles);
				for (int j = 0; j < numNewParticles; ++j)
					randSamples[j] = j * (1 - 1.0 / numNewParticles) / (numNewParticles - 1) + randNum[i];//测试
					//randSamples[j] = j * (1 - 1.0 / numNewParticles) / (numNewParticles - 1) + randNum / numNewParticles;

				int n = 0;
				int index = 0;
				// Q = cumsum(PF_w(2,:,i));
				vector<double> Q(PF_Count);
				for (int j = 0; j < PF_Count; ++j)
					Q[j] = accumulate(PF_W[i][1].begin(), PF_W[i][1].begin() + j + 1, 0.0);
				vector<vector<double>> PF_ParticleNew(3, vector<double>(numNewParticles));
				while (n < numNewParticles)
				{
					while (Q[index] <= randSamples[n])
						index = index % (PF_Count - 1) + 1;
					//得到新粒子集
					for (int j = 0; j < 3; ++j)
						PF_ParticleNew[j][n] = PF_ParticleRobot[i][j][index];
					++n;
				}
				//得到最后的重采样粒子集
				PF_ParticleRobot[i] = PF_ParticleNew;
				//求本次的最终权重
				for (auto& x : PF_W[i][1])
					x = 1.0 / PF_Count;
				PF_W_Slow = 0;
				PF_W_Fast = 0;
				PF_ReSample[i] = 1;
			}
			else
			{
				//不执行重采样
				PF_ReSample[i] = -1;
			}

			//定位误差
			robotXt[i] = mobileRobotPoseVisionRight[numberFlag][0];
			robotYt[i] = mobileRobotPoseVisionRight[numberFlag][1];
			robotTht[i] = mobileRobotPoseVisionRight[numberFlag][2];

			antennaLeft[0][i] = robotXt[i] - antennaHLeftError * cos(PI - antennaAlphaLeftError - robotTht[i]);
			antennaLeft[1][i] = robotYt[i] + antennaHLeftError * sin(PI - antennaAlphaLeftError - robotTht[i]);
			antennaRight[0][i] = robotXt[i] + antennaHRightError * cos(-antennaAlphaRightError + robotTht[i]);
			antennaRight[1][i] = robotYt[i] + antennaHRightError * sin(-antennaAlphaRightError + robotTht[i]);

			robotEstimationError[0][i] = PF_CenterMean[0][i] - robotXt[i];
			robotEstimationError[1][i] = PF_CenterMean[1][i] - robotYt[i];
			robotEstimationError[2][i] = sqrt(pow(robotEstimationError[0][i], 2) + pow(robotEstimationError[1][i], 2));

			robotEstimationError[3][i] = PF_CenterMean[2][i] - robotTht[i];

			robotEstimationError[4][i] = PF_CenterMean[3][i] - robotTht[i];

		}
		//绘图
		// 	figure(1)
		//     hold off
		//     plot(reference_tag(:,1), reference_tag(:,2), 'r.', 'markersize',50);  %系统状态位置
		//     hold on

		//     color = PF_w(2, :,i);
		//     scatter(PF_particle_robot(1, :, i), PF_particle_robot(2, :, i), 'filled','cdata',color);   %各个粒子位置
		//     plot(PF_center_mean(1, 1:i),PF_center_mean(2, 1:i),'b.', 'markersize',25);%粒子中心位置
		//     plot(robot_xt(1:i),robot_yt(1:i),'r.');%机器人的位置
		//     plot(antenna_left(1, 1:i),antenna_left(2, 1:i),'k.');%左天线位置
		//     plot(antenna_right(1, 1:i),antenna_right(2, 1:i),'k.');%右天线位置

		//     quiver(PF_center_mean(1, i),PF_center_mean(2,i), cos(PF_center_mean(3, i)), sin(PF_center_mean(3, i)),10, 'color', 'c', 'linewidth', 3);
		//     quiver(PF_particle_robot(1, :, i), PF_particle_robot(2, :, i), cos(PF_particle_robot(3, :, i)), sin(PF_particle_robot(3, :, i)),1, 'color', 'k', 'linewidth', 3);
		//     set(gca,'DataAspectRatio',[1 1 1])

		//     axis([PF_scope(1,1)-200 400 PF_scope(2,1)-100 PF_scope(2,2)+100]);

		//更新频率把控
		auto endTime = clock();
		//time[i] = (endTime - startTime) / CLOCKS_PER_SEC;
		//在Windows下时间按照毫秒计算，在Linux下时间按照秒计算(???)
		//Sleep((1 - time[i]) * 1000);
		// Sleep(1 - time[i]);
	}

	//销毁
	for (int i = 0; i < sizeof(myVisionTotal) / sizeof(double**); i++)
	{
		mxArray* pF = mxGetField(pS, 0, name[i]);
		auto row = mxGetM(pF);
		auto col = mxGetN(pF);
		for (int j = 0; j < row; j++)
			delete[] * (*(q + i) + j);
		delete[] * (q + i);
	}
	matClose(pM);

	/*数据验证：将外部数据导入并与之作差，将差值与阈值对比，若大于阈值则计数*/
	/*外部数据存储*/
	vector<vector<vector<double>>> PF_ParticleRobot0(500, vector<vector<double>>(3, vector<double>(PF_Count)));
	vector<vector<vector<double>>> PF_W0(Times, vector<vector<double>>(2, vector<double>(PF_Count)));
	vector<vector<double>> PF_CenterMean0(4, vector<double>(Times));
	vector<double> robotXt0(Times);
	vector<double> robotYt0(Times);
	vector<vector<double>> antennaLeft0(2, vector<double>(Times));
	vector<vector<double>> antennaRight0(2, vector<double>(Times));

	mxArray* pV1 = matGetVariable(pM0, "PF_particle_robot");
	double* pT1 = (double*)mxGetData(pV1);

	mxArray* pV2 = matGetVariable(pM0, "PF_w");
	double* pT2 = (double*)mxGetData(pV2);

	mxArray* pV3 = matGetVariable(pM0, "PF_center_mean");
	double* pT3 = (double*)mxGetData(pV3);

	mxArray* pV4 = matGetVariable(pM0, "robot_xt");
	double* pT4 = (double*)mxGetData(pV4);

	mxArray* pV5 = matGetVariable(pM0, "robot_yt");
	double* pT5 = (double*)mxGetData(pV5);

	mxArray* pV6 = matGetVariable(pM0, "antenna_left");
	double* pT6 = (double*)mxGetData(pV6);

	mxArray* pV7 = matGetVariable(pM0, "antenna_right");
	double* pT7 = (double*)mxGetData(pV7);

	for (int i = 0; i < Times; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			for (int k = 0; k < 500; ++k)
			{
				PF_ParticleRobot0[i][j][k] = pT1[i * 500 * 3 + k * 3 + j];
			}
		}
		for (int j = 0; j < 2; ++j)
		{
			for (int k = 0; k < 500; ++k)
			{
				PF_W0[i][j][k] = pT2[i * 500 * 2 + k * 2 + j];
			}
		}
		robotXt0[i] = pT4[i];
		robotYt0[i] = pT5[i];
	}
	for (int i = 0; i < 4; ++i)
	{
		for (int j = 0; j < Times; ++j)
		{
			PF_CenterMean0[i][j] = pT3[j * 4 + i];
		}
	}
	for (int i = 0; i < 2; ++i)
	{
		for (int j = 0; j < Times; ++j)
		{
			antennaLeft0[i][j] = pT6[j * 2 + i];
			antennaRight0[i][j] = pT7[j * 2 + i];
		}
	}

	/*与外部数据作差并计数*/
	int count = 0;
	for (int i = 0; i < Times; ++i)
	{
		//cout << i << endl;
		for (int j = 0; j < 3; ++j)
		{
			for (int k = 0; k < 500; ++k)
			{
				if (abs(PF_ParticleRobot0[i][j][k] - PF_ParticleRobot[i][j][k]) > ep)
				{
					++count;
					//cout << abs(PF_ParticleRobot0[i][j][k] - PF_ParticleRobot[i][j][k]) << ' ';
				}
					
			}
			//cout << endl;
		}
	}
	cout << count << endl;
	
	for (int i = 0; i < Times; ++i)
	{
		//cout << i << endl;
		for (int j = 0; j < 2; ++j)
		{
			for (int k = 0; k < 500; ++k)
			{
				if (abs(PF_W0[i][j][k] - PF_W[i][j][k]) > ep)
				{
					++count;
					//cout << abs(PF_W0[i][j][k] - PF_W[i][j][k]) << ' ';
				}
			}
			//cout << endl;
		}
		if (abs(robotXt0[i] - robotXt[i]) > ep)
		{
			++count;
			//cout << abs(robotXt0[i] - robotXt[i]) << ' ';
		}
			
		if (abs(robotYt0[i] - robotYt[i]) > ep)
			++count;
	}
	cout << count << endl;

	for (int i = 0; i < 4; ++i)
	{
		for (int j = 0; j < Times; ++j)
		{
			if (abs(PF_CenterMean0[i][j] - PF_CenterMean[i][j]) > ep)
				++count;
		}
	}
	cout << count << endl;

	for (int i = 0; i < 2; ++i)
	{
		for (int j = 0; j < Times; ++j)
		{
			if (abs(antennaLeft0[i][j] - antennaLeft[i][j]) > ep)
				++count;
			if (abs(antennaRight0[i][j] - antennaRight[i][j]) > ep)
				++count;
		}
	}
	cout << count << endl;

	return 0;
}
int mymax(const vector<int>& v)
{
	int m = v.front();
	for (auto x : v)
		m = max(m, x);
	return m;
}
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
double distance(const vector<double>& v1, const vector<double>& v2)
{
	double sum = 0;
	for (int i = 0; i < v1.size(); ++i)
		sum += (v1[i] - v2[i]) * (v1[i] - v2[i]);
	return sqrt(sum);
}
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
vector<bool> inpolygon(const vector<double>& xp, const vector<double>& yp, const vector<double>& xl, const vector<double>& yl)
{
	vector<bool> result(xp.size());
	for (int i = 0; i < xp.size(); ++i)
	{
		bool flag = false;	   //判断结果（true；点落在多边形内；false:点未落在多边形内）
		int k = xl.size() - 1; //是多边形的最后一个顶点
		for (int j = 0; j < xl.size(); ++j)
		{
			//判断点是否在线段的两侧
			if ((yl[j] < yp[i] && yl[k] >= yp[i]) || (yl[k] < yp[i] && yl[j] >= yp[i]))
			{
				//根据两点式方程计算出过点P且平行于X轴的直线与线段的交点，两点式方程：x = x1 +  (y - y1) * (x2 - x1) / (y2 - y1);
				if (xl[j] + (yp[i] - yl[j]) * (xl[k] - xl[j]) / (yl[k] - yl[j]) < xp[i])
					flag = !flag;
			}
			//进行下一线段判断
			k = j;
		}
		result[i] = flag;
	}
	return result;
}
