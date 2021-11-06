#pragma once
#include <vector>
#include <iostream>
#include <fstream>
using namespace std;

//PI
const double PI = 3.1415926535898;
//这个没有用到
const double Resolution = 0.0015;
//光速
const int vLight = 300000000;
//单位 raidan/s,标准差，应该是观测噪声的标准差
const double PhaseGauss = 0.1;

struct Data
{
	vector<int> no;				//1
	vector<double> phase;		//2
	vector<double> myThree;		//3
	vector<double> myFour;		//4
	vector<double> Xcoordinate; //5
	vector<double> Ycoordinate; //6
	vector<double> Zcoordinate; //7
	vector<double> readerTime;	//8
	vector<double> windowsTime; //9
	vector<double> myTen;		//10
	vector<int> RSSI;			//11
};
struct TextData
{
	vector<int> no;
	vector<string> epcData;
};
//struct MyDataResult
//{
//	Data data;
//	TextData textData;
//};
//class Data
//{
//	vector<int> no;				//1
//	vector<double> phase;		//2
//	vector<double> myThree;		//3
//	vector<double> myFour;		//4
//	vector<double> Xcoordinate; //5
//	vector<double> Ycoordinate; //6
//	vector<double> Zcoordinate; //7
//	vector<double> readerTime;	//8
//	vector<double> windowsTime; //9
//	vector<double> myTen;		//10
//	vector<int> RSSI;			//11
//};
//class TextData
//{
//	vector<int> no;
//	vector<string> epcData;
//};
class MyDataResult
{
public :
	MyDataResult(string textPath)
	{
		fstream fin;
		fin.open(textPath, ios::in);
		while (1)
		{
			int n;
			string s;
			double d;
			//到达文件结尾需再次读入，eof才会变为true
			fin >> n;
			if (fin.eof())
				break;
			this->textData.no.push_back(n);
			fin >> s;
			this->textData.epcData.push_back(s);
			fin >> n;
			this->data.no.push_back(n);
			fin >> d;
			this->data.phase.push_back(d);
			fin >> d;
			this->data.myThree.push_back(d);
			fin >> d;
			this->data.myFour.push_back(d);
			fin >> d;
			this->data.Xcoordinate.push_back(d);
			fin >> d;
			this->data.Ycoordinate.push_back(d);
			fin >> d;
			this->data.Zcoordinate.push_back(d);
			fin >> d;
			this->data.readerTime.push_back(d);
			fin >> d;
			this->data.windowsTime.push_back(d);
			fin >> d;
			this->data.myTen.push_back(d);
			fin >> n;
			this->data.RSSI.push_back(n);
		}
		fin.close();
	}
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