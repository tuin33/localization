#pragma once
#include <vector>
using std::vector;

//获取向量中最大值 的一个函数
int mymax(const vector<int>&);
//正态分布生成函数，生成期望为mu，标准差为sigma的正态分布
double gaussrand(const double& mu, const double& sigma);
//计算距离，v1与v2是两个大小相同的向量
double distance(const vector<double>& v1, const vector<double>& v2);
//实现k-临近算法
vector<int> knnsearch(vector<vector<double>> testData, vector<double> targetData, int k);
//判断点是否落入多边形内
vector<bool> inpolygon(const vector<vector<double>>& pointData, const vector<int>& optionalTagFlag, const vector<vector<double>>& testData,const int& numOfTag);
//找到里程计的时间在视觉时间序列中的前后两个最近的位姿点
int findpso(const vector<double>& visionTimeDiff, const vector<double>& odometerTimeDiff, const int& index, bool& flag);
int findpre(const vector<double>& visionTimeDiff, const vector<double>& odometerTimeDiff, const int& index, bool& flag);
//计算机器人所在大致范围（用于初始化粒子）
vector<vector<double>> calPF_Scope(vector<vector<double>>& referenceTag, vector<int>& randData, vector<int>& colLeft, vector<int>& colRight);
//读取.mat文件
//void readmat(const char *matPath,const char **fieldName,double ***q,int len);