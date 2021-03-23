#include "clipper.hpp"  
#include"clipper.cpp"
#include <gl/glut.h> 
#include<stdio.h>
#include <vector>
#define use_int32
struct point {double x, y;}p, p1;
point p = { (10.0, 0.0) };
using namespace ClipperLib;
using namespace std;

vector<int> ilist = { 1,2,4,5,6,7 };
vector<point>asd = { {10.0,10.0},{10.0,10.0} };//初始化例子

vector<vector<point> >model = { {{10.0,10.0},{10.0,60.0},{60.0,60.0},{60.0,10.0} }, {{ 20.0,30.0 },{ 20.0,40.0 },{ 30.0,40.0 },{ 30.0,30.0 },{-10000,-10000},{ 40.0,30.0 },{ 40.0,40.0 },{ 50.0,40.0 },{ 50.0,30.0 }} }; //slice  整个模型的，分层的轮廓线坐标，一维层数，二维此层的轮廓线点坐标

vector<Paths>model2;//从切片（model）转换类型得到的，做偏移之后画外轮廓。
vector<Paths>model3;//最小的那一圈，用来做取交。

//这里注意，model是2维（i层第j个点），而model2是3维（i层第m条path上的第n个点）
void main()
{
	ClipperOffset co;
	Paths solution;
	for (int i = 0; i < model.size(); i++)//遍历model，第i层
	{
		int m = 0, n = 0;

		for (int j = 0; j < model[i].size(); j++)//遍历model，第i层中的第j个点
		{
			if(model[i][j].x != -10000)
			{
				model2[i][m][n].X = cInt(model[i][j].x * 100);
				model2[i][m][n].Y = cInt(model[i][j].y * 100);
				n++;
			}
			else
			{
				m++;
				n = 0;
			}
		}

		for (int j = 0; j < model2[i].size(); j++)//此时遍历的是model2，某层里的第j条path做偏移
		{
			co.AddPath(model2[i][j], jtRound, etClosedPolygon);
			co.Execute(solution, -7.0);
			for (int k = 0; k < solution.size(); k++)//遍历solution，都加到第j条path的后面
			{
				model2[i][j].insert(model2[i][j].end(), solution[k].begin(), solution[k].end());
			}
			co.Clear();
			solution.clear();
		}
	}


	//draw solution ...
	printf("1");
	//DrawPolygon(solution);
}