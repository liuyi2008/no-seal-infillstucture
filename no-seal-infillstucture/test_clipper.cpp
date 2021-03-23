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
vector<point>asd = { {10.0,10.0},{10.0,10.0} };//��ʼ������

vector<vector<point> >model = { {{10.0,10.0},{10.0,60.0},{60.0,60.0},{60.0,10.0} }, {{ 20.0,30.0 },{ 20.0,40.0 },{ 30.0,40.0 },{ 30.0,30.0 },{-10000,-10000},{ 40.0,30.0 },{ 40.0,40.0 },{ 50.0,40.0 },{ 50.0,30.0 }} }; //slice  ����ģ�͵ģ��ֲ�����������꣬һά��������ά�˲�������ߵ�����

vector<Paths>model2;//����Ƭ��model��ת�����͵õ��ģ���ƫ��֮����������
vector<Paths>model3;//��С����һȦ��������ȡ����

//����ע�⣬model��2ά��i���j���㣩����model2��3ά��i���m��path�ϵĵ�n���㣩
void main()
{
	ClipperOffset co;
	Paths solution;
	for (int i = 0; i < model.size(); i++)//����model����i��
	{
		int m = 0, n = 0;

		for (int j = 0; j < model[i].size(); j++)//����model����i���еĵ�j����
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

		for (int j = 0; j < model2[i].size(); j++)//��ʱ��������model2��ĳ����ĵ�j��path��ƫ��
		{
			co.AddPath(model2[i][j], jtRound, etClosedPolygon);
			co.Execute(solution, -7.0);
			for (int k = 0; k < solution.size(); k++)//����solution�����ӵ���j��path�ĺ���
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