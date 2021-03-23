#include <iostream>  
#include <OpenMesh/Core/IO/MeshIO.hh>  
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>  
#include <OpenMesh/Core/Mesh/Handles.hh>
#include <OpenMesh\Core\Utils\Property.hh> 
#include <math.h>
#include <string>
#include <vector>
#include <algorithm>
#include<stdio.h>
#include <gl/glut.h>  
//#include<GL/glut.h>
//#include "GL\glut.h"
//#include "GL\freeglut.h"
using namespace std;

struct point {

	double x, y;

}p, p1;

typedef unsigned char boolean;

#define eps 1e-8
#define zero(x) (((x)>0?(x):-(x))<eps)

vector<point>coord;  //����һ���������
vector<point>polypoint; //һ�������������
vector<point>tripoint;  //һ��ı仯�� �߶����� �ɶ�
vector<point>new_tripoint;  //�µ�һ��ı仯�� �߶����� �ɶ�
vector<point>interpoint; //��������
vector<vector<point> >model; //slice  ����ģ�͵ģ��ֲ�����������꣬һά��������ά�˲�������ߵ�����
vector<vector<point> >modelfill(model); //��ɵ�model����ߣ�һά��������ά�˲�ı仯�� �߶�����


const string file_1 = "cube.obj";
const string file_2 = "bunnyr.obj";
const string file_3 = "check.obj";


//������ӽ��й�
double xs = 0;
double ys = 0;
double zs = 0;
double scale = 0;

//������ӵ��й�
double sx = 0;
double sy = 0;
double sz = 0;

typedef OpenMesh::TriMesh_ArrayKernelT<> MyMesh;//��������MyMesh�࣬�����������Զ�������
MyMesh mesh;

//���ڰ�Χ�������
double Bmax_x, Bmax_y, Bmax_z, Bmin_x, Bmin_y, Bmin_z, px, py, pz;

int countt = 0;


// ��ȡ�ļ��ĺ���
void readfile(string file) {//�ļ���ȡ����mesh����
	// ���󶥵㷨�� vertex normals
	mesh.request_vertex_normals();
	//��������ڶ��㷨�ߣ��򱨴� 
	if (!mesh.has_vertex_normals())
	{
		cout << "���󣺱�׼�������� �����ߡ�������" << endl;
		return;
	}
	// ����ж��㷢�����ȡ�ļ� 
	OpenMesh::IO::Options opt;
	if (!OpenMesh::IO::read_mesh(mesh, file, opt))
	{
		cout << "�޷���ȡ�ļ�:" << file << endl;
		return;
	}
	else cout << "�ɹ���ȡ�ļ�:" << file << endl;
	cout << endl; // Ϊ��ui��ʾ�ÿ�һЩ
				  //��������ڶ��㷨�ߣ�������
	if (!opt.check(OpenMesh::IO::Options::VertexNormal))
	{
		// ͨ���淨�߼��㶥�㷨��
		mesh.request_face_normals();
		// mesh��������㷨��
		mesh.update_normals();
		// �ͷ��淨��
		mesh.release_face_normals();
	}
}

//ͨ���򵥵ı�������ģ�Ͱ�Χ�У����ȡ�����޹أ�
void BoundingBox() {
	MyMesh::Point pt;
	int st = 0;
	for (auto it = mesh.vertices_begin(); it != mesh.vertices_end(); ++it) {
		MyMesh::VertexHandle vh_i = *it;
		pt = mesh.point(vh_i);
		px = pt.data()[0];
		py = pt.data()[1];
		pz = pt.data()[2];
		if (st == 0) {
			Bmax_x = Bmin_x = px;
			Bmax_y = Bmin_y = py;
			Bmax_z = Bmin_z = pz;
			st++;
		}
		else {
			if (px > Bmax_x)Bmax_x = px; else if (px < Bmin_x)Bmin_x = px;
			if (py > Bmax_y)Bmax_y = py; else if (py < Bmin_y)Bmin_y = py;
			if (pz > Bmax_z)Bmax_z = pz; else if (pz < Bmin_z)Bmin_z = pz;
		}
	}

}

//���溯��
bool IntersectPlane(MyMesh::Point pt, MyMesh::Point pnorm) //pt������һ�㣬pnorm���淨��
{
	const float ERR = 0.001;
	//�������� pt��pnorm��*pilist��pnum[]  ���庯��ԭ�� �� �����㷨.docx
	int starte, ne, ne1, nf;
	MyMesh::Point vt1, vt2;
	//MyMesh::Face f1;
	MyMesh::HalfedgeHandle nhe;
	MyMesh::FaceHandle nhf;
	float d1, d2, sd1, sd2;
	bool *flag, suc;
	float dist, mind = 1.0e+8;

	sd1 = sd2 = -10000;
	int esize = mesh.n_halfedges();
	flag = new bool[esize];

	suc = false;


	for (MyMesh::HalfedgeIter it = mesh.halfedges_begin(); it != mesh.halfedges_end(); ++it) //��������ģ�����еıߣ��н���İ�id��¼��flag��
	{

		MyMesh::HalfedgeHandle hh = *it;
		int id = hh.idx();
		flag[id] = false;

		auto fromVertex = mesh.from_vertex_handle(hh);
		auto toVertex = mesh.to_vertex_handle(hh);
		vt1 = mesh.point(fromVertex);
		vt2 = mesh.point(toVertex);
		//printf("$ %.3f %.3f $\n", vt1.data()[0],vt2.data()[0]);
		d1 = pnorm.data()[0] * (vt1.data()[0] - pt.data()[0]) + pnorm.data()[1] * (vt1.data()[1] - pt.data()[1])
			+ pnorm.data()[2] * (vt1.data()[2] - pt.data()[2]);//d1����ľ���
		d2 = pnorm.data()[0] * (vt2.data()[0] - pt.data()[0]) + pnorm.data()[1] * (vt2.data()[1] - pt.data()[1])
			+ pnorm.data()[2] * (vt2.data()[2] - pt.data()[2]);

		if (((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0) || d1 > 0 && d2 == 0 || d1 == 0 && d2 > 0))//�߶������ཻ
		{
			flag[id] = true;

			vt1.data()[0] = vt1.data()[0] - pt.data()[0];
			vt1.data()[1] = vt1.data()[1] - pt.data()[1];
			vt1.data()[2] = vt1.data()[2] - pt.data()[2];       // point date minus point date 
			dist = vt1.data()[0] * vt1.data()[0] + vt1.data()[1] * vt1.data()[1] + vt1.data()[2] * vt1.data()[2];

			if (dist < mind)
			{
				nhe = hh;	//��̱�
				mind = dist;//��̾���
				ne = id;    //������ڱߵı��               //  printf("ne:  %d  \n", ne);
				suc = true;
			}
		}
	}

	if (!suc)
	{
		delete[]flag;
		return false; //û�н��㣬����return false��������������
	}

	starte = ne;//���ѭ����ʼ�ı�

	suc = false;

	nhf = mesh.face_handle(nhe);//��̱�������

	while (!suc)
	{
		//printf("%%%%");	

		auto fromVertex = mesh.from_vertex_handle(nhe);
		auto toVertex = mesh.to_vertex_handle(nhe);

		vt1 = mesh.point(fromVertex);
		vt2 = mesh.point(toVertex);

		d1 = pnorm.data()[0] * (vt1.data()[0] - pt.data()[0]) + pnorm.data()[1] * (vt1.data()[1] - pt.data()[1])
			+ pnorm.data()[2] * (vt1.data()[2] - pt.data()[2]);
		d2 = pnorm.data()[0] * (vt2.data()[0] - pt.data()[0]) + pnorm.data()[1] * (vt2.data()[1] - pt.data()[1])
			+ pnorm.data()[2] * (vt2.data()[2] - pt.data()[2]);
		//printf("$$$%lf %lf \n", d1, d2);
		if ((sd1 == d1) && (sd2 == d2))
		{
			flag[ne] = false;
		}
		sd1 = d1; sd2 = d2;
		//pilist[num].data()[0] = (float)fabs(d1) / ((float)fabs(d1) + (float)fabs(d2))*(vt2.data()[0] - vt1.data()[0]) + vt1.data()[0];
		//pilist[num].data()[1] = (float)fabs(d1) / ((float)fabs(d1) + (float)fabs(d2))*(vt2.data()[1] - vt1.data()[1]) + vt1.data()[1];
		//pilist[num].data()[2] = (float)fabs(d1) / ((float)fabs(d1) + (float)fabs(d2))*(vt2.data()[2] - vt1.data()[2]) + vt1.data()[2];

		p = { (float)fabs(d1) / ((float)fabs(d1) + (float)fabs(d2))*(vt2.data()[0] - vt1.data()[0]) + vt1.data()[0] ,(float)fabs(d1) / ((float)fabs(d1) + (float)fabs(d2))*(vt2.data()[1] - vt1.data()[1]) + vt1.data()[1] };
		coord.push_back(p);


		do {
			for (auto it = mesh.fh_begin(nhf); it != mesh.fh_end(nhf); ++it) //nhf��̱�������,���������ıߣ�ֻ��3��
			{
				MyMesh::HalfedgeHandle halfnow = *it;

				const int ne1 = halfnow.idx();

				if (flag[ne1] == false || ne == ne1) continue;

				MyMesh::VertexHandle fromV = mesh.from_vertex_handle(halfnow);
				MyMesh::VertexHandle toV = mesh.to_vertex_handle(halfnow);

				vt1 = mesh.point(fromV);
				vt2 = mesh.point(toV);

				d1 = pnorm.data()[0] * (vt1.data()[0] - pt.data()[0]) + pnorm.data()[1] * (vt1.data()[1] - pt.data()[1])
					+ pnorm.data()[2] * (vt1.data()[2] - pt.data()[2]);
				d2 = pnorm.data()[0] * (vt2.data()[0] - pt.data()[0]) + pnorm.data()[1] * (vt2.data()[1] - pt.data()[1])
					+ pnorm.data()[2] * (vt2.data()[2] - pt.data()[2]);

				p = { (float)fabs(d1) / ((float)fabs(d1) + (float)fabs(d2))*(vt2.data()[0] - vt1.data()[0]) + vt1.data()[0] ,(float)fabs(d1) / ((float)fabs(d1) + (float)fabs(d2))*(vt2.data()[1] - vt1.data()[1]) + vt1.data()[1] };
				coord.push_back(p);

				MyMesh::HalfedgeHandle halfnext = mesh.opposite_halfedge_handle(halfnow);//��ȡ������

				nhf = mesh.face_handle(halfnext);//������������ڵ���

				int ne2 = halfnext.idx();

				flag[ne1] = flag[ne2] = false;//ne1,ne2�Ƕ����������ߣ��������һ���Ľ���Ͷ���Ϊfalse

				if (nhf.idx() == -1)
				{
					starte = ne;
					flag[ne] = false;
					break;
				}
				ne = ne2;//�ԶԱߵ��Ǹ���������һ�������ཻ����

				//pilist[num].data()[0] = (float)fabs(d1) / ((float)fabs(d1) + (float)fabs(d2))*(vt2.data()[0] - vt1.data()[0]) + vt1.data()[0];
				//pilist[num].data()[1] = (float)fabs(d1) / ((float)fabs(d1) + (float)fabs(d2))*(vt2.data()[1] - vt1.data()[1]) + vt1.data()[1];
				//pilist[num].data()[2] = (float)fabs(d1) / ((float)fabs(d1) + (float)fabs(d2))*(vt2.data()[2] - vt1.data()[2]) + vt1.data()[2];
				//printf("##%lf %lf %lf\n", pilist[num].data()[0], pilist[num].data()[1], pilist[num].data()[2]);

				break;
			}
		} while (ne != starte);

		suc = true;

		for (auto it = mesh.halfedges_begin(); it != mesh.halfedges_end(); ++it) //������û�еڶ�����
		{

			MyMesh::HalfedgeHandle hh = *it;

			int id = hh.idx();

			if (flag[id] == true)
			{
				ne = id;
				starte = ne;
				nhe = hh;
				nhf = mesh.face_handle(nhe);
				if (nhf.idx() == -1)
				{
					flag[ne] = false;
					continue;
				}
				//pilist[num].data()[0] = -10000;
				//pilist[num].data()[1] = -10000;
				//pilist[num].data()[2] = -10000;
				p = { -10000,-10000 };//�������м�ļ������
				coord.push_back(p);

				suc = false;
				break;
			}
		}


	};

	delete[]flag;
	return true;
}

void findIntersect(float z)
{
	for (double i = Bmin_z; i < Bmax_z; i += z)
	{
		MyMesh::Normal vf(0, 0, 1);//���淨��
		MyMesh::Point pt; //������һ��
		pt.data()[0] = 0; pt.data()[1] = 0; pt.data()[2] = i;
		IntersectPlane(pt, vf);//����һ�������ߣ�����coord

		model.push_back(coord);
		coord.clear();
	}
	countt = model.size();
}

//***********************************************************************������Ƭһ�����****************************************************************
//*************�ж��߶����߶εĹ�ϵ����*************************
double xmult(point p1, point p2, point p0) {
	return (p1.x - p0.x)*(p2.y - p0.y) - (p2.x - p0.x)*(p1.y - p0.y);
}

//�е��Ƿ����߶���,�����˵�
int dot_online_in(point p, point l1, point l2) {
	return zero(xmult(p, l1, l2)) && (l1.x - p.x)*(l2.x - p.x) < eps && (l1.y - p.y)*(l2.y - p.y) < eps;
}

//���������߶�ͬ��,�����߶��Ϸ���0
int same_side(point p1, point p2, point l1, point l2) {
	return xmult(l1, p1, l2)*xmult(l1, p2, l2) > eps;
}

//����ֱ��ƽ��
int parallel(point u1, point u2, point v1, point v2) {
	return zero((u1.x - u2.x)*(v1.y - v2.y) - (v1.x - v2.x)*(u1.y - u2.y));
}

//�����㹲��
int dots_inline(point p1, point p2, point p3) {
	return zero(xmult(p1, p2, p3));
}

//�����߶��ཻ,�����˵�Ͳ����غ�
int intersect_in(point u1, point u2, point v1, point v2) {
	if (!dots_inline(u1, u2, v1) || !dots_inline(u1, u2, v2))
		return !same_side(u1, u2, v1, v2) && !same_side(v1, v2, u1, u2);
	return dot_online_in(u1, v1, v2) || dot_online_in(u2, v1, v2) || dot_online_in(v1, u1, u2) || dot_online_in(v2, u1, u2);
}

//�������߶ν���,�����߶��Ƿ��ཻ(ͬʱ����Ҫ�ж��Ƿ�ƽ��!)
point intersection(point u1, point u2, point v1, point v2) {
	point ret = u1;
	double t = ((u1.x - v1.x)*(v1.y - v2.y) - (u1.y - v1.y)*(v1.x - v2.x))
		/ ((u1.x - u2.x)*(v1.y - v2.y) - (u1.y - u2.y)*(v1.x - v2.x));
	ret.x += (u2.x - u1.x)*t;
	ret.y += (u2.y - u1.y)*t;
	return ret;
}

//���������н���潻�㣬û�����޲���
void inter(point u1, point u2, point v1, point v2)
{
	point ans;
	if (parallel(u1, u2, v1, v2) || !intersect_in(u1, u2, v1, v2))
	{

	}
	else {
		ans = intersection(u1, u2, v1, v2);
		interpoint.push_back(ans);
		//printf("����Ϊ:(%lf,%lf)", ans.x, ans.y);
	}

}

//*************�ж��߶����߶εĹ�ϵ����*************************

int InOrOutPolygon(point a)  //�жϵ��Ƿ��ڶ������
{
	double x0 = a.x;
	double y0 = a.y;
	int crossings = 0;
	int n = polypoint.size();
	for (int i = 0; i < n; i++)
	{
		// ��������x֮�� ���Ե㴹ֱy������������
		double slope = (polypoint[(i + 1 + n) % n].y - polypoint[i].y) / (polypoint[(i + 1 + n) % n].x - polypoint[i].x);
		boolean cond1 = (polypoint[i].x <= x0) && (x0 < polypoint[(i + 1 + n) % n].x);
		boolean cond2 = (polypoint[(i + 1 + n) % n].x <= x0) && (x0 < polypoint[i].x);
		boolean above = (y0 < slope * (x0 - polypoint[i].x) + polypoint[i].y);
		if ((cond1 || cond2) && above) crossings++;
	}
	return (crossings % 2 != 0);    //�� �� ֵ:  0:�� 1:��
}

void getcoord(point A, point B) //��A��B����֮��ĵ㶼����tripoint������A��B
{
	float xmin, xmax, y, ymin, ymax;
	float x1, x2, y1, y2;
	x1 = A.x;
	y1 = A.y;
	x2 = B.x;
	y2 = B.y;

	float a = (y2 - y1) / (x2 - x1);
	float b = y1 - x1 * a;

	if (A.x > B.x)
	{
		xmin = x2;
		ymin = x2;
		xmax = x1;
		ymax = y1;
	}
	else
	{
		xmin = x1;
		ymin = y1;
		xmax = x2;
		ymax = y2;

	}
	//p = { xmin , ymin };
	//coord.push_back(p);
	for (; xmin < xmax;)
	{
		y = a * xmin + b;
		p = { xmin , y };
		tripoint.push_back(p);
		xmin++;

	}
	p = { xmax , ymax };
	tripoint.push_back(p);

}

point mindistance(point a, point b, point c) //��b �� c�з�����a ����ĵ�
{
	int dis1 = pow((a.x - b.x), 2) + pow((a.y - b.y), 2);
	int dis2 = pow((a.x - c.x), 2) + pow((a.y - c.y), 2);
	if (dis1 < dis2)return b;
	else return c;
}

double distance(point a, point b)
{
	return sqrt(pow((a.x - b.x), 2) + pow((a.y - b.y), 2));
}

//void triangle(int i)  //��i�㣬��1��ʼ����i��wһ�������A��B��C��λ��  1��
//{
//
//	float big_a = 17.0;
//	float small_a = 4.0;
//
//	float s = 0.1;//����,Ҫ��������ߵĿ�ȣ�0<s<=width/2��
//
//	float distance = (big_a - small_a) / (2 * sqrt(3));
//
//	while (i  > 4 * distance / s)
//	{
//		i = i - int(4 * distance / s);
//	}
//	int condition;
//	if (i  < 1 * distance / s) condition = 1;
//
//	else if (i  < 2 * distance / s) condition = 2;
//
//	else if (i  < 3 * distance / s) condition = 3;
//
//	else if (i  < 4 * distance / s) condition = 4;
//
//	switch (condition)
//	{
//	case 1:
//	{
//		for (int m = 0; m < 150 / big_a*sqrt(3)/2; m++)
//		{
//			for (int n = 0; n < 200 / big_a; n++)
//			{
//				point E = { m*0.5*big_a + n * big_a + 0 ,                                  m*0.5*sqrt(3)*big_a + 0 };
//				point B = { m*0.5*big_a + n * big_a + s * sqrt(3)*i  ,                       m*0.5*sqrt(3)*big_a + s * i };
//
//				point F = { m*0.5*big_a + n * big_a + big_a  ,                             m*0.5*sqrt(3)*big_a + 0 };
//				point C = { m*0.5*big_a + n * big_a + big_a - s * sqrt(3)*i ,              m*0.5*sqrt(3)*big_a + s * i };
//
//				point D = { m*0.5*big_a + n * big_a + big_a / 2 ,                          m*0.5*sqrt(3)*big_a + big_a * sqrt(3) / 2 };
//				point A = { m*0.5*big_a + n * big_a + big_a / 2 ,                          m*0.5*sqrt(3)*big_a + big_a * sqrt(3) / 2 - 2 * s * i };
//
//				tripoint.push_back(D); tripoint.push_back(A);
//				tripoint.push_back(A); tripoint.push_back(B);
//				tripoint.push_back(B); tripoint.push_back(E);
//				tripoint.push_back(B); tripoint.push_back(C);
//				tripoint.push_back(C); tripoint.push_back(F);
//				tripoint.push_back(C); tripoint.push_back(A);
//
//			}
//		}
//		break;
//	}
//
//	case 2:
//	{
//		i = i - int(distance / s);
//
//		for (int m = 0; m < 150 / big_a * sqrt(3) / 2; m++)
//		{
//			for (int n = 0; n < 200 / big_a; n++)
//			{
//				point E = { m*0.5*big_a + n * big_a + 0 ,                                  m*0.5*sqrt(3)*big_a + 0 };
//				point B = { m*0.5*big_a + n * big_a + (distance - s * i) * sqrt(3) ,             m*0.5*sqrt(3)*big_a + distance - s * i };
//
//				point F = { m*0.5*big_a + n * big_a + big_a  ,                             m*0.5*sqrt(3)*big_a + 0 };
//				point C = { m*0.5*big_a + n * big_a + big_a - (distance - s * i) * sqrt(3),m*0.5*sqrt(3)*big_a + distance - s * i };
//
//				point D = { m*0.5*big_a + n * big_a + big_a / 2 ,                          m*0.5*sqrt(3)*big_a + big_a * sqrt(3) / 2 };
//				point A = { m*0.5*big_a + n * big_a + big_a / 2 ,                          m*0.5*sqrt(3)*big_a + distance + small_a * sqrt(3) / 2 + 2 * s*i };
//
//				tripoint.push_back(D); tripoint.push_back(A);
//				tripoint.push_back(A); tripoint.push_back(B);
//				tripoint.push_back(B); tripoint.push_back(E);
//				tripoint.push_back(B); tripoint.push_back(C);
//				tripoint.push_back(C); tripoint.push_back(F);
//				tripoint.push_back(C); tripoint.push_back(A);
//			}
//		}
//		break;
//	}
//
//	case 3://��1��ȣ������겻��������ȡ��
//	{
//		i = i - int(2 * distance / s);
//
//		for (int m = 0; m < 150 / big_a * sqrt(3) / 2; m++)
//		{
//			for (int n = 0; n < 200 / big_a; n++)
//			{
//				point E = { m*0.5*big_a + n * big_a + 0 ,                                  m*0.5*sqrt(3)*big_a + 0 };
//				point B = { m*0.5*big_a + n * big_a + s * sqrt(3)*i  ,                     m*0.5*sqrt(3)*big_a - s * i };
//
//				point F = { m*0.5*big_a + n * big_a + big_a  ,                             m*0.5*sqrt(3)*big_a + 0 };
//				point C = { m*0.5*big_a + n * big_a + big_a - s * sqrt(3)*i ,              m*0.5*sqrt(3)*big_a - s * i };
//
//				point D = { m*0.5*big_a + n * big_a + big_a / 2 ,                          m*0.5*sqrt(3)*big_a - big_a * sqrt(3) / 2 };
//				point A = { m*0.5*big_a + n * big_a + big_a / 2 ,                          m*0.5*sqrt(3)*big_a - big_a * sqrt(3) / 2 + 2 * s * i };
//
//				tripoint.push_back(D); tripoint.push_back(A);
//				tripoint.push_back(A); tripoint.push_back(B);
//				tripoint.push_back(B); tripoint.push_back(E);
//				tripoint.push_back(B); tripoint.push_back(C);
//				tripoint.push_back(C); tripoint.push_back(F);
//				tripoint.push_back(C); tripoint.push_back(A);
//
//			}
//		}
//		break;
//	}
//
//	case 4://��2��ȣ������겻�䣬������ȡ��
//	{
//		i = i - int(3 * distance / s);
//
//		for (int m = 0; m < 150 / big_a * sqrt(3) / 2; m++)
//		{
//			for (int n = 0; n < 200 / big_a; n++)
//			{
//				point E = { m*0.5*big_a + n * big_a + 0 ,                                  m*0.5*sqrt(3)*big_a + 0 };
//				point B = { m*0.5*big_a + n * big_a + (distance - s * i) * sqrt(3) ,       m*0.5*sqrt(3)*big_a - distance + s * i };
//
//				point F = { m*0.5*big_a + n * big_a + big_a  ,                             m*0.5*sqrt(3)*big_a + 0 };
//				point C = { m*0.5*big_a + n * big_a + big_a - (distance - s * i) * sqrt(3),m*0.5*sqrt(3)*big_a - distance + s * i };
//
//				point D = { m*0.5*big_a + n * big_a + big_a / 2 ,                          m*0.5*sqrt(3)*big_a - big_a * sqrt(3) / 2 };
//				point A = { m*0.5*big_a + n * big_a + big_a / 2 ,                          m*0.5*sqrt(3)*big_a - distance - small_a * sqrt(3) / 2 - 2 * s*i };
//
//				tripoint.push_back(D); tripoint.push_back(A);
//				tripoint.push_back(A); tripoint.push_back(B);
//				tripoint.push_back(B); tripoint.push_back(E);
//				tripoint.push_back(B); tripoint.push_back(C);
//				tripoint.push_back(C); tripoint.push_back(F);
//				tripoint.push_back(C); tripoint.push_back(A);
//			}
//		}
//		break;
//	}
//
//	default: printf("������");
//		break;
//	}
//}

void triangle(int i)  //û��
{

	float big_a = 17.0;
	float small_a = 6.0;

	float s = 0.1;//����,Ҫ��������ߵĿ�ȣ�0<s<=width/2��

	float distance = (big_a - small_a) / (2 * sqrt(3));

	while (i > 4 * distance / s)
	{
		i = i - int(4 * distance / s);
	}
	int condition;
	if (i < 1 * distance / s) condition = 1;

	else if (i < 2 * distance / s) condition = 2;

	else if (i < 3 * distance / s) condition = 3;

	else if (i < 4 * distance / s) condition = 4;

	switch (condition)
	{
	case 1:
	{
		for (int m = 0; m < 150*2 / big_a * sqrt(3) ; m++)
		{
			for (int n = 0; n < 200 / big_a; n++) //��������꣬���Դ�߳������Ǻ����ܹ��м���������
			{
				point E = { m*0.5*big_a + n * big_a + 0 ,                                  m*0.5*sqrt(3)*big_a + 0 };
				point B = { m*0.5*big_a + n * big_a + s * sqrt(3)*i  ,                       m*0.5*sqrt(3)*big_a + s * i };

				point F = { m*0.5*big_a + n * big_a + big_a  ,                             m*0.5*sqrt(3)*big_a + 0 };
				point C = { m*0.5*big_a + n * big_a + big_a - s * sqrt(3)*i ,              m*0.5*sqrt(3)*big_a + s * i };

				point D = { m*0.5*big_a + n * big_a + big_a / 2 ,                          m*0.5*sqrt(3)*big_a + big_a * sqrt(3) / 2 };
				point A = { m*0.5*big_a + n * big_a + big_a / 2 ,                          m*0.5*sqrt(3)*big_a + big_a * sqrt(3) / 2 - 2 * s * i };


				tripoint.push_back(A); tripoint.push_back(B); 
				tripoint.push_back(B); tripoint.push_back(C);
				tripoint.push_back(C); tripoint.push_back(A);


			}
		}
		break;
	}

	case 2:
	{
		i = i - int(distance / s);

		for (int m = 0; m < 150 / big_a * sqrt(3) / 2; m++)
		{
			for (int n = 0; n < 200 / big_a; n++)
			{
				point E = { m*0.5*big_a + n * big_a + 0 ,                                  m*0.5*sqrt(3)*big_a + 0 };
				point B = { m*0.5*big_a + n * big_a + (distance - s * i) * sqrt(3) ,             m*0.5*sqrt(3)*big_a + distance - s * i };

				point F = { m*0.5*big_a + n * big_a + big_a  ,                             m*0.5*sqrt(3)*big_a + 0 };
				point C = { m*0.5*big_a + n * big_a + big_a - (distance - s * i) * sqrt(3),m*0.5*sqrt(3)*big_a + distance - s * i };

				point D = { m*0.5*big_a + n * big_a + big_a / 2 ,                          m*0.5*sqrt(3)*big_a + big_a * sqrt(3) / 2 };
				point A = { m*0.5*big_a + n * big_a + big_a / 2 ,                          m*0.5*sqrt(3)*big_a + distance + small_a * sqrt(3) / 2 + 2 * s*i };


				tripoint.push_back(A); tripoint.push_back(B);
				tripoint.push_back(B); tripoint.push_back(C);
				tripoint.push_back(C); tripoint.push_back(A);
			}
		}
		break;
	}

	case 3://��1��ȣ������겻��������ȡ��
	{
		i = i - int(2 * distance / s);

		for (int m = 0; m < 150 / big_a * sqrt(3) / 2; m++)
		{
			for (int n = 0; n < 200 / big_a; n++)
			{
				point E = { m*0.5*big_a + n * big_a + 0 ,                                  m*0.5*sqrt(3)*big_a + 0 };
				point B = { m*0.5*big_a + n * big_a + s * sqrt(3)*i  ,                     m*0.5*sqrt(3)*big_a - s * i };

				point F = { m*0.5*big_a + n * big_a + big_a  ,                             m*0.5*sqrt(3)*big_a + 0 };
				point C = { m*0.5*big_a + n * big_a + big_a - s * sqrt(3)*i ,              m*0.5*sqrt(3)*big_a - s * i };

				point D = { m*0.5*big_a + n * big_a + big_a / 2 ,                          m*0.5*sqrt(3)*big_a - big_a * sqrt(3) / 2 };
				point A = { m*0.5*big_a + n * big_a + big_a / 2 ,                          m*0.5*sqrt(3)*big_a - big_a * sqrt(3) / 2 + 2 * s * i };


				tripoint.push_back(A); tripoint.push_back(B);
				tripoint.push_back(B); tripoint.push_back(C);
				tripoint.push_back(C); tripoint.push_back(A);
			}
		}
		break;
	}

	case 4://��2��ȣ������겻�䣬������ȡ��
	{
		i = i - int(3 * distance / s);

		for (int m = 0; m < 150 / big_a * sqrt(3) / 2; m++)
		{
			for (int n = 0; n < 200 / big_a; n++)
			{
				point E = { m*0.5*big_a + n * big_a + 0 ,                                  m*0.5*sqrt(3)*big_a + 0 };
				point B = { m*0.5*big_a + n * big_a + (distance - s * i) * sqrt(3) ,       m*0.5*sqrt(3)*big_a - distance + s * i };

				point F = { m*0.5*big_a + n * big_a + big_a  ,                             m*0.5*sqrt(3)*big_a + 0 };
				point C = { m*0.5*big_a + n * big_a + big_a - (distance - s * i) * sqrt(3),m*0.5*sqrt(3)*big_a - distance + s * i };

				point D = { m*0.5*big_a + n * big_a + big_a / 2 ,                          m*0.5*sqrt(3)*big_a - big_a * sqrt(3) / 2 };
				point A = { m*0.5*big_a + n * big_a + big_a / 2 ,                          m*0.5*sqrt(3)*big_a - distance - small_a * sqrt(3) / 2 - 2 * s*i };


				tripoint.push_back(A); tripoint.push_back(B);
				tripoint.push_back(B); tripoint.push_back(C);
				tripoint.push_back(C); tripoint.push_back(A);
			}
		}
		break;
	}

	default: printf("������");
		break;
	}
}

//void triangle(int i)  //û��,��ͼ�����˿���⣬�����Ȼ����еĺᣬ�ٻ����е�د���ٻ����е���
//{
//
//	float big_a = 17.0;
//	float small_a = 6.0;
//
//	int m = 0, n = 0;
//	int M, N;
//	N = 10;
//	M = N;
//
//	float s = 0.1;//����,Ҫ��������ߵĿ�ȣ�0<s<=width/2��
//
//	float distance = (big_a - small_a) / (2 * sqrt(3));
//
//	while (i > 4 * distance / s)
//	{
//		i = i - int(4 * distance / s);
//	}
//	int condition;
//	if (i < 1 * distance / s) condition = 1;
//
//	else if (i < 2 * distance / s) condition = 2;
//
//	else if (i < 3 * distance / s) condition = 3;
//
//	else if (i < 4 * distance / s) condition = 4;
//
//	switch (condition)
//	{
//	case 1:
//	{
//
//		for (; m <= M; m++) //һ
//		{
//			if (m % 2 == 0)
//			{
//				for (int n = 0; n <= N; n++)
//				{
//					point B = { m*0.5*big_a + n * big_a + s * sqrt(3)*i  ,                       m*0.5*sqrt(3)*big_a + s * i };
//					point C = { m*0.5*big_a + n * big_a + big_a - s * sqrt(3)*i ,              m*0.5*sqrt(3)*big_a + s * i };
//					tripoint.push_back(B); tripoint.push_back(C);
//				}
//			}
//			else
//			{
//				for (int n = N; n >= 0; n--)
//				{
//					point B = { m*0.5*big_a + n * big_a + s * sqrt(3)*i  ,                       m*0.5*sqrt(3)*big_a + s * i };
//					point C = { m*0.5*big_a + n * big_a + big_a - s * sqrt(3)*i ,              m*0.5*sqrt(3)*big_a + s * i };
//					tripoint.push_back(C); tripoint.push_back(B);
//				}
//			}
//		}
//		m = 0, n = 0;
//		for (; n <= N; n++) //د
//		{
//			if (n % 2 == 0)
//			{
//				for (int m = 0; m <= M; m++)
//				{
//					point A = { m*0.5*big_a + n * big_a + big_a / 2 ,                          m*0.5*sqrt(3)*big_a + big_a * sqrt(3) / 2 - 2 * s * i };
//					point B = { m*0.5*big_a + n * big_a + s * sqrt(3)*i  ,                       m*0.5*sqrt(3)*big_a + s * i };
//					tripoint.push_back(B); tripoint.push_back(A);
//				}
//			}
//			else
//			{
//				for (int m = M ; m >= 0; m--)
//				{
//					point A = { m*0.5*big_a + n * big_a + big_a / 2 ,                          m*0.5*sqrt(3)*big_a + big_a * sqrt(3) / 2 - 2 * s * i };
//					point B = { m*0.5*big_a + n * big_a + s * sqrt(3)*i  ,                       m*0.5*sqrt(3)*big_a + s * i };
//					tripoint.push_back(A); tripoint.push_back(B);
//				}
//			}
//		}
//	
//		for (int j = 1; j <= 2*M+1; j++)//������б�ű�����
//		{
//			if (j <= M + 1)
//			{
//				if (j % 2 == 0)
//				{
//					n = N - j + 1;
//					m = M;
//					goto even1;
//				}
//				else
//				{
//					n = N;
//					m = M - j + 1;
//					goto odd1;
//				}
//			}
//			else
//			{
//				if (j % 2 == 0)
//				{
//					n = 0;
//					m = 2 * M + 1 - j;
//					goto even1;
//				}
//				else
//				{
//					n = 2 * M + 1 - j;
//					m = 0;
//					goto odd1;
//				}
//			}
//
//		even1:
//			do//������Ҫ���࣬��ǰ������Ӧ��
//			{
//				point A = { m*0.5*big_a + n * big_a + big_a / 2 ,                          m*0.5*sqrt(3)*big_a + big_a * sqrt(3) / 2 - 2 * s * i };
//				point C = { m*0.5*big_a + n * big_a + big_a - s * sqrt(3)*i ,              m*0.5*sqrt(3)*big_a + s * i };
//				tripoint.push_back(A); tripoint.push_back(C);
//				
//				m--;
//				n++;
//			} while (m >= 0 && n <= N);
//			continue;
//		odd1:
//			do
//			{
//				point A = { m*0.5*big_a + n * big_a + big_a / 2 ,                          m*0.5*sqrt(3)*big_a + big_a * sqrt(3) / 2 - 2 * s * i };
//				point C = { m*0.5*big_a + n * big_a + big_a - s * sqrt(3)*i ,              m*0.5*sqrt(3)*big_a + s * i };
//				tripoint.push_back(C); tripoint.push_back(A);
//				
//				m++;
//				n--;
//			} while (m <= M && n >= 0);
//
//		}
//
//		break;
//	}
//
//	case 2:
//	{
//		i = i - int(distance / s);
//
//		/*for (int m = 0; m < 150 / big_a * sqrt(3) / 2; m++)
//		{
//			for (int n = 0; n < 200 / big_a; n++)
//			{
//				point E = { m*0.5*big_a + n * big_a + 0 ,                                  m*0.5*sqrt(3)*big_a + 0 };
//				point B = { m*0.5*big_a + n * big_a + (distance - s * i) * sqrt(3) ,             m*0.5*sqrt(3)*big_a + distance - s * i };
//
//				point F = { m*0.5*big_a + n * big_a + big_a  ,                             m*0.5*sqrt(3)*big_a + 0 };
//				point C = { m*0.5*big_a + n * big_a + big_a - (distance - s * i) * sqrt(3),m*0.5*sqrt(3)*big_a + distance - s * i };
//
//				point D = { m*0.5*big_a + n * big_a + big_a / 2 ,                          m*0.5*sqrt(3)*big_a + big_a * sqrt(3) / 2 };
//				point A = { m*0.5*big_a + n * big_a + big_a / 2 ,                          m*0.5*sqrt(3)*big_a + distance + small_a * sqrt(3) / 2 + 2 * s*i };
//
//
//				tripoint.push_back(A); tripoint.push_back(B);
//				tripoint.push_back(B); tripoint.push_back(C);
//				tripoint.push_back(C); tripoint.push_back(A);
//			}
//		}*/
//
//
//		for (; m <= M; m++) //һ
//		{
//			if (m % 2 == 0)
//			{
//				for (int n = 0; n <= N; n++)
//				{
//					point B = { m*0.5*big_a + n * big_a + (distance - s * i) * sqrt(3) ,             m*0.5*sqrt(3)*big_a + distance - s * i };
//					point C = { m*0.5*big_a + n * big_a + big_a - (distance - s * i) * sqrt(3),m*0.5*sqrt(3)*big_a + distance - s * i };
//					tripoint.push_back(B); tripoint.push_back(C);
//				}
//			}
//			else
//			{
//				for (int n = N; n >= 0; n--)
//				{
//					point B = { m*0.5*big_a + n * big_a + (distance - s * i) * sqrt(3) ,             m*0.5*sqrt(3)*big_a + distance - s * i };
//					point C = { m*0.5*big_a + n * big_a + big_a - (distance - s * i) * sqrt(3),m*0.5*sqrt(3)*big_a + distance - s * i };
//					tripoint.push_back(C); tripoint.push_back(B);
//				}
//			}
//		}
//		m = 0, n = 0;
//		for (; n <= N; n++) //د
//		{
//			if (n % 2 == 0)
//			{
//				for (int m = 0; m <= M; m++)
//				{
//					point A = { m*0.5*big_a + n * big_a + big_a / 2 ,                          m*0.5*sqrt(3)*big_a + distance + small_a * sqrt(3) / 2 + 2 * s*i };
//					point B = { m*0.5*big_a + n * big_a + (distance - s * i) * sqrt(3) ,             m*0.5*sqrt(3)*big_a + distance - s * i };
//					tripoint.push_back(B); tripoint.push_back(A);
//				}
//			}
//			else
//			{
//				for (int m = M; m >= 0; m--)
//				{
//					point A = { m*0.5*big_a + n * big_a + big_a / 2 ,                          m*0.5*sqrt(3)*big_a + distance + small_a * sqrt(3) / 2 + 2 * s*i };
//					point B = { m*0.5*big_a + n * big_a + (distance - s * i) * sqrt(3) ,             m*0.5*sqrt(3)*big_a + distance - s * i };
//					tripoint.push_back(A); tripoint.push_back(B);
//				}
//			}
//		}
//
//		for (int j = 1; j <= 2 * M + 1; j++)//������б�ű�����
//		{
//			if (j <= M + 1)
//			{
//				if (j % 2 == 0)
//				{
//					n = N - j + 1;
//					m = M;
//					goto even2;
//				}
//				else
//				{
//					n = N;
//					m = M - j + 1;
//					goto odd2;
//				}
//			}
//			else
//			{
//				if (j % 2 == 0)
//				{
//					n = 0;
//					m = 2 * M + 1 - j;
//					goto even2;
//				}
//				else
//				{
//					n = 2 * M + 1 - j;
//					m = 0;
//					goto odd2;
//				}
//			}
//
//		even2:
//			do//������Ҫ���࣬��ǰ������Ӧ��
//			{
//				point A = { m*0.5*big_a + n * big_a + big_a / 2 ,                          m*0.5*sqrt(3)*big_a + distance + small_a * sqrt(3) / 2 + 2 * s*i };
//				point C = { m*0.5*big_a + n * big_a + big_a - (distance - s * i) * sqrt(3),m*0.5*sqrt(3)*big_a + distance - s * i };
//				tripoint.push_back(A); tripoint.push_back(C);
//
//				m--;
//				n++;
//			} while (m >= 0 && n <= N);
//			continue;
//		odd2:
//			do
//			{
//				point A = { m*0.5*big_a + n * big_a + big_a / 2 ,                          m*0.5*sqrt(3)*big_a + distance + small_a * sqrt(3) / 2 + 2 * s*i };
//				point C = { m*0.5*big_a + n * big_a + big_a - (distance - s * i) * sqrt(3),m*0.5*sqrt(3)*big_a + distance - s * i };
//				tripoint.push_back(C); tripoint.push_back(A);
//
//				m++;
//				n--;
//			} while (m <= M && n >= 0);
//
//		}
//		break;
//	}
//
//	case 3://��1��ȣ������겻��������ȡ��
//	{
//		i = i - int(2 * distance / s);
//
//		/*for (int m = 0; m < 150 / big_a * sqrt(3) / 2; m++)
//		{
//			for (int n = 0; n < 200 / big_a; n++)
//			{
//				point E = { m*0.5*big_a + n * big_a + 0 ,                                  m*0.5*sqrt(3)*big_a + 0 };
//				point B = { m*0.5*big_a + n * big_a + s * sqrt(3)*i  ,                     m*0.5*sqrt(3)*big_a - s * i };
//
//				point F = { m*0.5*big_a + n * big_a + big_a  ,                             m*0.5*sqrt(3)*big_a + 0 };
//				point C = { m*0.5*big_a + n * big_a + big_a - s * sqrt(3)*i ,              m*0.5*sqrt(3)*big_a - s * i };
//
//				point D = { m*0.5*big_a + n * big_a + big_a / 2 ,                          m*0.5*sqrt(3)*big_a - big_a * sqrt(3) / 2 };
//				point A = { m*0.5*big_a + n * big_a + big_a / 2 ,                          m*0.5*sqrt(3)*big_a - big_a * sqrt(3) / 2 + 2 * s * i };
//
//
//				tripoint.push_back(A); tripoint.push_back(B);
//				tripoint.push_back(B); tripoint.push_back(C);
//				tripoint.push_back(C); tripoint.push_back(A);
//			}
//		}*/
//
//
//		for (; m <= M; m++) //һ
//		{
//			if (m % 2 == 0)
//			{
//				for (int n = 0; n <= N; n++)
//				{
//					point B = { m*0.5*big_a + n * big_a + s * sqrt(3)*i  ,                     m*0.5*sqrt(3)*big_a - s * i };
//					point C = { m*0.5*big_a + n * big_a + big_a - s * sqrt(3)*i ,              m*0.5*sqrt(3)*big_a - s * i };
//					tripoint.push_back(B); tripoint.push_back(C);
//				}
//			}
//			else
//			{
//				for (int n = N; n >= 0; n--)
//				{
//					point B = { m*0.5*big_a + n * big_a + s * sqrt(3)*i  ,                     m*0.5*sqrt(3)*big_a - s * i };
//					point C = { m*0.5*big_a + n * big_a + big_a - s * sqrt(3)*i ,              m*0.5*sqrt(3)*big_a - s * i };
//					tripoint.push_back(C); tripoint.push_back(B);
//				}
//			}
//		}
//		m = 0, n = 0;
//		for (; n <= N; n++) //د
//		{
//			if (n % 2 == 0)
//			{
//				for (int m = 0; m <= M; m++)
//				{
//					point A = { m*0.5*big_a + n * big_a + big_a / 2 ,                          m*0.5*sqrt(3)*big_a - big_a * sqrt(3) / 2 + 2 * s * i };
//					point C = { m*0.5*big_a + n * big_a + big_a - s * sqrt(3)*i ,              m*0.5*sqrt(3)*big_a - s * i };
//					tripoint.push_back(A); tripoint.push_back(C);
//				}
//			}
//			else
//			{
//				for (int m = M; m >= 0; m--)
//				{
//					point A = { m*0.5*big_a + n * big_a + big_a / 2 ,                          m*0.5*sqrt(3)*big_a - big_a * sqrt(3) / 2 + 2 * s * i };
//					point C = { m*0.5*big_a + n * big_a + big_a - s * sqrt(3)*i ,              m*0.5*sqrt(3)*big_a - s * i };
//					tripoint.push_back(C); tripoint.push_back(A);
//				}
//			}
//		}
//
//		for (int j = 1; j <= 2 * M + 1; j++)//������б�ű�����
//		{
//			if (j <= M + 1)
//			{
//				if (j % 2 == 0)
//				{
//					n = N - j + 1;
//					m = M;
//					goto even3;
//				}
//				else
//				{
//					n = N;
//					m = M - j + 1;
//					goto odd3;
//				}
//			}
//			else
//			{
//				if (j % 2 == 0)
//				{
//					n = 0;
//					m = 2 * M + 1 - j;
//					goto even3;
//				}
//				else
//				{
//					n = 2 * M + 1 - j;
//					m = 0;
//					goto odd3;
//				}
//			}
//
//		even3:
//			do//������Ҫ���࣬��ǰ������Ӧ��
//			{
//				point A = { m*0.5*big_a + n * big_a + big_a / 2 ,                          m*0.5*sqrt(3)*big_a - big_a * sqrt(3) / 2 + 2 * s * i };
//				point B = { m*0.5*big_a + n * big_a + s * sqrt(3)*i  ,                     m*0.5*sqrt(3)*big_a - s * i };
//				tripoint.push_back(B); tripoint.push_back(A);
//
//				m--;
//				n++;
//			} while (m >= 0 && n <= N);
//			continue;
//		odd3:
//			do
//			{
//				point A = { m*0.5*big_a + n * big_a + big_a / 2 ,                          m*0.5*sqrt(3)*big_a - big_a * sqrt(3) / 2 + 2 * s * i };
//				point B = { m*0.5*big_a + n * big_a + s * sqrt(3)*i  ,                     m*0.5*sqrt(3)*big_a - s * i };
//				tripoint.push_back(A); tripoint.push_back(B);
//
//				m++;
//				n--;
//			} while (m <= M && n >= 0);
//
//		}
//
//		break;
//	}
//
//	case 4://��2��ȣ������겻�䣬������ȡ��
//	{
//		i = i - int(3 * distance / s);
//
//		/*for (int m = 0; m < 150 / big_a * sqrt(3) / 2; m++)
//		{
//			for (int n = 0; n < 200 / big_a; n++)
//			{
//				point E = { m*0.5*big_a + n * big_a + 0 ,                                  m*0.5*sqrt(3)*big_a + 0 };
//				point B = { m*0.5*big_a + n * big_a + (distance - s * i) * sqrt(3) ,       m*0.5*sqrt(3)*big_a - distance + s * i };
//
//				point F = { m*0.5*big_a + n * big_a + big_a  ,                             m*0.5*sqrt(3)*big_a + 0 };
//				point C = { m*0.5*big_a + n * big_a + big_a - (distance - s * i) * sqrt(3),m*0.5*sqrt(3)*big_a - distance + s * i };
//
//				point D = { m*0.5*big_a + n * big_a + big_a / 2 ,                          m*0.5*sqrt(3)*big_a - big_a * sqrt(3) / 2 };
//				point A = { m*0.5*big_a + n * big_a + big_a / 2 ,                          m*0.5*sqrt(3)*big_a - distance - small_a * sqrt(3) / 2 - 2 * s*i };
//
//
//				tripoint.push_back(A); tripoint.push_back(B);
//				tripoint.push_back(B); tripoint.push_back(C);
//				tripoint.push_back(C); tripoint.push_back(A);
//			}
//		}*/
//
//
//		for (; m <= M; m++) //һ
//		{
//			if (m % 2 == 0)
//			{
//				for (int n = 0; n <= N; n++)
//				{
//					point B = { m*0.5*big_a + n * big_a + (distance - s * i) * sqrt(3) ,       m*0.5*sqrt(3)*big_a - distance + s * i };
//					point C = { m*0.5*big_a + n * big_a + big_a - (distance - s * i) * sqrt(3),m*0.5*sqrt(3)*big_a - distance + s * i };
//					tripoint.push_back(B); tripoint.push_back(C);
//				}
//			}
//			else
//			{
//				for (int n = N; n >= 0; n--)
//				{
//					point B = { m*0.5*big_a + n * big_a + (distance - s * i) * sqrt(3) ,       m*0.5*sqrt(3)*big_a - distance + s * i };
//					point C = { m*0.5*big_a + n * big_a + big_a - (distance - s * i) * sqrt(3),m*0.5*sqrt(3)*big_a - distance + s * i };
//					tripoint.push_back(C); tripoint.push_back(B);
//				}
//			}
//		}
//		m = 0, n = 0;
//		for (; n <= N; n++) //د
//		{
//			if (n % 2 == 0)
//			{
//				for (int m = 0; m <= M; m++)
//				{
//					point A = { m*0.5*big_a + n * big_a + big_a / 2 ,                          m*0.5*sqrt(3)*big_a - distance - small_a * sqrt(3) / 2 - 2 * s*i };
//					point C = { m*0.5*big_a + n * big_a + big_a - (distance - s * i) * sqrt(3),m*0.5*sqrt(3)*big_a - distance + s * i };
//					tripoint.push_back(A); tripoint.push_back(C);
//				}
//			}
//			else
//			{
//				for (int m = M; m >= 0; m--)
//				{
//					point A = { m*0.5*big_a + n * big_a + big_a / 2 ,                          m*0.5*sqrt(3)*big_a - distance - small_a * sqrt(3) / 2 - 2 * s*i };
//					point C = { m*0.5*big_a + n * big_a + big_a - (distance - s * i) * sqrt(3),m*0.5*sqrt(3)*big_a - distance + s * i };
//					tripoint.push_back(C); tripoint.push_back(A);
//				}
//			}
//		}
//
//		for (int j = 1; j <= 2 * M + 1; j++)//������б�ű�����
//		{
//			if (j <= M + 1)
//			{
//				if (j % 2 == 0)
//				{
//					n = N - j + 1;
//					m = M;
//					goto even4;
//				}
//				else
//				{
//					n = N;
//					m = M - j + 1;
//					goto odd4;
//				}
//			}
//			else
//			{
//				if (j % 2 == 0)
//				{
//					n = 0;
//					m = 2 * M + 1 - j;
//					goto even4;
//				}
//				else
//				{
//					n = 2 * M + 1 - j;
//					m = 0;
//					goto odd4;
//				}
//			}
//
//		even4:
//			do//������Ҫ���࣬��ǰ������Ӧ��
//			{
//				point A = { m*0.5*big_a + n * big_a + big_a / 2 ,                          m*0.5*sqrt(3)*big_a - distance - small_a * sqrt(3) / 2 - 2 * s*i };
//				point B = { m*0.5*big_a + n * big_a + (distance - s * i) * sqrt(3) ,       m*0.5*sqrt(3)*big_a - distance + s * i };
//				tripoint.push_back(B); tripoint.push_back(A);
//
//				m--;
//				n++;
//			} while (m >= 0 && n <= N);
//			continue;
//		odd4:
//			do
//			{
//				point A = { m*0.5*big_a + n * big_a + big_a / 2 ,                          m*0.5*sqrt(3)*big_a - distance - small_a * sqrt(3) / 2 - 2 * s*i };
//				point B = { m*0.5*big_a + n * big_a + (distance - s * i) * sqrt(3) ,       m*0.5*sqrt(3)*big_a - distance + s * i };
//				tripoint.push_back(A); tripoint.push_back(B);
//
//				m++;
//				n--;
//			} while (m <= M && n >= 0);
//
//		}
//		break;
//	}
//
//	default: printf("������");
//		
//	}
//}

//void triangle(int i)  //��ͼ��ͨ���еĶ�
//{
//
//	float big_a = 10.0;
//	float small_a = 6.0;
//
//	float s = 0.1;//����,Ҫ��������ߵĿ�ȣ�0<s<=width/2��
//
//	float distance = (big_a - small_a) / (2 * sqrt(3));
//
//	while (i > 6 * distance / s)
//	{
//		i = i - int(6 * distance / s);
//	}
//	int condition;
//	if (i < 1 * distance / s) condition = 1;
//
//	else if (i < 2 * distance / s) condition = 2;
//
//	else if (i < 3 * distance / s) condition = 3;
//
//	else if (i < 4 * distance / s) condition = 4;
//
//	else if (i < 5 * distance / s) condition = 5;
//
//	else if (i < 6 * distance / s) condition = 6;
//
//	switch (condition)
//	{
//	case 1:
//	{
//		for (int m = 0; m < 150 / big_a * sqrt(3) / 2; m++)
//		{
//			for (int n = 0; n < 200 / big_a; n++)
//			{
//				point E = { m*0.5*big_a + n * big_a + 0 ,                                  m*0.5*sqrt(3)*big_a + 0 };
//				point B = { m*0.5*big_a + n * big_a + s * sqrt(3)*i  ,                       m*0.5*sqrt(3)*big_a + s * i };
//
//				point F = { m*0.5*big_a + n * big_a + big_a  ,                             m*0.5*sqrt(3)*big_a + 0 };
//				point C = { m*0.5*big_a + n * big_a + big_a - s * sqrt(3)*i ,              m*0.5*sqrt(3)*big_a + s * i };
//
//				point D = { m*0.5*big_a + n * big_a + big_a / 2 ,                          m*0.5*sqrt(3)*big_a + big_a * sqrt(3) / 2 };
//				point A = { m*0.5*big_a + n * big_a + big_a / 2 ,                          m*0.5*sqrt(3)*big_a + big_a * sqrt(3) / 2 - 2 * s * i };
//
//
//				tripoint.push_back(A); tripoint.push_back(B);
//				tripoint.push_back(B); tripoint.push_back(C);
//				tripoint.push_back(C); tripoint.push_back(A);
//
//
//			}
//		}
//		break;
//	}
//
//	case 2:
//	{
//		i = i - int(distance / s);
//
//		for (int m = 0; m < 150 / big_a * sqrt(3) / 2; m++)
//		{
//			for (int n = 0; n < 200 / big_a; n++)
//			{
//				point E = { m*0.5*big_a + n * big_a + 0 ,                                  m*0.5*sqrt(3)*big_a + 0 };
//				point B = { m*0.5*big_a + n * big_a + (distance - s * i) * sqrt(3) ,       m*0.5*sqrt(3)*big_a + distance - s * i };
//
//				point F = { m*0.5*big_a + n * big_a + big_a  ,                             m*0.5*sqrt(3)*big_a + 0 };
//				point C = { m*0.5*big_a + n * big_a + big_a - (distance - s * i) * sqrt(3),m*0.5*sqrt(3)*big_a + distance - s * i };
//
//				point D = { m*0.5*big_a + n * big_a + big_a / 2 ,                          m*0.5*sqrt(3)*big_a + big_a * sqrt(3) * 0.5 };
//				point A = { m*0.5*big_a + n * big_a + big_a / 2 ,                          m*0.5*sqrt(3)*big_a + distance + small_a * sqrt(3) * 0.5 + 2 * s*i };
//
//
//				tripoint.push_back(A); tripoint.push_back(B);
//				tripoint.push_back(B); tripoint.push_back(C);
//				tripoint.push_back(C); tripoint.push_back(A);
//			}
//		}
//		break;
//	}
//
//	case 3://��1��ȣ������겻��������ȡ��
//	{
//		i = i - int(2 * distance / s);
//
//		for (int m = 0; m < 150 / big_a * sqrt(3) / 2; m++)
//		{
//			for (int n = 0; n < 200 / big_a; n++)
//			{
//				point E = { m*0.5*big_a + n * big_a + 0 ,                                  m*0.5*sqrt(3)*big_a + 0 };
//				point B = { m*0.5*big_a + n * big_a + s * sqrt(3)*i  ,                     m*0.5*sqrt(3)*big_a - s * i };
//
//				point F = { m*0.5*big_a + n * big_a + big_a  ,                             m*0.5*sqrt(3)*big_a + 0 };
//				point C = { m*0.5*big_a + n * big_a + big_a - s * sqrt(3)*i ,              m*0.5*sqrt(3)*big_a - s * i };
//
//				point D = { m*0.5*big_a + n * big_a + big_a / 2 ,                          m*0.5*sqrt(3)*big_a - big_a * sqrt(3) / 2 };
//				point A = { m*0.5*big_a + n * big_a + big_a / 2 ,                          m*0.5*sqrt(3)*big_a - big_a * sqrt(3) / 2 + 2 * s * i };
//
//
//				tripoint.push_back(A); tripoint.push_back(B);
//				tripoint.push_back(B); tripoint.push_back(C);
//				tripoint.push_back(C); tripoint.push_back(A);
//			}
//		}
//		break;
//	}
//
//	case 4://��2��ȣ������겻�䣬������ȡ��
//	{
//		i = i - int(3 * distance / s);
//
//		for (int m = 0; m < 150 / big_a * sqrt(3) / 2; m++)
//		{
//			for (int n = 0; n < 200 / big_a; n++)
//			{
//				point E = { m*0.5*big_a + n * big_a + 0 ,                                  m*0.5*sqrt(3)*big_a + 0 };
//				point B = { m*0.5*big_a + n * big_a + (distance - s * i) * sqrt(3) ,       m*0.5*sqrt(3)*big_a - distance + s * i };
//
//				point F = { m*0.5*big_a + n * big_a + big_a  ,                             m*0.5*sqrt(3)*big_a + 0 };
//				point C = { m*0.5*big_a + n * big_a + big_a - (distance - s * i) * sqrt(3),m*0.5*sqrt(3)*big_a - distance + s * i };
//
//				point D = { m*0.5*big_a + n * big_a + big_a / 2 ,                          m*0.5*sqrt(3)*big_a - big_a * sqrt(3) / 2 };
//				point A = { m*0.5*big_a + n * big_a + big_a / 2 ,                          m*0.5*sqrt(3)*big_a - distance - small_a * sqrt(3) / 2 - 2 * s*i };
//
//
//				tripoint.push_back(A); tripoint.push_back(B);
//				tripoint.push_back(B); tripoint.push_back(C);
//				tripoint.push_back(C); tripoint.push_back(A);
//			}
//		}
//		break;
//	}
//
//	case 5:
//	{
//		i = i - int(4 * distance / s);
//
//		for (int m = 0; m < 150 / big_a * sqrt(3) / 2; m++)
//		{
//			for (int n = 0; n < 200 / big_a; n++)
//			{
//				point E = { n * big_a + s * sqrt(3)*i ,                     m*sqrt(3)*big_a + sqrt(3)*big_a - s * i };
//				point B = { n * big_a + s * sqrt(3)*i ,                     m*sqrt(3)*big_a + s * i };
//
//				point F = { n * big_a + big_a - s * sqrt(3)*i ,             m*sqrt(3)*big_a + sqrt(3)*big_a - s * i };
//				point C = { n * big_a + big_a - s * sqrt(3)*i ,             m*sqrt(3)*big_a + s * i };
//
//				point D = { n * big_a + big_a * 0.5 ,                       m*sqrt(3)*big_a + big_a * sqrt(3) * 0.5 + 2 * s * i };
//				point A = { n * big_a + big_a * 0.5 ,                       m*sqrt(3)*big_a + big_a * sqrt(3) * 0.5 - 2 * s * i };
//
//
//
//				//point K = { n * big_a + big_a * 0.5 + s * sqrt(3)*i,        m*sqrt(3)*big_a + big_a * sqrt(3) * 0.5 + s * i };
//				//point H = { n * big_a + big_a * 0.5 + s * sqrt(3)*i,        m*sqrt(3)*big_a + big_a * sqrt(3) * 0.5 - s * i };
//
//				//point J = { n * big_a + big_a ,                             m*sqrt(3)*big_a + sqrt(3)*big_a - 2 * s * i };
//				//point G = { n * big_a + big_a ,                             m*sqrt(3)*big_a + 2 * s * i };
//
//				//point L = { n * big_a + big_a * 1.5 - s * sqrt(3)*i ,       m*sqrt(3)*big_a + big_a * sqrt(3) * 0.5 + s * i };
//				//point I = { n * big_a + big_a * 1.5 - s * sqrt(3)*i,        m*sqrt(3)*big_a + big_a * sqrt(3) * 0.5 - s * i };
//
//
//
//
//				tripoint.push_back(D); tripoint.push_back(E);
//				tripoint.push_back(E); tripoint.push_back(F);
//				tripoint.push_back(F); tripoint.push_back(D);
//				tripoint.push_back(A); tripoint.push_back(B);
//				tripoint.push_back(B); tripoint.push_back(C);
//				tripoint.push_back(C); tripoint.push_back(A);
//
//				/*tripoint.push_back(K); tripoint.push_back(J);
//				tripoint.push_back(J); tripoint.push_back(L);
//				tripoint.push_back(L); tripoint.push_back(K);
//				tripoint.push_back(H); tripoint.push_back(G);
//				tripoint.push_back(G); tripoint.push_back(I);
//				tripoint.push_back(I); tripoint.push_back(H);
//*/
//
//			}
//		}
//		break;
//	}
//
//	case 6:
//	{
//		i = i - int(5 * distance / s);
//
//		for (int m = 0; m < 150 / big_a * sqrt(3) / 2; m++)
//		{
//			for (int n = 0; n < 200 / big_a; n++)
//			{
//				point E = { n * big_a + (distance - s * i) * sqrt(3) ,           m*sqrt(3)*big_a + sqrt(3)*big_a - distance + s * i };
//				point B = { n * big_a + (distance - s * i) * sqrt(3) ,           m*sqrt(3)*big_a + distance - s * i };
//
//				point F = { n * big_a + big_a - (distance - s * i) * sqrt(3) ,   m*sqrt(3)*big_a + sqrt(3)*big_a - distance + s * i };
//				point C = { n * big_a + big_a - (distance - s * i) * sqrt(3) ,   m*sqrt(3)*big_a + distance - s * i };
//
//				point D = { n * big_a + big_a / 2 ,                              m*sqrt(3)*big_a + 5*distance + small_a * sqrt(3) * 0.5 - 2 * s*i };
//				point A = { n * big_a + big_a / 2 ,                              m*sqrt(3)*big_a + distance + small_a * sqrt(3) * 0.5 + 2 * s*i };
//
//				/*point K = { n * big_a + big_a - 0.5*small_a - s * i * sqrt(3),   m*sqrt(3)*big_a + big_a * sqrt(3) * 0.5 +distance - s * i };
//				point H = { n * big_a + big_a - 0.5*small_a - s * i * sqrt(3),   m*sqrt(3)*big_a + small_a * sqrt(3) * 0.5 + 2*distance + s * i };
//
//				point J = { n * big_a + big_a ,                                  m*sqrt(3)*big_a + sqrt(3)*big_a - 2 * distance + 2 * s * i };
//				point G = { n * big_a + big_a ,                                  m*sqrt(3)*big_a + 2 * distance - 2 * s * i };
//
//				point L = { n * big_a + big_a + 0.5*small_a + s * sqrt(3)*i,     m*sqrt(3)*big_a + big_a * sqrt(3) * 0.5 + distance - s * i };
//				point I = { n * big_a + big_a + 0.5*small_a + s * sqrt(3)*i,     m*sqrt(3)*big_a + small_a * sqrt(3) * 0.5 + 2 * distance + s * i };
//*/
//
//				tripoint.push_back(D); tripoint.push_back(E);
//				tripoint.push_back(E); tripoint.push_back(F);
//				tripoint.push_back(F); tripoint.push_back(D);
//				tripoint.push_back(A); tripoint.push_back(B);
//				tripoint.push_back(B); tripoint.push_back(C);
//				tripoint.push_back(C); tripoint.push_back(A);
//
//				/*tripoint.push_back(K); tripoint.push_back(J);
//				tripoint.push_back(J); tripoint.push_back(L);
//				tripoint.push_back(L); tripoint.push_back(K);
//				tripoint.push_back(H); tripoint.push_back(G);
//				tripoint.push_back(G); tripoint.push_back(I);
//				tripoint.push_back(I); tripoint.push_back(H);*/
//
//			}
//		}
//		break;
//	}
//
//	default: printf("������");
//		break;
//	}
//}

void intersection()           //�����������Ҫ��Ҫ�仯
{
	for (int i = 0; i < model.size(); i++)
	{
		tripoint.clear(); //����һ��ı仯�������
		interpoint.clear(); //����һ��Ľ��������
		polypoint.clear();//����һ��������������
		triangle(i); //������i��ı仯�ߴ���tripoint        0���������ǣ�48���������Σ� i�仯��

		polypoint = model[i];//�ѵ�i��������߸�polypoint

		for (int m = 0; m < tripoint.size(); m += 2)  //�����ÿһ���߶Σ���β�����߶ζ˵��ظ�������������
		{
			interpoint.clear(); //��һ���߶εĽ��������
			for (int n = 0; n < polypoint.size(); n++) //����ε�ÿһ���߶Σ���β�����߶ζ˵㲻�ظ���������һ��
			{
				inter(tripoint[m], tripoint[m + 1], polypoint[n], polypoint[(n + 1 + polypoint.size()) % polypoint.size()]);
			}
			int a = InOrOutPolygon(tripoint[m]);
			int b = InOrOutPolygon(tripoint[m + 1]);  //�ж��ǲ����ڵ�,0:�� 1:��

			if (interpoint.size() == 0)   // �޽���
			{
				if ((a == 0) && (b == 0)) //���ڵ㣬��
				{
					continue;
				}
				if ((a == 1) || (b == 1)) //�����ڵ㣬ȡ
				{
					new_tripoint.push_back(tripoint[m]);
					new_tripoint.push_back(tripoint[m + 1]);
					continue;
				}
			}

			if (interpoint.size() == 1)   //��һ����
			{
				if ((a == 0) && (b == 0))   //��һ���������ڵ㣬�� �����У�
				{
					continue;
				}
				if ((a == 1) && (b == 0))//��һ��������һ�ڵ㣬ȡ�ڵ㵽���� ���ཻ��
				{
					new_tripoint.push_back(tripoint[m]);
					new_tripoint.push_back(interpoint[0]);
					continue;
				}
				if ((a == 0) && (b == 1))//��һ��������һ�ڵ㣬ȡ�ڵ㵽���� ���ཻ��
				{
					new_tripoint.push_back(interpoint[0]);
					new_tripoint.push_back(tripoint[m + 1]);
					continue;
				}
				if ((a == 1) && (b == 1)) //��һ�����������ڵ㣬ȡ ���ڲ�������
				{
					new_tripoint.push_back(tripoint[m]);
					new_tripoint.push_back(tripoint[m + 1]);
					continue;
				}

			}

			if (interpoint.size() == 2)   //��������
			{
				if ((a == 0) && (b == 0))   //�������������ڵ㣬ȡ������
				{
					new_tripoint.push_back(interpoint[0]);
					new_tripoint.push_back(interpoint[1]);
					continue;
				}
				if ((a == 1) && (b == 0))//����������һ���ڵ� ����������ע��˳�� 
				{
					new_tripoint.push_back(tripoint[m]);
					new_tripoint.push_back(mindistance(tripoint[m], interpoint[0], interpoint[1]));
					continue;
				}
				if ((a == 0) && (b == 1))//����������һ���ڵ� ����������ע��˳��
				{
					new_tripoint.push_back(tripoint[m + 1]);
					new_tripoint.push_back(mindistance(tripoint[m + 1], interpoint[0], interpoint[1]));
					continue;
				}
				if ((a == 1) && (b == 1))//���������������ڵ� ����������ע��˳�� 
				{
					new_tripoint.push_back(tripoint[m]);
					new_tripoint.push_back(mindistance(tripoint[m], interpoint[0], interpoint[1]));
					new_tripoint.push_back(mindistance(tripoint[m + 1], interpoint[0], interpoint[1]));
					new_tripoint.push_back(tripoint[m + 1]);
					continue;
				}

			}

		}

		modelfill.push_back(new_tripoint);
		new_tripoint.clear();
		printf("%.2lf%%\r", i * 100.0 / model.size());
	}


}

void writegcode_ultimaker3(float z) //cura
{
	float t = 0.017776;//��Ӧ0.2�Ĳ��,���0.42
	double E = 0;

	FILE* fp;

	errno_t err;     //�жϴ��ļ����Ƿ���� ���ڷ���1

	err = fopen_s(&fp, "cura.gcode", "a"); //��return 1 , ��ָ������ļ����ļ�����fp1

	fprintf(fp, "T0\n");
	fprintf(fp, "M82\n");
	fprintf(fp, "G92 E0\n");
	fprintf(fp, "M109 S210\n");
	fprintf(fp, "G280 S1\n");
	fprintf(fp, "G0 Z20.001\n");
	fprintf(fp, "G1 F1500 E-6.5\n");
	fprintf(fp, "M107\n");
	fprintf(fp, ";LAYER_COUNT: %d\n", modelfill.size());

	fprintf(fp, ";LAYER:0\n");

	fprintf(fp, "G0 F2000 Z%.1f\n", 0.1);
	for (int i = 0; Bmin_x + 2 * i * 0.42 + 0.42 < Bmax_x; i++)  //��1��2���ʵ����һ��߶���0.27
	{
		fprintf(fp, "G0 F2000 X%f Y%f \n", Bmin_x + 2 * i * 0.42,        Bmin_y);
		fprintf(fp, "G1 F1000 X%f Y%f E%f\n", Bmin_x + 2 * i * 0.42,        Bmax_y, E += (Bmax_y - Bmin_y)*t);
		fprintf(fp, "G0 F2000 X%f Y%f \n",    Bmin_x + 2 * i * 0.42 + 0.42, Bmax_y );
		fprintf(fp, "G1 F1000 X%f Y%f E%f\n", Bmin_x + 2 * i * 0.42 + 0.42, Bmin_y, E += (Bmax_y - Bmin_y)*t);
	}
	fprintf(fp, ";LAYER:1\n");

	fprintf(fp, "G0 F2000 Z%.1f\n", 0.3);
	for (int i = 0; Bmin_y + 2 * i * 0.42 + 0.42 < Bmax_y; i++)  //���ǵڶ��㣬�����Ժ���0.2
	{
		fprintf(fp, "G0 F2000 X%f Y%f \n", Bmin_x , Bmin_y + 2 * i * 0.42);
		fprintf(fp, "G1 F1000 X%f Y%f E%f\n", Bmax_x , Bmin_y + 2 * i * 0.42,        E += (Bmax_y - Bmin_y)*t);
		fprintf(fp, "G0 F2000 X%f Y%f \n",    Bmax_x , Bmin_y + 2 * i * 0.42 + 0.42);
		fprintf(fp, "G1 F1000 X%f Y%f E%f\n", Bmin_x , Bmin_y + 2 * i * 0.42 + 0.42, E += (Bmax_y - Bmin_y)*t);
	}

	
	for (int i = 0; i < model.size()-4; i++) //i���������Ҫ�Ӽ����ʵ�ڵģ����޾ͼ������ٸĸ߶�   
	{
		fprintf(fp, ";LAYER:%d\n", i+2);//���˼���ף�����ͼӼ�
		fprintf(fp, ";TYPE:OUTLINE\n");
		//fprintf(fp, "G0 F2000 Z%.3f\n", 0.500 + i * z);//ȥ��������ʱ�õ���ȷ��Z

		//������
		fprintf(fp, "G0 F2000 X%.3f Y%.3f Z%.1f\n", model[i][0].x, model[i][0].y, 0.3 + i * z);//3��߶ȼ�i*z��i��0��ʼ

		fprintf(fp, "G1 F1200 X%.3f Y%.3f E%.5f\n", model[i][1].x, model[i][1].y, E += distance(model[i][0], model[i][1])*t);
		for (int j = 2; j < model[i].size(); j++)
		{
			fprintf(fp, "G1 X%.3f Y%.3f E%.5f\n", model[i][j].x, model[i][j].y, E += distance(model[i][j - 1], model[i][j])*t);
		}
		fprintf(fp, "G1 X%.3f Y%.3f E%.5f\n", model[i][0].x, model[i][0].y, E += distance(model[i][model[i].size() - 1], model[i][0])*t); //��һ��Ȧ��Ҫ��ԭ��
		fprintf(fp, "G1 X%.3f Y%.3f E%.5f\n", model[i][1].x, model[i][1].y, E += distance(model[i][model[i].size() - 1], model[i][0])*0.01);//�ص�ԭ���ٶ��ȥһ�㣬�����ڲ��̱�

		//�����
		fprintf(fp, ";TYPE:FILL\n");
		for (int k = 1; k < modelfill[i].size(); k += 2)
		{
			fprintf(fp, "G0 F2000 X%.3f Y%.3f\n", modelfill[i][k - 1].x, modelfill[i][k - 1].y);
			fprintf(fp, "G1 F1200 X%.3f Y%.3f E%.5f\n", modelfill[i][k].x, modelfill[i][k].y, E += distance(modelfill[i][k - 1], modelfill[i][k])*t);
		}
	}

	fprintf(fp, ";LAYER:%d\n",model.size() - 2);

	fprintf(fp, "G0 F2000 Z%.1f\n", 0.3 + (model.size() - 3) * z);
	for (int i = 0; Bmin_x + 2 * i * 0.42 + 0.42 < Bmax_x; i++)  //���������ʵ
	{
		fprintf(fp, "G0 F2000 X%f Y%f \n", Bmin_x + 2 * i * 0.42, Bmin_y);
		fprintf(fp, "G1 F1000 X%f Y%f E%f\n", Bmin_x + 2 * i * 0.42, Bmax_y, E += (Bmax_y - Bmin_y)*t);
		fprintf(fp, "G0 F2000 X%f Y%f \n", Bmin_x + 2 * i * 0.42 + 0.42, Bmax_y);
		fprintf(fp, "G1 F1000 X%f Y%f E%f\n", Bmin_x + 2 * i * 0.42 + 0.42, Bmin_y, E += (Bmax_y - Bmin_y)*t);
	}
	fprintf(fp, ";LAYER:%d\n", model.size() - 1);

	fprintf(fp, "G0 F2000 Z%.1f\n", 0.3 + (model.size() - 2) * z);
	for (int i = 0; Bmin_y + 2 * i * 0.42 + 0.42 < Bmax_y; i++)  
	{
		fprintf(fp, "G0 F2000 X%f Y%f Z%.1f\n", Bmin_x, Bmin_y + 2 * i * 0.42, 0.47 + (model.size() - 2) * z);
		fprintf(fp, "G1 F1000 X%f Y%f E%f\n", Bmax_x, Bmin_y + 2 * i * 0.42, E += (Bmax_y - Bmin_y)*t);
		fprintf(fp, "G0 F2000 X%f Y%f \n", Bmax_x, Bmin_y + 2 * i * 0.42 + 0.42);
		fprintf(fp, "G1 F1000 X%f Y%f E%f\n", Bmin_x, Bmin_y + 2 * i * 0.42 + 0.42, E += (Bmax_y - Bmin_y)*t);
	}

	fprintf(fp, "G1 F1500 E%f\n",E-=6);
	fprintf(fp, "M104 S0                     ;extruder heater off\n");
	fprintf(fp, "M107\n");
	fprintf(fp, "G91\n");
	fprintf(fp, "G90                         ;Disable relative movement\n");
	fprintf(fp, "M82                         ;absolute extrusion mode\n");
	fprintf(fp, "M104 S0\n");
	fprintf(fp, "M104 T1 S0\n");
	fprintf(fp, ";End GCode\n");

	fclose(fp);
}

void writegcode(float z) //���Ƶĸ�ʽ������cura
{
	float t = 0.060;//��Ӧ0.2�Ĳ��
	double E = 0;

	FILE* fp;

	errno_t err;     //�жϴ��ļ����Ƿ���� ���ڷ���1

	err = fopen_s(&fp, "cura.gcode", "a"); //��return 1 , ��ָ������ļ����ļ�����fp1

	fprintf(fp, "M104 S190\n");
	fprintf(fp, "M105\n");
	fprintf(fp, "M109 S200\n");
	fprintf(fp, "M82;set extruder to absolute mode\n");
	fprintf(fp, "G28;move X/Y to min endstops\n");
	fprintf(fp, "G1 Z15.0 F6000 ;move the platform down 15mm\n");
	fprintf(fp, "G92 E0                  ;zero the extruded length\n");
	fprintf(fp, "G1 F200 E3             ;extrude 3mm of feed stock\n");
	fprintf(fp, "G92 E0\n");
	fprintf(fp, "G92 E0\n");
	fprintf(fp, "G1 F1500 E-6.5\n");

	for (int i = 0; Bmin_x + i * 4 < Bmax_x; i++)  //����
	{
		fprintf(fp, "G1 F2000 X%f Y%f Z%f F600\n", Bmin_x + i * 4, Bmin_y, 0.3);
		fprintf(fp, "G1 F1000 X%f Y%f E%f\n", Bmin_x + i * 4, Bmax_y, E += (Bmax_y - Bmin_y)*(t + 0.03));
		fprintf(fp, "G1 F300  X%f Y%f E%f\n", Bmin_x + i * 4 + 2, Bmax_y, E += (2)*(t + 0.03));
		fprintf(fp, "G1 F1000 X%f Y%f E%f\n", Bmin_x + i * 4 + 2, Bmin_y, E += (Bmax_y - Bmin_y)*(t + 0.03));
		fprintf(fp, "G1 F300  X%f Y%f E%f\n", Bmin_x + i * 4 + 4, Bmin_y, E += (2)*(t + 0.03));
	}

	fprintf(fp, ";LAYER_COUNT: %d\n", modelfill.size() + 1);
	for (int i = 1; i < model.size(); i++) //ÿһ�� i�������      z��� t������
	{
		fprintf(fp, ";LAYER:%d\n", i);

		fprintf(fp, "G0 F2000 Z%.3f\n", 0.500 + i * z);//ȥ��������ʱ�õ���ȷ��Z

		//������
		fprintf(fp, "G0 F2000 X%.3f Y%.3f Z%.3f\n", model[i][0].x, model[i][0].y, 0.500 + i * z);
		fprintf(fp, ";TYPE:OUTLINE\n");
		fprintf(fp, "G1 F1200 X%.3f Y%.3f E%.5f\n", model[i][1].x, model[i][1].y, E += distance(model[i][0], model[i][1])*t);
		for (int j = 2; j < model[i].size(); j++)
		{
			fprintf(fp, "G1 X%.3f Y%.3f E%.5f\n", model[i][j].x, model[i][j].y, E += distance(model[i][j - 1], model[i][j])*t);
		}
		fprintf(fp, "G1 X%.3f Y%.3f E%.5f\n", model[i][0].x, model[i][0].y, E += distance(model[i][model[i].size() - 1], model[i][0])*t); //��һ��Ȧ��Ҫ��ԭ��
		fprintf(fp, "G1 X%.3f Y%.3f E%.5f\n", model[i][1].x, model[i][1].y, E += distance(model[i][model[i].size() - 1], model[i][0])*0.01);//�ص�ԭ���ٶ��ȥһ�㣬�����ڲ��̱�

		//�����
		fprintf(fp, ";TYPE:FILL\n");
		for (int k = 1; k < modelfill[i].size(); k += 2)
		{
			fprintf(fp, "G0 F2000 X%.3f Y%.3f\n", modelfill[i][k - 1].x, modelfill[i][k - 1].y);
			fprintf(fp, "G1 F1200 X%.3f Y%.3f E%.5f\n", modelfill[i][k].x, modelfill[i][k].y, E += distance(modelfill[i][k - 1], modelfill[i][k])*t);
		}
	}

	fprintf(fp, "M107\n");
	fprintf(fp, "M104 S0                     ;extruder heater off\n");
	fprintf(fp, "G92 E1                      ;relative positioning\n");
	fprintf(fp, "G1 E-1 F300                 ;retract the filament a bit before lifting the nozzle, to release some of the pressure\n");
	fprintf(fp, "G1 X120 Y120 Z120 F9000     ;move to max so the object can take out\n");
	fprintf(fp, "M84                         ;steppers off\n");
	fprintf(fp, "M82                         ;absolute extrusion mode\n");
	fprintf(fp, "M104 S0\n");
	fprintf(fp, ";End GCode\n");

	fclose(fp);
}

void writegcode_makerbot(float z)//������cura�ĸ�ʽȻ����x3gת����
{
	int layers = model.size();//��ӡ����
	float layer_height = z;//���,�Ĵ˲���Ҫ�Ĳ�����t
	float t = 0.0349595312;
	//if(layer_height == 0.2)      t = 0.0349595312;//�����ʲ��Ϊ0.2
	//else if(layer_height == 0.3) t = 0.0522797843;//�����ʲ��Ϊ0.3
	int F = 1800;//��ӡ�ٶ� mm/min

	double A = 0.0;//������

	FILE* fp;

	errno_t err;     //�жϴ��ļ����Ƿ���� ���ڷ���1

	err = fopen_s(&fp, "makerbot.gcode", "a"); //��return 1 , ��ָ������ļ����ļ�����fp1

	//��һ����һ��
	fprintf(fp, "M73 P%d\n", 0);//��ӡ���Ȱٷֱ�
	fprintf(fp, "G1 X%f Y%f Z%f F3000;travel move\n", Bmin_x - 0.01, Bmin_y - 0.01, 0.3);
	//for (int i = 0; Bmin_x+i*4<Bmax_x; i++)  //����
	//{
	//	fprintf(fp, "G1 X%f Y%f Z%f F600\n",     Bmin_x + i * 4,   Bmin_y, 0.3);
	//	fprintf(fp, "G1 X%f Y%f Z%f F600 A%f\n", Bmin_x + i * 4,   Bmax_y, 0.3, A += (Bmax_y- Bmin_y)*(t+0.03));
	//	fprintf(fp, "G1 X%f Y%f Z%f F600 A%f\n", Bmin_x + i * 4+2, Bmax_y, 0.3, A += (2)*(t+0.03));
	//	fprintf(fp, "G1 X%f Y%f Z%f F600 A%f\n", Bmin_x + i * 4+2, Bmin_y, 0.3, A += (Bmax_y - Bmin_y)*(t+0.03));
	//	fprintf(fp, "G1 X%f Y%f Z%f F600 A%f\n", Bmin_x + i * 4+4, Bmin_y, 0.3, A += (2)*(t+0.03));
	//}
	fprintf(fp, ";LAYER:0\n");
	for (int i = 0; Bmin_x + 2 * i * 0.42 + 0.42 < Bmax_x; i++)  //��1��2���ʵ����һ��߶���0.27
	{
		fprintf(fp, "G0 F600 X%f Y%f Z%f\n", Bmin_x + 2 * i * 0.42, Bmin_y, 0.27);
		fprintf(fp, "G1 F600 X%f Y%f Z%f A%f\n", Bmin_x + 2 * i * 0.42, Bmax_y, 0.27,A += (Bmax_y - Bmin_y)*t);
		fprintf(fp, "G0 F600 X%f Y%f Z%f\n", Bmin_x + 2 * i * 0.42 + 0.42, Bmax_y, 0.27);
		fprintf(fp, "G1 F600 X%f Y%f Z%f A%f\n", Bmin_x + 2 * i * 0.42 + 0.42, Bmin_y, 0.27,A += (Bmax_y - Bmin_y)*t);
	}
	fprintf(fp, ";LAYER:1\n");
	for (int i = 0; Bmin_y + 2 * i * 0.42 + 0.42 < Bmax_y; i++)  //���ǵڶ��㣬�����Ժ���0.2
	{
		fprintf(fp, "G0 F600 X%f Y%f Z%f\n", Bmin_x, Bmin_y + 2 * i * 0.42, 0.47);
		fprintf(fp, "G1 F600 X%f Y%f Z%f A%f\n", Bmax_x, Bmin_y + 2 * i * 0.42, 0.47,A += (Bmax_y - Bmin_y)*t);
		fprintf(fp, "G0 F600 X%f Y%f Z%f\n", Bmax_x, Bmin_y + 2 * i * 0.42 + 0.42, 0.47);
		fprintf(fp, "G1 F600 X%f Y%f Z%f A%f\n", Bmin_x, Bmin_y + 2 * i * 0.42 + 0.42, 0.47,A += (Bmax_y - Bmin_y)*t);
	}

	for (int i = 0; i < model.size()-4; i++) //��ӡ���������,��ӡ����ʵ�ľͼ������ٸĸ߶�
	{
		fprintf(fp, ";LAYER:%d\n",i+2);
		fprintf(fp, "M73 P%d\n", int(i / layers) * 100);//��ӡ���Ȱٷֱ�
		fprintf(fp, "G0 X%f Y%f Z%f F2000;travel move\n", modelfill[i][0].x-0.01, modelfill[i][0].y-0.01, 0.67 + i * z);

		fprintf(fp, "G0 F2000 X%.3f Y%.3f Z%.3f\n", model[i][0].x, model[i][0].y, 0.67 + i * z);
		fprintf(fp, ";TYPE:OUTLINE\n");
		fprintf(fp, "G1 F1200 X%.3f Y%.3f Z%f A%.5f\n", model[i][1].x, model[i][1].y, 0.67 + i * z, A += distance(model[i][0], model[i][1])*t);
		for (int j = 2; j < model[i].size(); j++)
		{
			fprintf(fp, "G1 X%.3f Y%.3f Z%f A%.5f\n", model[i][j].x, model[i][j].y, 0.67 + i * z, A += distance(model[i][j - 1], model[i][j])*t);
		}
		fprintf(fp, "G1 X%.3f Y%.3f Z%f A%.5f\n", model[i][0].x, model[i][0].y, 0.67 + i * z, A += distance(model[i][model[i].size() - 1], model[i][0])*t); //��һ��Ȧ��Ҫ��ԭ��
		fprintf(fp, "G1 X%.3f Y%.3f Z%f A%.5f\n", model[i][1].x, model[i][1].y, 0.67 + i * z, A += distance(model[i][model[i].size() - 1], model[i][0])*0.01);//�ص�ԭ���ٶ��ȥһ�㣬�����ڲ��̱�

		fprintf(fp, ";TYPE:FILL\n");
		for (int k = 1; k < modelfill[i].size(); k += 2)
		{
			fprintf(fp, "G1 X%f Y%f Z%f F2000\n", modelfill[i][k - 1].x, modelfill[i][k - 1].y, 0.67 + i * z);
			fprintf(fp, "G1 X%f Y%f Z%f F1800 A%f\n", modelfill[i][k].x, modelfill[i][k].y, 0.67 + i * z, A += distance(modelfill[i][k - 1], modelfill[i][k])*t);
		}
	}

	//for (int i = 0; Bmin_x + 2 * i * 0.42 + 0.42 < Bmax_x; i++)  //���������ʵ
	//{
	//	fprintf(fp, "G0 X%f Y%f Z%f F2000\n", Bmin_x + 2 * i * 0.42, Bmin_y, 0.47 + (model.size() - 3) * z);
	//	fprintf(fp, "G1 X%f Y%f Z%f A%f F1000\n", Bmin_x + 2 * i * 0.42, Bmax_y, 0.47 + (model.size() - 3) * z, A += (Bmax_y - Bmin_y)*t);
	//	fprintf(fp, "G0 X%f Y%f Z%f F2000\n", Bmin_x + 2 * i * 0.42 + 0.42, Bmax_y,0.47 + (model.size() - 3) * z);
	//	fprintf(fp, "G1 X%f Y%f Z%f A%f F1000\n", Bmin_x + 2 * i * 0.42 + 0.42, Bmin_y, 0.47 + (model.size() - 3) * z, A += (Bmax_y - Bmin_y)*t);
	//}
	//fprintf(fp, ";LAYER:%d\n", model.size() - 1);
	//for (int i = 0; Bmin_y + 2 * i * 0.42 + 0.42 < Bmax_y; i++)
	//{
	//	fprintf(fp, "G0 X%f Y%f Z%f F2000\n", Bmin_x, Bmin_y + 2 * i * 0.42, 0.47 + (model.size() - 2) * z);
	//	fprintf(fp, "G1 X%f Y%f Z%f A%f F1000 \n", Bmax_x, Bmin_y + 2 * i * 0.42, 0.47 + (model.size() - 2) * z, A += (Bmax_y - Bmin_y)*t);
	//	fprintf(fp, "G0 X%f Y%f Z%f F2000\n", Bmax_x, Bmin_y + 2 * i * 0.42 + 0.42, 0.47 + (model.size() - 2) * z);
	//	fprintf(fp, "G1 X%f Y%f Z%f A%f F1000 \n", Bmin_x, Bmin_y + 2 * i * 0.42 + 0.42, 0.47 + (model.size() - 2) * z, A += (Bmax_y - Bmin_y)*t);
	//}
	fclose(fp);
}

void InitGL(GLvoid)
{
	glShadeModel(GL_SMOOTH);
	glClearColor(1.0, 1.0, 1.0, 1.0);
	glClearDepth(1.0f);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);
	glEnable(GL_COLOR_MATERIAL);
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
}

void lines() //���ߵĺ���
{
	//glColor3f(1.0, 0.0, 0.0);
	//glBegin(GL_LINES);
	//for (int j = 0; j < 50; j++)
	//{
	//	for (int i = 0; i < tripoint.size(); i += 2)
	//	{
	//		glVertex3f(tripoint[i].x, tripoint[i].y, 0);  glVertex3f(tripoint[i + 1].x, tripoint[i + 1].y, 0);
	//	}
	//}
	//glEnd();
	glBegin(GL_LINES);//������
	glColor3f(1.0, 0.0, 0.0);  glVertex3f(0, 0, 0); glVertex3f(100, 0, 0);
	glColor3f(0.0, 1.0, 0.0); glVertex3f(0, 0, 0); glVertex3f(0, 100, 0);
	glColor3f(0.0, 0.0, 1.0); glVertex3f(0, 0, 0); glVertex3f(0, 0, 100);
	glEnd();

	glColor3f(0.945, 0.945, 0.945);
	glBegin(GL_LINES);
	for (int i = 1; i < model.size(); i++)
	{
		for (int j = 1; j < model[i].size(); j++)
		{
			glVertex3f(model[i][j - 1].x, model[i][j - 1].y, i*0.2); glVertex3f(model[i][j].x, model[i][j].y, i*0.2);
		}
		glVertex3f(model[i][0].x, model[i][0].y, i*0.2); glVertex3f(model[i][model[i].size() - 1].x, model[i][model[i].size() - 1].y, i*0.2);
	}
	glEnd();

	glColor3f(0.615, 0.615, 0.615);
	glBegin(GL_LINES);
	for (int i = 1; i < modelfill.size(); i += 19)
	{
		for (int j = 1; j < modelfill[i].size(); j += 2)
		{
			glVertex3f(modelfill[i][j - 1].x, modelfill[i][j - 1].y, i*0.2); glVertex3f(modelfill[i][j].x, modelfill[i][j].y, i*0.2);
		}
	}
	glEnd();

}

void display(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();
	glPushMatrix();

	gluLookAt(30 + xs, -30 + ys, 100 + zs, 30 + sx, 30 + sy, 30 + sz, 0, 0, 1);
	//glTranslatef(0.0f, 0.0f, -1.0f);	//ƽ��
	//glRotatef(rquad, 1.0f, 0.0f, 0.0f);	//��תһ���Ƕ�
	lines();

	//glColor3f(1.0, 0.0, 0.0);
	//glPointSize(5);
	//glColor3f(1.0,0.0,0.0);
	//glBegin(GL_POINTS);
	////glVertex3f(0, 0, 0);
	////glVertex3f(-50, 0, 1);
	//glVertex3f(60,60,60);
	//glEnd();
	glPopMatrix();

	//rquad -= 0.15f;	//�޸����������ת�Ƕ�
	glutSwapBuffers();
}

void reshape(int width, int height)
{
	if (height == 0)
		height = 1;
	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	//gluPerspective(90.0f, (GLfloat)width / (GLfloat)height, 0.0f, 100.0f);
	gluPerspective(90.0f, (GLfloat)width / (GLfloat)height, 0.0f, 100.0f);
	//if (width <= height)
	//	glOrtho(-2.0, 2.0, -2.0*(GLfloat)height / (GLfloat)width, 2.0*(GLfloat)height / (GLfloat)width, 1.0, 20.0);
	//else
	//	glOrtho(-2.0*(GLfloat)width / (GLfloat)height, 2.0*(GLfloat)width / (GLfloat)height, -2.0, 2.0, 1.0, 20.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void myMouse(int button, int state, int x, int y)
{
	//if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) //����������
	//{
	//	mousetate = 1;
	//	Oldx = x;
	//	Oldy = y;
	//}
	//if (button == GLUT_LEFT_BUTTON && state == GLUT_UP) //�������ɿ�
	//{
	//	mousetate = 0;
	//}
	//�����¼���С
	if (state == GLUT_UP && button == 3) {

		scale -= 0.1f;
	}//�Ŵ�
	if (state == GLUT_UP && button == 4) {

		scale += 0.1f;
	}
	glutPostRedisplay();//����ػ�ý���
}

void onMouseMove(int x, int y)//��������
{

}

void keyboard(unsigned char key, int x, int y)
{
	switch (key)
	{
	case 'w':
		sy += 5;
		glutPostRedisplay();
		break;
	case 's':
		sy -= 5;
		glutPostRedisplay();
		break;
	case 'a':
		sx -= 5;
		glutPostRedisplay();
		break;
	case 'd':
		sx += 5;
		glutPostRedisplay();
		break;
	case 'q':
		sz += 5;
		glutPostRedisplay();
		break;
	case 'e':
		sz -= 5;
		glutPostRedisplay();
		break;
	default:
		break;
	}
}

void SpecialKey(GLint key, GLint x, GLint y)
{
	if (key == GLUT_KEY_UP)
	{
		ys += 10;
	}
	if (key == GLUT_KEY_LEFT)
	{
		xs -= 10;
	}
	if (key == GLUT_KEY_DOWN)
	{
		ys -= 10;
	}
	if (key == GLUT_KEY_RIGHT)
	{
		xs += 10;
	}
	glutPostRedisplay();
}

void main(int argc, char** argv)
{
	float z = 0.2;
	//************************************************************************��Ƭ
	printf("slicing...\n");
	readfile(file_1);
	BoundingBox();
	findIntersect(z);  //ǰ����������Ƭ������model
	printf("slice complete,layer count: %d\n", countt);

	//************************************************************************ȡ��
	printf("path planning...\n");
	intersection();
	printf("path planning complete\n"); //modelfill���

	//************************************************************************��дGcode
	printf("Gcode writing...\n");
	//writegcode(z);
	//writegcode_makerbot(z);
	writegcode_ultimaker3(z);
	//writpointcloud();
	printf("Gcode writing complete,file save as \"test gcode.gcode\"\n");

	//************************************************************************Ԥ��
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
	glutInitWindowSize(600, 600);
	glutCreateWindow("Hello Cube");
	InitGL();
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutIdleFunc(display);
	glutMouseFunc(myMouse);
	glutKeyboardFunc(keyboard);
	glutSpecialFunc(&SpecialKey);
	glutMainLoop();

}


