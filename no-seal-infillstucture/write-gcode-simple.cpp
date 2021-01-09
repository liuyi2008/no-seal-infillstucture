#include <math.h>
#include <vector>

std::vector<point>tripoint;  //һ��ı仯�� �߶����� �ɶ�

struct point {

	double x, y;

};

void main()
{
	float big_a = 10.0; 
	float small_a = 5.0;
	float a = 0;

	float s = 0.2;//����,Ҫ��������ߵĿ�ȣ�0<s<width/2��

	float distance = (big_a - small_a) / (2 * sqrt(3));

	int ii;
	for (int i = 0; i < 250; i++)
	{
		ii = i;//����iȥ�����㣬��Ӱ��iʼ�մ�����ʵ����

		while (ii - 1 > 4 * distance / s)
		{
			ii = ii - int(4 * distance / s);
		}
		int condition;
		if      (ii - 1 < 1 * distance / s) condition = 1;
		
		else if (ii - 1 < 2 * distance / s) condition = 2;

		else if (ii - 1 < 3 * distance / s) condition = 3;

		else if (ii - 1 < 4 * distance / s) condition = 4;

		switch (condition)
		{
		case 1:
		{

			point A = { 0,0 };
			point B = { big_a,0 };
			point D = { s*sqrt(3)*ii,s*ii };
			point E = { big_a - s * sqrt(3)*ii,s*ii };

			tripoint.push_back(A); tripoint.push_back(D); tripoint.push_back(E); tripoint.push_back(B);

		}

		case 2:
		{
			ii = ii - int(sqrt(3)*0.333*a / s);

			point A = { 0,0 };
			point B = { big_a,0 };
			point D = { (distance-s*ii)*sqrt(3),distance - s*ii };
			point E = { big_a - (distance-s*ii)*sqrt(3),distance - s * ii };

			tripoint.push_back(A); tripoint.push_back(D); tripoint.push_back(E); tripoint.push_back(B);

		}

		case 3://��1��ȣ������겻��������ȡ��
		{
			ii = ii - int(2 * sqrt(3)*0.333*a / s);

			point A = { 0,0 };
			point B = { big_a,0 };
			point D = { s*sqrt(3)*ii, -s*ii };
			point E = { big_a - s * sqrt(3)*ii, -s*ii };

			tripoint.push_back(A); tripoint.push_back(D); tripoint.push_back(E); tripoint.push_back(B);

		}

		case 4:
		{
			ii = ii - int(3 * sqrt(3)*0.333*a / s);

			point A = { 0,0 };
			point B = { big_a,0 };
			point D = { (distance - s * ii)*sqrt(3),-distance + s * ii };
			point E = { big_a - (distance - s * ii)*sqrt(3),-distance + s * ii };

			tripoint.push_back(A); tripoint.push_back(D); tripoint.push_back(E); tripoint.push_back(B);

		}

		}

	}
	
}