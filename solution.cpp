//Do not edit the code below (unless you know what you are doing)

#include"solution.h"

int solution::f_calls = 0;
int solution::g_calls = 0;
int solution::H_calls = 0;

solution::solution(double L)
{
	x = matrix(L);
	g = NAN;
	H = NAN;
	y = NAN;
}

solution::solution(const matrix& A)
{
	x = A;
	g = NAN;
	H = NAN;
	y = NAN;
}

solution::solution(double* A, int n)
{
	x = matrix(A, n);
	g = NAN;
	H = NAN;
	y = NAN;
}

void solution::clear_calls()
{
	f_calls = 0;
	g_calls = 0;
	H_calls = 0;
}

ostream& operator<<(ostream& S, const solution& A)
{
	S << "x = " << A.x << endl;
	S << "y = " << A.y << endl;
	S << "f_calls = " << solution::f_calls << endl;
	S << "g_calls = " << solution::g_calls << endl;
	S << "H_calls = " << solution::H_calls << endl;
	return S;
}

//You can edit the following code

void solution::fit_fun(matrix O)
{
#if LAB_NO == 2 && LAB_PART==1
	y = -cos(0.1 * x(0)) * exp(-pow(0.1 * x(0) - 2 * 3.14, 2)) + 0.002 * pow(0.1 * x(0), 2);
#elif LAB_NO == 2 && LAB_PART==2
	double target_temp = 50;
	matrix Y0(new double[3]{ 5.0,1.0,10.0 }, 3);
	matrix* Y = solve_ode(0, 1, 1000, Y0, x);
	int* w = get_size(Y[1]);
	double max = Y[1](0, 2);
	for (int i = 1; i < w[0]; i++)
	{
		if (max < Y[1](i, 2))
		{
			max = Y[1](i, 2);
		}
	}

	y = abs(max - target_temp);
#elif LAB_NO== 3 && LAB_PART==1
	y = x(0) * x(0) + x(1) * x(1) - cos(2.5 * 3.14 * x(0)) - cos(2.5 * 3.14 * x(1)) + 2;
#elif LAB_NO== 3 && LAB_PART==2
	double a_ref = 3.14, o_ref = 0;
	matrix Y0(2, 1);
	double dt = 0.1;
	matrix* Y = solve_ode(0, dt, 100, Y0, x);
	int* n = get_size(Y[1]);
	y(0) = 0;
	for (int i = 0; i < n[0]; i++)
	{
		y = y + 10 * pow(a_ref - Y[1](i, 0), 2) +
			pow(o_ref - Y[1](i, 1), 2) +
			pow(x(0) * (a_ref - Y[1](i, 0)) + x(1) * (o_ref - Y[1](i, 1)), 2);
		cout << Y[1](i, 0) << " " << Y[1](i, 1) << " " << i * dt << endl;
	}
	y = y * dt;
#elif LAB_NO == 4
#if LAB_PART == 1
	double arg = 3.14 * sqrt(pow(x(0) / 3.14, 2));
	y = sin(arg) / arg;
	if (-x(0) + 1 > 0)
		y = y + O(0) * pow(-x(0) + 1, 2);
	if (-x(1) + 1 > 0)
		y = y + O(0) * pow(-x(1) + 1, 2);
	if (sqrt(pow(x(0), 2) + pow(x(1), 2)) - O(1) > 0)
		y = y + O(0) * pow(sqrt(pow(x(0), 2) + pow(x(1), 2)) - O(1), 2);

#elif LAB_PART == 2
	double arg = 3.14 * sqrt(pow(x(0) / 3.14, 2));
	y = sin(arg) / arg;
	if (-x(0) + 1 > 0)
	{
		y = 1e10;
		return;
	}
	else
		y = y - O(0) / (-x(0) + 1);
	if (-x(1) + 1 > 0)
	{
		y = 1e10;
		return;
	}
	else
		y = y - O(0) / (-x(1) + 1);
	if (sqrt(pow(x(0), 2) + pow(x(1), 2)) - O(1) > 0)
	{
		y = 1e10;
		return;
	}
	else
		y = y - O(0) / (sqrt(pow(x(0), 2) + pow(x(1), 2)) - O(1));


#elif LAB_PART == 3
#endif
#elif LAB_NO == 5
	int* n = get_size(O);
	if (n[1] == 1)
	{
		y = pow(x(0) + 2 * x(1) - 7, 2) + pow(2 * x(0) + x(1) - 5, 2);
	}
	else
	{
		solution temp;
		temp.x = O[0] + x * O[1];
		temp.fit_fun();
		y = temp.y;
		--f_calls;
	}
#elif LAB_NO==6

#if LAB_PART == 1
	int* n = get_size(O);
	if (n[1] == 1)
	{
		int a = 1;
		y = matrix(2, 1);
		y(0) = a * (pow(x(0) - 5, 2) + pow(x(1) - 5, 2));
		y(1) = 1.0 / a * (pow(x(0) + 5, 2) + pow(x(1) + 5, 2));
	}
	else
	{
		solution temp;
		temp.x = O[0] + x * O[1];
		temp.fit_fun();
		y = O(0, 2) * temp.y(0) + (1 - O(0, 2)) * temp.y(1);
		--f_calls;
	}
#elif LAB_PART == 2
	int* n = get_size(O);
	if (n[1] == 1)
	{
		double ro = 7800, P = 1e3, E = 207e9;
		y = matrix(3, 1);
		y(0) = ro * x(0) * 3.14 * pow(x(1), 2) / 4;
		y(1) = 64 * P * pow(x(0), 3) / (3 * E * 3.14 * pow(x(1), 4));
		y(2) = 32 * P * x(0) / (3.14 * pow(x(1), 3));
	}
	else
	{
		solution temp;
		temp.x = O[0] + x * O[1];
		temp.fit_fun();
		y = O(0, 2) * temp.y(0) + (1 - O(0, 2)) * temp.y(1);
		if (temp.y(1) > 0.005)
		{
			y = y + 1e6 * pow(temp.y(1) - 0.005, 2);
		}
		if (temp.y(2) > 300e6)
		{
			y = y + 1e6 * pow(temp.y(2) - 300e6, 2);
		}
		--f_calls;
	}
#endif
#endif
	++f_calls;
}

void solution::grad(matrix O)
{
	g = matrix(2, 1);
	g(0) = 10 * x(0) + 8 * x(1) - 34;
	g(1) = 8 * x(0) + 10 * x(1) - 38;
	++g_calls;
}

void solution::hess(matrix O)
{
	H = matrix(2, 2);
	H(0, 0) = H(1, 1) = 10;
	H(0, 1) = H(1, 0) = 8;
	++H_calls;
}
