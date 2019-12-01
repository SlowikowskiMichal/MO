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

solution::solution(const matrix &A)
{
	x = A;
	g = NAN;
	H = NAN;
	y = NAN;
}

solution::solution(double *A, int n)
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

ostream &operator<<(ostream &S, const solution &A)
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
	y = -cos(0.1*x(0))*exp(-pow(0.1*x(0) - 2 * 3.14, 2)) + 0.002*pow(0.1*x(0), 2);
	#elif LAB_NO == 2 && LAB_PART==2
	double target_temp = 50;
	matrix Y0(new double[3]{ 5.0,1.0,10.0 }, 3);
	matrix *Y = solve_ode(0, 1, 1000, Y0, x);
	int *w = get_size(Y[1]);
	double max = Y[1](0, 2);
	for (int i = 1; i < w[0]; i++)
	{
		if (max < Y[1](i, 2))
		{
			max = Y[1](i, 2);
		}
	}

	y = abs(max - target_temp);

	#endif
	++f_calls;
}

void solution::grad(matrix O)
{
	++g_calls;
}

void solution::hess(matrix O)
{
	++H_calls;
}
