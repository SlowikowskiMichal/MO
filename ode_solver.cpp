//Do not edit the code below (unless you know what you are doing)

#include"ode_solver.h"

matrix *solve_ode(double t0, double dt, double tend, const matrix &Y0, matrix P)
{
	int N = static_cast<int>(floor((tend - t0) / dt) + 1);
	if (N < 2)
		throw "The time interval is defined incorrectly";
	int* s = get_size(Y0);
	if (s[1] != 1)
		throw "Initial condition must be a vector";
	int n = s[0];
	delete[]s;
	matrix* S = new matrix[2]{ matrix(N), matrix(n,N) };
	S[0](0) = t0;
	for (int i = 0; i < n; ++i)
		S[1](i, 0) = Y0(i);
	matrix k1(n), k2(n), k3(n), k4(n);
	for (int i = 1; i < N; ++i)
	{
		S[0](i) = S[0](i - 1) + dt;
		k1 = dt * diff(S[0](i - 1), S[1][i - 1], P);
		k2 = dt * diff(S[0](i - 1) + 0.5 * dt, S[1][i - 1] + 0.5 * k1, P);
		k3 = dt * diff(S[0](i - 1) + 0.5 * dt, S[1][i - 1] + 0.5 * k2, P);
		k4 = dt * diff(S[0](i - 1) + dt, S[1][i - 1] + k3, P);
		for (int j = 0; j < n; ++j)
			S[1](j, i) = S[1](j, i - 1) + (k1(j) + 2 * k2(j) + 2 * k3(j) + k4(j)) / 6;
	}
	S[1] = trans(S[1]);
	return S;
}

//You can edit the following code

matrix diff(double t, const matrix &Y, matrix P)
{
#if LAB_NO==1
	matrix dY(3);
	double a = 0.98, b = 0.63, g = 9.81, PA = 1.0, PB = 1.0,
		TA = 90.0, DA = 0.005, DB = 0.00365665, Fin = 0.01, Tin = 10.0;

	double FAout = Y(0) > 0.0 ?  a*b*DA*sqrt(2 * g*Y(0) / PA) : 0.0;
	double FBout = Y(1) > 0.0 ? a*b*DB*sqrt(2 * g*Y(1) / PB) : 0.0;

	dY(0) = -FAout;
	dY(1) = FAout + Fin - FBout;
	dY(2) = Fin / Y(1)*(Tin - Y(2)) + FAout / Y(1)*(TA - Y(2));
#elif LAB_NO==2 && LAB_PART==2 && TASK == 4 
	matrix dY(3);
	double a = 0.98, b = 0.63, g = 9.81, PA = 1.0, PB = 1.0, // DA nie jest sta��, teraz pr�bujemy znale�� �rednic� otworu zbiornika A
		TA = 90.0, DB = 0.00365665, Fin = 0.01, Tin = 10.0;

	double FAout = Y(0) > 0.0 ? a*b*P(0)*sqrt(2 * g*Y(0) / PA) : 0.0; // przekazujemy x za pomoc� parametru P
	double FBout = Y(1) > 0.0 ? a*b*DB*sqrt(2 * g*Y(1) / PB) : 0.0;

	dY(0) = -FAout;
	dY(1) = FAout + Fin - FBout;
	dY(2) = Fin / Y(1)*(Tin - Y(2)) + FAout / Y(1)*(TA - Y(2));
#elif LAB_NO==2 && LAB_PART==2 && TASK == 5 //policzenie warto�ci z LAB_NO = 1 dla nowego DA
	matrix dY(3);
	double a = 0.98, b = 0.63, g = 9.81, PA = 1.0, PB = 1.0,
		TA = 90.0, DA = 0.002578566, DB = 0.00365665, Fin = 0.01, Tin = 10.0;

	double FAout = Y(0) > 0.0 ? a * b*DA*sqrt(2 * g*Y(0) / PA) : 0.0;
	double FBout = Y(1) > 0.0 ? a * b*DB*sqrt(2 * g*Y(1) / PB) : 0.0;

	dY(0) = -FAout;
	dY(1) = FAout + Fin - FBout;
	dY(2) = Fin / Y(1)*(Tin - Y(2)) + FAout / Y(1)*(TA - Y(2));
#elif LAB_NO==3&&LAB_PART==2
	double mr = 1, l = 0.5, mc = 10, b = 0.5, a_ref = 3.14, o_ref = 0;
	double I = mr * l * l / 3 + mc * l * l;
	matrix dY(Y);
	dY(0) = Y(1);
	dY(1) = (P(0) * (a_ref - Y(0)) + (P(1) * (o_ref - Y(1))) - b * Y(1)) / I;
#elif LAB_NO==4 && LAB_PART==3
	double C = 0.47, r = 0.12, m = 0.6, ro = 1.2, g = 9.81;
	double S = 3.14 * r * r,
		Dx = 0.5 * C * ro * S * Y(1) * abs(Y(1)),
		Dy = 0.5 * C * ro * S * Y(3) * abs(Y(3)),
		FMx = 3.14 * ro * Y(3) * P(0) * pow(r, 3),
		FMy = 3.14 * ro * Y(1) * P(0) * pow(r, 3);
	matrix dY(Y);
	dY(0) = Y(1);
	dY(1) = (-Dx - FMx) / m;
	dY(2) = Y(3);
	dY(3) = (-Dy - FMy - m * g) / m;
	return dY
#else 
	matrix dY(Y);
	
#endif
	return dY;
}