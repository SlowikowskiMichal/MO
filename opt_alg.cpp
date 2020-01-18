#include"opt_alg.h"
#if LAB_NO>1
double* expansion(double x0, double d, double alfa, int Nmax, matrix O)
{
	double* p = new double[2];
	solution X0(x0), X1(x0 + d);
	X0.fit_fun();
	X1.fit_fun();
	if (X0.y == X1.y)
	{
		p[0] = X0.x(0);
		p[1] = X1.x(0);
		return p;
	}
	if (X0.y < X1.y)
	{
		d *= -1;
		X1.x = X0.x + d;
		X1.fit_fun();
		if (X0.y <= X1.y)
		{
			p[0] = X1.x(0);
			p[1] = X0.x(0) - d;
			return p;
		}
	}
	solution X2;
	int i = 1;
	while (true)
	{
		X2.x = x0 + pow(alfa, i) * d;
		X2.fit_fun();
		if (X2.y >= X1.y || solution::f_calls > Nmax)
			break;
		X0 = X1;
		X1 = X2;
		++i;
	}
	d > 0 ? p[0] = X0.x(0), p[1] = X2.x(0) : (p[0] = X2.x(0), p[1] = X0.x(0));
	return p;
}

solution fib(double a, double b, double epsilon, matrix O)
{
#if LAB_NO == 2 && LAB_PART == 1 && TASK == 2
	fstream file("wynik1.txt", ios::app);
#endif
	int n = static_cast<int>(ceil(log2(sqrt(5) * (b - a) / epsilon) / log2((1 + sqrt(5)) / 2)));
	int* F = new int[n] {1, 1};
	for (int i = 2; i < n; ++i)
		F[i] = F[i - 2] + F[i - 1];
	solution A(a), B(b), C, D;
	C.x = B.x - (1.0 * F[n - 2] / F[n - 1]) * (B.x - A.x);
	D.x = A.x + B.x - C.x;
	C.fit_fun();
	D.fit_fun();
	for (int i = 0; i <= n - 4; ++i)
	{
#if LAB_NO == 2 && LAB_PART == 1 && TASK == 2
		file << B.x - A.x << endl;
#endif
		if (C.y < D.y)
			B = D;
		else
			A = C;
		C.x = B.x - 1.0 * F[n - i - 2] / F[n - i - 1] * (B.x - A.x);
		D.x = A.x + B.x - C.x;
		C.fit_fun();
		D.fit_fun();

	}
#if LAB_NO == 2 && LAB_PART == 1 && TASK == 2
	file.close();
#endif
	return C;
}

solution lag(double a, double b, double epsilon, double gamma, int Nmax, matrix O)
{ //w ksi¹¿ce 2 razy F(c) powinno byæ raz F(D) i m < 0, tutaj musi byæ m > 0 przez odwrócony znak
#if LAB_NO == 2 && LAB_PART == 1 && TASK == 2
	fstream file("wynik2.txt", ios::app);
#endif
	solution A(a), B(b), C, D;
	C.x = (a + b) / 2.0;
	A.fit_fun();
	B.fit_fun();
	C.fit_fun();
	double l, m;
	while (true)
	{
#if LAB_NO == 2 && LAB_PART == 1 && TASK == 2
		file << B.x - A.x << endl;
#endif
		l = A.y(0) * (pow(B.x(0), 2) - pow(C.x(0), 2)) + B.y(0) * (pow(C.x(0), 2) - pow(A.x(0), 2)) + C.y(0) * (pow(A.x(0), 2) - pow(B.x(0), 2));
		m = A.y(0) * (B.x(0) - C.x(0)) + B.y(0) * (C.x(0) - A.x(0)) + C.y(0) * (A.x(0) - B.x(0));
		if (m <= 0)
		{
			C.x = NAN;
			C.y = NAN;
			return C;
		}
		D.x = 0.5 * l / m;
		D.fit_fun();
		if (A.x <= D.x && D.x < C.x) //A D C B
		{
			if (D.y < C.y) {
				B = C;
				C = D;
			}
			else
			{
				A = D;
			}
		}
		else {
			if (C.x <= D.x && D.x <= B.x) // A C D B
			{
				if (D.y < C.y) {
					A = C;
					C = D;
				}
				else
				{
					B = D;
				}
			}
			else
			{
				C.x = NAN;
				C.y = NAN;
				return C;
			}
		}
		if (B.x - A.x < epsilon ||
			abs(C.x(0) - D.x(0)) <= gamma
			|| solution::f_calls > Nmax) {
			//return C;
#if LAB_NO == 2 && LAB_PART == 1 && TASK == 2
			file.close();
#endif
			return D;
		}
	}
}
#endif
#if LAB_NO>2
solution HJ(matrix x0, double s, double alfa, double epsilon, int Nmax, matrix O)
{
	solution XB, XB_old, X;
	XB.x = x0;
	XB.fit_fun();
	cout << XB.x(0) << " " << XB.x(1) << endl;
	while (true)
	{
		X = HJ_trial(XB, s);
		if (X.y < XB.y)
		{
			while (true)
			{
				XB_old = XB;
				XB = X;
				X.x = 2.0 * XB.x - XB_old.x;
				X.fit_fun();
				X = HJ_trial(X, s);

				cout << XB.x(0) << " " << XB.x(1) << endl;
				if (X.y >= XB.y)
				{

					break;
				}
				if (solution::f_calls > Nmax)
				{
					cout << XB.x(0) << " " << XB.x(1) << endl;
					return XB;
				}
			}
		}
		else
			s *= alfa;
		cout << XB.x(0) << " " << XB.x(1) << endl;
		if (s < epsilon || solution::f_calls>Nmax)
			return XB;
	}
}

solution HJ_trial(solution XB, double s, matrix O)
{
	int* n = get_size(XB.x);
	matrix D = unit_mat(n[0]);
	solution X;
	for (int i = 0; i < n[0]; ++i)
	{
		X.x = XB.x + s * D[i];
		X.fit_fun();
		if (X.y < XB.y)
			XB = X;
		else
		{
			X.x = XB.x - s * D[i];
			X.fit_fun();
			if (X.y < XB.y)
				XB = X;
		}
	}
	return XB;
}

solution Rosen(matrix x0, matrix s0, double alfa, double beta, double epsilon, int Nmax, matrix O)
{
	int* n = get_size(x0);
	matrix l(n[0], 1), p(n[0], 1), s(s0), D = unit_mat(n[0]);
	solution X, Xt;
	X.x = x0;
	X.fit_fun();
	while (true)
	{
		for (int i = 0; i < n[0]; ++i)
		{
			Xt.x = X.x + s(i) * D[i];
			Xt.fit_fun();
			if (Xt.y < X.y)
			{
				X = Xt;
				l(i) += s(i);
				s(i) *= alfa;
			}
			else
			{
				++p(i);
				s(i) *= -beta;
			}
		}
		bool change = true;
		for (int i = 0; i < n[0]; ++i)
			if (l(i) == 0 || p(i) == 0)
			{
				change = false;
				break;
			}
		if (change)
		{
			matrix Q(n[0], n[0]), v(n[0], 1);
			for (int i = 0; i < n[0]; ++i)
				for (int j = 0; j <= i; ++j)
					Q(i, j) = l(i); //l to lambda
			Q = D * Q;
			v = Q[0] / norm(Q[0]);
			D = set_col(D, v, 0);
			for (int i = 1; i < n[0]; ++i)
			{
				matrix temp(n[0], 1);
				for (int j = 0; j < i; ++j)
					temp = temp + trans(Q[i]) * D[j] * D[j];
				v = (Q[i] - temp) / norm(Q[i] - temp);
				D = set_col(D, v, i);
			}
			s = s0;
			l = matrix(n[0], 1);
			p = matrix(n[0], 1);
		}
		double max_s = abs(s(0));
		for (int i = 1; i < n[0]; ++i)
			if (max_s < abs(s(i)))
				max_s = abs(s(i));
		//cout << X.x(0) << " " << X.x(1) << endl;
		if (max_s < epsilon || solution::f_calls>Nmax)
			return X;
	}
}
#endif
#if LAB_NO>3
solution pen(matrix x0, double c, double dc, double epsilon, int Nmax, matrix O)
{
	double alfa = 1, beta = 0.5, gama = 2, delta = 0.5, s = 0.5;
	matrix A(new double[2]{ c,O(0) }, 2);
	solution X, X1;
	X.x = x0;
	while (true)
	{
		X1 = sym_NM(X.x, s, alfa, beta, gama, delta, epsilon, Nmax, A);
		if (norm(X.x - X1.x) < epsilon || solution::f_calls > Nmax)
			return X1;
		A(0) *= dc;
		X = X1;
	}
}

solution sym_NM(matrix x0, double s, double alfa, double beta, double gama, double delta, double epsilon, int Nmax, matrix O)
{
	int* n = get_size(x0);
	matrix D = unit_mat(n[0]);
	int N = n[0] + 1;
	solution* S = new solution[N];
	S[0].x = x0;
	S[0].fit_fun(O);
	for (int i = 1; i < N; ++i)
	{
		S[i].x = S[0].x + s * D[i - 1];
		S[i].fit_fun(O);
	}
	solution p_o, p_e, p_z;
	matrix p_sr;
	int i_min, i_max;
	while (true)
	{
		i_min = i_max = 0;
		for (int i = 1; i < N; ++i)
		{
			if (S[i_min].y > S[i].y)
				i_min = i;
			if (S[i_max].y < S[i].y)
				i_max = i;
		}
		p_sr = matrix(n[0], 1);
		for (int i = 0; i < N; ++i)
			if (i != i_max)
				p_sr = S[i].x;
		p_sr = p_sr / (N - 1.0);
		p_o.x = p_sr + alfa * (p_sr - S[i_max].x);
		p_o.fit_fun(O);
		if (S[i_min].y <= p_o.y && p_o.y < S[i_max].y)
			S[i_max] = p_o;
		else if (p_o.y < S[i_min].y)
		{
			p_e.x = p_sr + gama * (p_o.x - p_sr);
			p_e.fit_fun(O);
			if (p_e.y < p_o.y)
				S[i_max] = p_e;
			else
				S[i_max] = p_o;
		}
		else
		{
			p_z.x = p_sr + beta * (S[i_max].x - p_sr);
			p_z.fit_fun(O);
			if (p_z.y < S[i_max].x - p_sr)
				S[i_max] = p_z;
			else
			{
				for (int i = 0; i < N; ++i)
					if (i != i_min)
					{
						S[i].x = delta * (S[i].x - S[i_min].x);
						S[i].fit_fun(O);
					}
			}
		}
		double max_s = norm(S[0].x - S[i_min].x);
		for (int i = 1; i < N; ++i)
			if (max_s < norm(S[i].x - S[i_min].x))
				max_s = norm(S[i].x - S[i_min].x);
		if (max_s<epsilon || solution::f_calls>Nmax)
			return S[i_min];
	}
}
#endif
#if LAB_NO>4
solution SD(matrix x0, double h0, double epsilon, int Nmax, matrix O)
{
	int* n = get_size(x0);
	solution X, X1;
	X.x = x0;
	matrix d(n[0], 1), P(n[0], 2), limits = O;
	solution h;
	double b;
	while (true)
	{
		X.grad();
		d = -X.g;
		P = set_col(P, X.x, 0);
		P = set_col(P, d, 1);
		if (h0 < 0)
		{
			b = compute_b(X.x, d, limits);
			h = golden(0, b, epsilon, Nmax, P);
			X1.x = X.x + h.x * d;
		}
		else
		{
			X1.x = X.x + h0 * d;
		}
		if (solution::f_calls > Nmax || solution::g_calls > Nmax || norm(X1.x - X.x) < epsilon)
		{
			X1.fit_fun();
			return X1;
		}
		X = X1;
	}
}

solution CG(matrix x0, double h0, double epsilon, int Nmax, matrix O)
{
	int* n = get_size(x0);
	solution X, X1;
	X.x = x0;
	matrix d(n[0], 1), P(n[0], 2), limits = O;
	solution h;
	double b, beta;
	X.grad();
	d = -X.g;
	while (true)
	{
		P = set_col(P, X.x, 0);
		P = set_col(P, d, 1);
		if (h0 < 0) {
			b = compute_b(X.x, d, limits);
			h = golden(0, b, epsilon, Nmax, P);
			X1.x = X.x + h.x * d;
		}
		else
		{
			X1.x = X.x + h0 * d;
		}
		if (solution::f_calls > Nmax || solution::g_calls > Nmax || norm(X1.x - X.x) < epsilon)
		{
			X1.fit_fun();
			return X1;
		}
		X1.grad();
		beta = pow(norm(X1.g), 2) / pow(norm(X.g), 2);
		d = -X1.g + beta * d;
		X = X1;
	}
}

solution Newton(matrix x0, double h0, double epsilon, int Nmax, matrix O)
{
	int* n = get_size(x0);
	solution X, X1;
	X.x = x0;
	matrix d(n[0], 1), P(n[0], 2), limits = O;
	solution h;
	double b;
	while (true)
	{
		X.grad();
		X.hess();
		d = -inv(X.H) * X.g;
		P = set_col(P, X.x, 0);
		P = set_col(P, d, 1);
		if (h0 < 0)
		{
			b = compute_b(X.x, d, limits);
			h = golden(0, b, epsilon, Nmax, P);
			X1.x = X.x + h.x * d;
		}
		else
		{
			X1.x = X.x + h0 * d;
		}
		if (solution::f_calls > Nmax || solution::g_calls > Nmax || norm(X1.x - X.x) < epsilon)
		{
			X1.fit_fun();
			return X1;
		}
		X = X1;
	}
}

solution golden(double a, double b, double epsilon, int Nmax, matrix O)
{
	double alfa = (sqrt(5.0) - 1) / 2;
	solution A, B, C, D;
	A.x = a;
	B.x = b;
	C.x = B.x - alfa * (B.x - A.x);
	C.fit_fun(O);
	D.x = A.x + alfa * (B.x - A.x);
	D.fit_fun(O);
	while (true)
	{
		if (C.y < D.y)
		{
			B = D;
			D = C;
			C.x = B.x - alfa * (B.x - A.x);
			C.fit_fun(O);
		}
		else
		{
			A = C;
			C = D;
			D.x = A.x + alfa * (B.x - A.x);
			D.fit_fun(O);
		}
		if (B.x - A.x<epsilon || solution::f_calls>Nmax)
		{
			A.x = (A.x + B.x) / 2.0;
			A.fit_fun(O);
			return A;
		}
	}
}

double compute_b(matrix x, matrix d, matrix limits)
{
	int* n = get_size(x);
	double b = 1e9, bi;
	for (int i = 0; i < n[0]; ++i)
	{
		if (d(i) == 0)
			bi = 1e9;
		else if (d(i) > 0)
			bi = (limits(i, 1) - x(i)) / d(i);
		else
			bi = (limits(i, 0) - x(i)) / d(i);
		if (b > bi)
			b = bi;
	}
	return b;
}
#endif
#if LAB_NO>5
solution Powell(matrix x0, double epsilon, int Nmax, matrix O)
{
	int* n = get_size(x0);
	matrix D = unit_mat(n[0]), A(n[0], 3), limits(n[0], 2);
	limits = set_col(limits, O[0], 0);
	limits = set_col(limits, O[1], 1);
	A(0, 2) = O(0, 2);
	solution X, P, h;
	X.x = x0;
	double* ab;
	while (true)
	{
		P = X;
		for (int i = 0; i < n[0]; ++i)
		{
			A = set_col(A, P.x, 0);
			A = set_col(A, D[i], 1);
			ab = compute_ab(P.x, D[i], limits);
			h = golden(ab[0], ab[1], epsilon, Nmax, A);
			P.x = P.x + h.x * D[i];
		}
		if (norm(P.x - X.x) < epsilon || solution::f_calls > Nmax)
		{
			P.fit_fun();
			return P;
		}
		for (int i = 0; i < n[0] - 1; ++i) //wszystkie kierunki przestawiamy w lewo
			D = set_col(D, D[i + 1], 1);
		D = set_col(D, P.x - X.x, n[0] - 1); //ostatni kierunek to kierunek średni
		A = set_col(A, P.x, 0);
		A = set_col(A, D[n[0] - 1], 1);
		ab = compute_ab(P.x, D[n[0] - 1], limits);
		h = golden(ab[0], ab[1], epsilon, Nmax, A);
		X.x = P.x + h.x * D[n[0] - 1];
	}
}

double* compute_ab(matrix x, matrix d, matrix limits)
{
	int* n = get_size(x);
	double* ab = new double[2]{ -1e9,1e9 };
	double ai, bi;
	for (int i = 0; i < n[0]; ++i)
	{
		if (d(i) == 0)
		{
			ai = -1e9;
			bi = 1e9;
		}
		else if (d(i) > 0)
		{
			ai = (limits(i, 0) - x(i)) / d(i);
			bi = (limits(i, 1) - x(i)) / d(i);
		}
		else
		{
			ai = (limits(i, 1) - x(i)) / d(i);
			bi = (limits(i, 0) - x(i)) / d(i);
		}
		if (ab[0] < ai)
			ab[0] = ai;
		if (ab[1] > bi)
			ab[1] = bi;
	}
	return ab;
}
#endif
#if LAB_NO>6
solution EA(int N, matrix limits, double epsilon, int Nmax, matrix O)
{
	//N - liczba zmiennych decyzyjnych
	//limits - ograniczenia
	//epsilon - dok�adno��
	//Nmax - maksymalna liczba wywo�a� funckji celu
	//O - macierz 1x1 zawieraj�ca warto�� sigma

	int mi = 20, lambda = 40;
	solution *P = new solution[mi + lambda];
	solution *Pm = new solution[mi];
	random_device rd;
	default_random_engine gen;
	gen.seed(static_cast<unsigned int>(chrono::system_clock::now().time_since_epoch().count()));
	normal_distribution<double> distr(0.0, 1.0);
	matrix IFF(mi, 1), temp(N, 2);
	double r, s, s_IFF;
	double tau = pow(2 * N, -0.5), tau1 = pow(2 * pow(N, 0.5), -0.5);
	int j_min;
	for (int i = 0; i < mi; ++i)
	{
		P[i].x = matrix(N, 2);
		for (int j = 0; j < N; ++j)
		{
			P[i].x(j, 0) = (limits(j, 1) - limits(j, 0))*rd() / rd.max() + limits(j, 0);
			P[i].x(j, 1) = O(0);
		}
		P[i].fit_fun();
		if (P[i].y < epsilon)
			return P[i];
	}
	while (true)
	{
		s_IFF = 0;
		for (int i = 0; i < mi; ++i)
		{
			IFF(i) = 1 / P[i].y(0);
			s_IFF += IFF(i);
		}
		for (int i = 0; i < lambda; ++i)
		{
			r = s_IFF * rd() / rd.max();
			s = 0;
			for (int j = 0; j < mi; ++j)
			{
				s += IFF(j);
				if (r <= s)
				{
					P[mi + i] = P[j];
					break;
				}
			}
		}
		for (int i = 0; i < lambda; ++i)
		{
			r = distr(gen);
			for (int j = 0; j < N; ++j)
			{
				P[mi + i].x(j, 1) *= exp(tau1*r + tau * distr(gen));
				P[mi + i].x(j, 0) += P[mi + i].x(j, 1)*distr(gen);
			}
		}
		for (int i = 0; i < lambda; i += 2)
		{
			r = 1.0*rd() / rd.max();
			temp = P[mi + i].x;
			P[mi + i].x = r * P[mi + i].x + (1 - r)*P[mi + i + 1].x;
			P[mi + i + 1].x = r * P[mi + i + 1].x + (1 - r)*temp;
		}
		for (int i = 0; i < lambda; ++i)
		{
			P[mi + i].fit_fun();
			if (P[mi + i].y < epsilon)
				return P[mi + i];
		}
		for (int i = 0; i < mi; ++i)
		{
			j_min = 0;
			for (int j = 1; j < mi + lambda; ++j)
				if (P[j_min].y > P[j].y)
					j_min = j;
			Pm[i] = P[j_min];
			P[j_min].y = 1e10;
		}
		for (int i = 0; i < mi; ++i)
			P[i] = Pm[i];
		if (solution::f_calls > Nmax)
			return P[0];
	}
}
#endif
