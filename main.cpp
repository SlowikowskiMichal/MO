#include<iostream>
#include<random>
#include<chrono>
#include<fstream>
#include"opt_alg.h"
#include"ode_solver.h"

using namespace std;

int main()
{
	try
	{
		cout << "LAB NUMBER " << LAB_NO << endl;
		cout << "LAB PART " << LAB_PART << endl;
		cout << "TASK " << TASK << endl << endl;
#if LAB_NO==1
		matrix Y0 = matrix(new double[3]{ 5,1,10 }, 3);
		matrix *Y = solve_ode(0.0, 1.0, 1000.0, Y0);
		ofstream S("sim_t.csv");
		S << Y[0];
		S.close();
		S.open("sim_v.csv");
		S << Y[1];
		S.close();
#elif LAB_NO==2
#if LAB_PART == 1
#if TASK == 1 // wywo�ania po 100 razy dla danego alfa
		int counter = 0;
		fstream file("wynik.txt", ios::app);
		while (counter < 100) {
			double x0, d = 1.0, alfa = 1.5, epsilon = 1e-5, gamma = 1e-200;
			int Nmax = 1000;
			random_device R;
			x0 = 200.0*R() / R.max() - 100;
			file << x0 << ";";
			cout << x0 << endl << endl;
			double *ab = expansion(x0, d, alfa, Nmax);
			cout << ab[0] << '\t' << ab[1] << endl << endl;
			file << ab[0] << ";" << ab[1] << ";" << solution::f_calls << ";";
			solution::clear_calls();
			solution opt_F = fib(ab[0], ab[1], epsilon);
			cout << opt_F << endl << endl;
			file << opt_F.x(0) << ";" << opt_F.y(0) << ";" << solution::f_calls << ";";
			solution::clear_calls();
			solution opt_L = lag(ab[0], ab[1], epsilon, gamma, Nmax);
			cout << opt_L << endl << endl;
			file << opt_L.x(0) << ";" << opt_L.y(0) << ";" << solution::f_calls;
			solution::clear_calls();
			file << endl;
			counter++;
		}
		file.close();
#elif TASK == 2 // wywo�anie bez zaw�ania przedzia�u
		double d = 1.0, alfa = 1.5, epsilon = 1e-5, gamma = 1e-200;
		int Nmax = 1000;
		random_device R;
		solution::clear_calls();
		solution opt_F = fib(-100, 100, epsilon);
		cout << opt_F << endl << endl;
		solution::clear_calls();
		solution opt_L = lag(-100, 100, epsilon, gamma, Nmax);
		cout << opt_L << endl << endl;
		solution::clear_calls();
#elif TASK == 3 // kod z �wicze� bez wi�kszych zmian
		double x0, d = 1.0, alfa = 1.5, epsilon = 1e-5, gamma = 1e-200;
		int Nmax = 1000;
		random_device R;
		x0 = 200.0*R() / R.max() - 100;
		cout << x0 << endl << endl;
		double *ab = expansion(x0, d, alfa, Nmax);
		cout << ab[0] << '\t' << ab[1] << endl << endl;
		solution::clear_calls();
		solution opt_F = fib(ab[0], ab[1], epsilon);
		cout << opt_F << endl << endl;
		solution::clear_calls();
		solution opt_L = lag(ab[0], ab[1], epsilon, gamma, Nmax);
		cout << opt_L << endl << endl;
		solution::clear_calls();
#endif
#elif LAB_PART == 2 
#if TASK == 4 // szukanie pola przekroju DA
		double x0, d = 1e-4, alfa = 1.5, epsilon = 1e-10, gamma = 1e-200;
		int Nmax = 1000;
		random_device R;
		x0 = (1e-2 - 1e-4)*R() / R.max() + 1e-4;
		cout << x0 << endl << endl;
		double *ab = expansion(x0, d, alfa, Nmax);
		cout << ab[0] << '\t' << ab[1] << endl << endl;
		solution::clear_calls();
		solution opt_F = fib(ab[0], ab[1], epsilon);
		cout << opt_F << endl << endl;
		solution::clear_calls();
		solution opt_L = lag(ab[0], ab[1], epsilon, gamma, Nmax);
		cout << opt_L << endl << endl;
		solution::clear_calls();
#elif TASK == 5 // symlacja dla wyszukanego DA w TASK == 4
		matrix Y0 = matrix(new double[3]{ 5,1,10 }, 3);
		matrix *Y = solve_ode(0.0, 1.0, 1000.0, Y0);
		ofstream S("sim_t2.csv");
		S << Y[0];
		S.close();
		S.open("sim_v2.csv");
		S << Y[1];
		S.close();
#endif
#endif

#elif LAB_NO==3

#elif LAB_NO==4

#elif LAB_NO==5
		
#elif LAB_NO==6

#elif LAB_NO==7

#endif
	}
	catch (char * EX_INFO)
	{
		cout << EX_INFO << endl;
	}
	system("pause");
	return 0;
}