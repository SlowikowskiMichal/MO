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
		matrix* Y = solve_ode(0.0, 1.0, 1000.0, Y0);
		ofstream S("sim_t.csv");
		S << Y[0];
		S.close();
		S.open("sim_v.csv");
		S << Y[1];
		S.close();
#elif LAB_NO==2
#if LAB_PART == 1
#if TASK == 1 // wywo³ania po 100 razy dla danego alfa
		int counter = 0;
		fstream file("wynik.txt", ios::app);
		while (counter < 100) {
			double x0, d = 1.0, alfa = 1.5, epsilon = 1e-5, gamma = 1e-200;
			int Nmax = 1000;
			random_device R;
			x0 = 200.0 * R() / R.max() - 100;
			file << x0 << ";";
			cout << x0 << endl << endl;
			double* ab = expansion(x0, d, alfa, Nmax);
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
#elif TASK == 2 // wywo³anie bez zawê¿ania przedzia³u
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
#elif TASK == 3 // kod z æwiczeñ bez wiêkszych zmian
		double x0, d = 1.0, alfa = 1.5, epsilon = 1e-5, gamma = 1e-200;
		int Nmax = 1000;
		random_device R;
		x0 = 200.0 * R() / R.max() - 100;
		cout << x0 << endl << endl;
		double* ab = expansion(x0, d, alfa, Nmax);
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
		x0 = (1e-2 - 1e-4) * R() / R.max() + 1e-4;
		cout << x0 << endl << endl;
		double* ab = expansion(x0, d, alfa, Nmax);
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
		matrix * Y = solve_ode(0.0, 1.0, 1000.0, Y0);
		ofstream S("sim_t2.csv");
		S << Y[0];
		S.close();
		S.open("sim_v2.csv");
		S << Y[1];
		S.close();
#endif
#endif

#elif LAB_NO==3
#if LAB_PART==1
		matrix x0(2, 1);
		double alfa, beta, epsilon = 1e-3, s = 0.1;
		random_device R;
		int Nmax = 100;
		x0(0) = -0.203577;


			;//2.0 * R() / R.max() - 1;
		x0(1) = 0.272492;//2.0 * R() / R.max() - 1;
		//cout << x0 << endl << endl;
		alfa = 0.5;
		//cout << "----------------------------"<<endl;
		//cout << x0(0) << " " << x0(1) <<endl;
		solution opt_HJ = HJ(x0, s, alfa, epsilon, Nmax);
		//cout << opt_HJ << endl << endl;
		//cout << "----------------------------" << endl;
		solution::clear_calls();

		alfa = 2.0;
		beta = 0.5;
		matrix s0(2, 1);
		s0(0) = s;
		s0(1) = s;
		//cout << "----------------------------" << endl;
		//cout << x0(0) << " " << x0(1) << endl;
		solution opt_R = Rosen(x0, s0, alfa, beta, epsilon, Nmax);
		//cout << opt_R << endl;
		//cout << "----------------------------" << endl;
		solution::clear_calls();
#elif LAB_PART==2
		/*matrix x0(2, 1);
		double alfa, beta, epsilon = 1e-3, s = 0.5;
		random_device R;
		int Nmax = 100;
		x0(1) = 3.49888;//10.0 * R() / R.max();
		x0(0) = 8.263;//10.0 * R() / R.max();
		cout << x0 << endl << endl;
		alfa = 0.5;
		solution opt_HJ = HJ(x0, s, alfa, epsilon, Nmax);
		cout << opt_HJ << endl << endl;
		solution::clear_calls();

		alfa = 2.0;
		beta = 0.5;
		matrix s0(2, 1);
		s0(0) = s;
		s0(1) = s;
		solution opt_R = Rosen(x0, s0, alfa, beta, epsilon, Nmax);
		cout << opt_R << endl;
		solution::clear_calls();
		*/
		matrix x0(2, 1);
		double alfa, beta, epsilon = 1e-3, s = 0.5;
		x0(1) = 2.77232;//10.0 * R() / R.max();
		x0(0) = 3.18488;//10.0 * R() / R.max();
		solution X;
		X.x = x0;
		X.fit_fun();

#endif
#elif LAB_NO==4

#elif LAB_NO==5

#elif LAB_NO==6

#elif LAB_NO==7

#endif
	}
	catch (char* EX_INFO)
	{
		cout << EX_INFO << endl;
	}
	system("pause");
	return 0;
}