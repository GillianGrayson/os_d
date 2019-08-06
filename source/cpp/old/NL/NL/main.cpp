#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<fstream>
#include<time.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <complex>
#include <mkl.h>
#include <omp.h>
#include <numeric>
#include <algorithm>
#include <math.h>
using namespace std;

long const N = 1; //number of units
long const n = 3 * N; // dimension of the system
double const pi = 3.141592653589793238463643;
int Nstep = 1000;

double gamma = 0.1, J = 1.0, U = 0.0, E0 = 1.0, w = 1, A0 = 1.5, phi = 0;
double h = 2 * pi / w / Nstep, h1, t_tr = 2000 * 2 * pi / w, t_fin = 2000 * 2 * pi / w + t_tr;
double t = 0;

inline void fun(double[n], double[n]);
void stepmake4(double[n], double[n]);
double k1[n], k2[n], k3[n], k4[n], newpoint[n];

string file_name_suffix(int seed, double U)
{
	stringstream fns;
	fns << "_U(" << setprecision(4) << fixed << U << ")";
	fns << "_seed(" << seed << ")";
	fns << ".txt";

	return fns.str();
}

//--------------------------------------------
//  main function
//--------------------------------------------
int main()
{
	long i, j, k, l;
	double flag, flag1;
	double x[n], dx[n];
	int count = 0;

	double U_real = 0.50;
	U = 4.0 * U_real;

	clock_t t0, t1;
	double cputime;
	time_t tstart, tend;

	time_t tm;
	srand((unsigned)time(&tm));

	tstart = time(NULL);
	printf("Experiment started at %s\n", ctime(&tstart));
	tstart = time(NULL);
	t0 = clock();

	for (U = 0.36; U <= 0.3600001; U += 0.03)
	{
		U_real = U / 4.0;

		printf("U: %0.4le\n", U_real);

		for (int seed = 0; seed < 10; seed++)
		{
			t = 0;
			count = 0;

			string file_name = "data" + file_name_suffix(seed, U_real);
			ofstream ofs(file_name);
			if (ofs.is_open())
			{
				ofs << setprecision(16) << scientific;

				VSLStreamStatePtr stream;
				vslNewStream(&stream, VSL_BRNG_MCG31, 77778888);
				vslLeapfrogStream(stream, seed, 1000000);
				vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, n, x, -pi, pi);

				for (i = 0; i < n; i++)
				{
					dx[i] = 0.;
				}

				x[2] = 0; // phi / w;

				while (t < t_fin)
				{
					stepmake4(x, dx);
					for (i = 0; i < n; i++)
					{
						x[i] += dx[i];
					}

					count++;
					t += h;

					if (count == Nstep)
					{
						count = 0;

						if (t > t_tr)
						{
							for (i = 0; i < n; i++)
							{
								ofs << x[i] << " ";
							}
							ofs << endl;
						}
					}
				}

				ofs.close();
			}
		}
	}

	tend = time(NULL);
	printf("Experiment finished at %s\n", ctime(&tend));

	return 0;
}
//--------------------------------------------
// end of main function
//--------------------------------------------
//--------------------------------------------
//  kernel function
//--------------------------------------------
void fun(double ff[], double x[])
{
	double ft, rem;

	ft = sin(w * x[2] + phi);
	ft = A0 * ft;
	//ft = A0 * (2 * (ft > 0) - 1);
	ff[0] = 2.0 * J*sin(x[1]) + 4.0 * gamma*cos(x[1])*cos(x[0]);
	ff[1] = 2.0 * J*cos(x[0])*cos(x[1]) / sin(x[0]) - 2.0 * E0 - 2.0 * ft + U*cos(x[0]) - 4.0 * gamma*sin(x[1]) / sin(x[0]);
	ff[2] = 1.0;

}

//--------------------------------------------
// end of kernel function
//--------------------------------------------

//--------------------------------------------
// RK 4-th order
//--------------------------------------------
void stepmake4(double x[], double dx[])
{

	long i;


	fun(k1, x);
	for (i = 0; i < n; i++) newpoint[i] = x[i] + h*k1[i] / 2;

	fun(k2, newpoint);
	for (i = 0; i < n; i++) newpoint[i] = x[i] + h*k2[i] / 2;

	fun(k3, newpoint);
	for (i = 0; i < n; i++) newpoint[i] = x[i] + h*k3[i];

	fun(k4, newpoint);

	for (i = 0; i < n; i++) dx[i] = h*(k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6;
}
//--------------------------------------------
// end of RK
//--------------------------------------------


int old_main()
{
	long i, j, k, l;
	double count = 0, flag, flag1;
	double x[n], dx[n];

	clock_t t0, t1;
	double cputime;
	time_t tstart, tend;

	char char_out[] = "out.txt";
	char char_out1[] = "out1.txt";
	char data_out[] = "data.txt";
	char inifile[] = "ini.dat";
	fstream initial;
	ofstream output_file;
	ofstream output_file1;
	ofstream data_file;

	data_file.open(data_out, ios::out);


	output_file.open(char_out, ios::out);
	output_file1.open(char_out1, ios::out);
	initial.precision(16);
	output_file.precision(8);

	time_t tm;
	srand((unsigned)time(&tm));

	tstart = time(NULL);
	printf("Experiment started at %s\n", ctime(&tstart));
	tstart = time(NULL);
	t0 = clock();

	for (U = 0.0; U < 3; U += 0.02)
	{
		t = 0;
		count = 0;

		for (i = 0; i < n; i++)
		{
			x[i] = ((double)rand() / RAND_MAX - 0.5);
			dx[i] = 0.;
		}

		x[2] = 0; // phi / w;

		while (t < t_fin)
		{
			stepmake4(x, dx);
			for (i = 0; i < n; i++)
			{
				x[i] += dx[i];
			}

			count++;
			t += h;

			if (count == Nstep)
			{
				count = 0;

				if (t > t_tr)
				{
					output_file << J << " " << E0 << " " << U << " " << gamma << " " << A0 << " " << w << " " << t << " ";
					for (i = 0; i < n; i++)
					{
						output_file << x[i] << " ";
					}
					output_file << endl;
				}
			}
		}
	}

	tend = time(NULL);
	printf("Experiment finished at %s\n", ctime(&tend));

	output_file.close();
	output_file1.close();
	data_file.close();

	return 0;

}