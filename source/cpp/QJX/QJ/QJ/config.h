#pragma once
#include "config.h"
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

#define PI 3.1415926535897932384626433832795
#define EPS 1.0e-14

using namespace std;

struct RunParam
{
	int		task;				// Task id	

	int		debug;				// Dump debug info?

	int		issmtx;				// Save matrixes?

	int		ipp;				// Print progress?

	string	fn_data;			// File name with main data

	string	fn_aux;				// File name with aux data

	string	path;				// Path to files (if nessesary)

	RunParam(
		int _task = 0,

		int	_debug = 0,

		int _issmtx = 0,

		int _ipp = 0,

		string _init_file = "",

		string _path = ""
	)
	{
		task = _task;

		debug = _debug;

		issmtx = _issmtx;

		ipp = _ipp;

		init_file = _init_file;

		path = _path;
	}
};

struct ConfigParam
{
	int		N;						// Num of spins

	int		dt;						// Dissipator type
	double	dp;						// Dissipator phase
	double	g;						// Dissipation

	double J;						// Disorder J
	double h_z;						// Disorder z
	double h_x;						// Disorder x

	int num_steps_t_0;				// Number of integration steps per_period
	int num_steps_t_1;				// Number of integration steps per_period
	int num_periods;				// Number of periods
	double t_0;						// Subperiod t_0
	double t_1;						// Subperiod t_1

	int seed;						// Seed
	int max_num_seeds;				// Maximum number of seeds

	int			int_ist;			// Init state type
	int			int_isi;			// Init state index
	int			int_dt;				// Dump type
	int			int_dn;				// Num dumps

	ConfigParam(
		int _N = 6,

		int	_dt = 1,
		double _dp = 0.0,
		double _g = 0.1,

		double _J = 1.0,
		double _h_z = 1.0,
		double _h_x = 0.0,

		int _num_steps_t_0 = 1000,
		int _num_steps_t_1 = 1000,
		int _num_periods = 1000,
		double _t_0 = 1.5 * PI,
		double _t_1 = 0.5 * PI,

		int _seed = 0,
		int _max_num_seeds = 1000000,

		int	_int_ist = 0,
		int	_int_isi = 0,
		int	_int_dt = 1,
		int	_int_dn = 1
	)
	{
		N = _N;

		dt = _dt;
		dp = _dp;
		g = _g;

		J = _J;
		h_z = _h_z;
		h_x = _h_x;

		num_steps_t_0 = _num_steps_t_0;
		num_steps_t_1 = _num_steps_t_1;
		num_periods = _num_periods;
		t_0 = _t_0;
		t_1 = _t_1;

		seed = _seed;
		max_num_seeds = _max_num_seeds;

		int_ist = _int_ist;
		int_isi = _int_isi;
		int_dt = _int_dt;
		int_dn = _int_dn;
	}
};


void set_param(RunParam &rp, ConfigParam &cp, string str, string val);

void init_params(RunParam &rp, ConfigParam &cp, char * file_name);

void output_params(RunParam &rp, ConfigParam &param);
