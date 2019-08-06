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

#define TASK_ID_STD 0
#define TASK_ID_DEEP 1


using namespace std;

struct RunParam
{
	int		task;				// Task id	
	int		debug;				// Dump debug info?
	int		issmtx;				// Save matrixes?
	int		ipp;				// Print progress?

	string	init_file;			// Init file (if nessesary)

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
	int		N;						// System size

	int		dt;						// Dissipator type
	double	dp;						// Dissipator phase
	double	g;						// Dissipation
	double	g_add;					// Additional dissipation

	double drv_ampl;				// drv_ampl
	double prm_alpha;				// prm_alpha

	int num_steps_t_0;				// Number of integration steps per_period
	int num_steps_t_1;				// Number of integration steps per_period
	int num_periods_trans;			// Number of periods in trans process
	int num_periods_obser;			// Number of periods in obser process
	double t_0;						// Subperiod t_0
	double t_1;						// Subperiod t_1

	int seed;						// Seed
	int max_num_seeds;				// Maximum number of seeds

	int			int_ist;			// Init state type
	int			int_isi;			// Init state index
	int			int_dt;				// Dump type
	int			int_dn;				// Num dumps

	ConfigParam(
		int _N = 200,

		int	_dt = 1,
		double _dp = 0.0,
		double _g = 0.1,
		double _g_add = 0.01,

		double _drv_ampl = 3.2,
		double _prm_alpha = 5.0,

		int _num_steps_t_0 = 1000,
		int _num_steps_t_1 = 1000,
		int _num_periods_trans = 1000,
		int _num_periods_obser = 1000,
		double _t_0 = 0.98,
		double _t_1 = 1.0,

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
		g_add = _g_add;

		drv_ampl = _drv_ampl;
		prm_alpha = _prm_alpha;

		num_steps_t_0 = _num_steps_t_0;
		num_steps_t_1 = _num_steps_t_1;
		num_periods_trans = _num_periods_trans;
		num_periods_obser = _num_periods_obser;
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
