#pragma once
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
#include <map>

#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

#define DIMER_SYS_ID 0

#define LPN_TASK_ID 0 
#define STD_TASK_ID 1 
#define CD_TASK_ID 2 
#define SIGMA_TASK_ID 3

#define PI 3.1415926535897932384626433832795
#define EPS 1.0e-14

using namespace std;

struct RunParam
{
	int		sys_id;				// System id
	int		task_id;			// Task id

	int		is_debug;			// Dump debug info?
	int		is_pp;				// Print progress?

	string	init_fn;			// Name of init file (if nessesary)
	string	path;				// Path to files (if nessesary)

	int		num_threads;		// Number of threads

	RunParam(
		int _sys_id = 0,
		int	_task_id = 0,

		int _is_debug = 0,
		int _is_pp = 0,

		string _init_fn = "",
		string _path = "",

		int _num_threads = 1
	)
	{
		sys_id = _sys_id;
		task_id = _task_id;

		is_debug = _is_debug;
		is_pp = _is_pp;

		init_fn = _init_fn;
		path = _path;

		num_threads = _num_threads;
	}
};

struct ConfigParam
{
	map<string, double> params;					// System params

	string fn_suffix;							// File name suffixes

	int	qj_num_tp_periods;						// Number of periods in trans process
	int	qj_num_obs_periods;						// Number of observable propagation periods
	int	qj_deep;								// Deep
	int qj_num_trajectories;					// Number of trajectories
	int qj_seed;								// Seed 
	int qj_mns;									// Max number of seeds

	int	dump_type;								// Dump type
	int	dump_num;								// Number of dumps
	
	ConfigParam(
		int _qj_num_tp_periods = 0,
		int _qj_num_obs_periods = 10,
		int _qj_deep = 16,
		int _qj_num_trajectories = 1,
		int _qj_seed = 0,
		int _qj_mns = 1000000
	)
	{
		qj_num_tp_periods = _qj_num_tp_periods;
		qj_num_obs_periods = _qj_num_obs_periods;
		qj_deep = _qj_deep;
		qj_num_trajectories = _qj_num_trajectories;
		qj_seed = _qj_seed;
		qj_mns = _qj_mns;
	}
};


void set_param(RunParam * rp, ConfigParam * cp, string str, string val);

void init_params(RunParam * rp, ConfigParam * cp, char * fn_config, char * fn_param);

void output_params(RunParam * rp, ConfigParam * cp);
