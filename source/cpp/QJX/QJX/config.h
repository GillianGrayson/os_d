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
#include <utility>

#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

#define DIMER_SYS_ID 0
#define JCS_SYS_ID 1
#define PS_SYS_ID 2

#define LPN_TASK_ID 0
#define STD_TASK_ID 1
#define CD_TASK_ID 2
#define SIGMA_TASK_ID 3
#define STD_DEEP_TASK_ID 4
#define LPN_DEEP_TASK_ID 5
#define LPN_ALL_TASK_ID 6

#define QJ_PROP_TYPE 0
#define RK_PROP_TYPE 1

#define PI 3.1415926535897932384626433832795
#define EPS 1.0e-14

using namespace std;

struct RunParam
{
	int		sys_id;				// System id
	int		task_id;			// Task id
	int		prop_id;			// Propagator id

	int		is_debug;			// Dump debug info?
	int		is_pp;				// Print progress?

	string	init_fn;			// Name of init file (if nessesary)
	string	path;				// Path to files (if nessesary)

	int		num_threads;		// Number of threads

	RunParam(
		int _sys_id = 0,
		int	_task_id = 0,
		int _prop_id = 0,

		int _is_debug = 0,
		int _is_pp = 0,

		string _init_fn = "",
		string _path = "",

		int _num_threads = 1
	)
	{
		sys_id = _sys_id;
		task_id = _task_id;
		prop_id = _prop_id;

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

	int	num_tp_periods;							// Number of periods in trans process
	int	num_obs_periods;						// Number of observable propagation periods
	int num_trajectories;						// Number of trajectories
	int seed;									// Seed 
	int mns;									// Max number of seeds

	int	ex_deep;								// Deep for quantum jumps

	int rk_ns;									// Num steps per (sub-)period for RK

	int	dump_type;								// Dump type
	int	dump_num;								// Number of dumps
	
	ConfigParam(
		int _num_tp_periods = 0,
		int _num_obs_periods = 10,
		int _num_trajectories = 1,
		int _seed = 0,
		int _mns = 1000000,
		int _ex_deep = 16,
		int _rk_ns = 1000
	)
	{
		num_tp_periods = _num_tp_periods;
		num_obs_periods = _num_obs_periods;
		num_trajectories = _num_trajectories;
		seed = _seed;
		mns = _mns;
		ex_deep = _ex_deep;
		rk_ns = _rk_ns;
	}
};


void set_param(RunParam * rp, ConfigParam * cp, string str, string val);

void init_params(RunParam * rp, ConfigParam * cp, char * fn_config, char * fn_param);

void output_params(RunParam * rp, ConfigParam * cp);
