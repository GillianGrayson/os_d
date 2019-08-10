#ifndef __CONFIG__
#define __CONFIG__

#pragma once
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <complex>
#include <omp.h>
#include <numeric>
#include <algorithm>
#include <math.h>

#define TASK_ID_STD 0
#define FLOQUET_ID_BASE 1
#define FLOQUET_ID_F 2

using namespace std;

struct RunParam
{
	int task; // Task id	
	int debug; // Dump debug info?
	int issmtx; // Save matrixes?
	int ipp; // Print progress?

	RunParam(
		int _task = 0,
		int	_debug = 0,
		int _issmtx = 0,
		int _ipp = 0
	)
	{
		task = _task;
		debug = _debug;
		issmtx = _issmtx;
		ipp = _ipp;
	}
};


struct ConfigParam
{
	int N;

	int diss_type;
	double diss_gamma;
	double diss_phase;

	int drv_type;
	double drv_ampl;
	double drv_freq;
	double drv_phase;

	double prm_E;
	double prm_U;
	double prm_J;

	int num_steps;
	int num_periods_trans;
	int num_periods_obser;

	int seed;
	int max_num_seeds;

	int int_ist; // Init state type
	int int_isi; // Init state index
	int int_dt; // Dump type
	int int_dn; // Num dumps

	ConfigParam(
		int _N = 10,
		int _diss_type = 0,
		double _diss_gamma = 0.1,
		double _diss_phase = 0.0,

		int _drv_type = 1,
		double _drv_ampl = 3.4,
		double _drv_freq = 1.0,
		double _drv_phase = 0.0,

		double _prm_E = 0.0,
		double _prm_U = 0.15,
		double _prm_J = 1.0,

		int _num_steps = 10000,
		int _num_periods_trans = 100,
		int _num_periods_obser = 100,

		int _seed = 0,
		int _max_num_seeds = 1000000,

		int _int_ist = 0,
		int _int_isi = 0,
		int _int_dt = 0,
		int _int_dn = 100
	)
	{
		N = _N;

		diss_type = _diss_type;
		diss_gamma = _diss_gamma;
		diss_phase = _diss_phase;

		drv_type = _drv_type;
		drv_ampl = _drv_ampl;
		drv_freq = _drv_freq;
		drv_phase = _drv_phase;

		prm_E = _prm_E;
		prm_U = _prm_U;
		prm_J = _prm_J;

		num_steps = _num_steps;
		num_periods_trans = _num_periods_trans;
		num_periods_obser = _num_periods_obser;

		seed = _seed;
		max_num_seeds = _max_num_seeds;

		int_ist = _int_ist;
		int_isi = _int_isi;
		int_dt = _int_dt;
		int_dn = _int_dn;
	}
};

#endif