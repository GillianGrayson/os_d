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

#define BASIC_EXP_ID 0
#define LPN_FIN_EXP_ID 1
#define BASIC_AND_LPN_FIN_EXP_ID 2
#define CD_EXP_ID 3
#define CD_D_EXP_ID 4
#define CD_SD_EXP_ID 5

using namespace std;

struct RunParam
{
	int		task;				// Task id	

	string	path;				// Path to files (if nessesary)

	double	U_start;			// Start value of interaction
	double	U_shift;			// Shift of interaction
	int		U_num;				// Number of slices U

	int		seed_start;			// Start seed for rng
	int		seed_num;			// Number of seeds
	int		max_num_seeds;		// Maximum number of seeds

	double	cd_eps;				// Start value of eps
	int	cd_eps_num;				// Number of eps
	int	cd_eps_ndpd;			// Number of dumps per decade

	RunParam(
		int _task = 0,

		string _path = "",

		double _U_start = 0.5,
		double _U_shift = 0.1,
		int _U_num = 1,

		int _seed_start = 0,
		int _seed_num = 1,
		int _max_num_seeds = 1000000,

		double _cd_eps = 0.0001,
		int _cd_eps_num = 1,
		int _cd_eps_ndpd = 10
	)
	{
		task = _task;

		path = _path;

		U_start = _U_start;
		U_shift = _U_shift;
		U_num = _U_num;

		seed_start = _seed_start;
		seed_num = _seed_num;
		max_num_seeds = _max_num_seeds;

		cd_eps = _cd_eps;
		cd_eps_num = _cd_eps_num;
		cd_eps_ndpd = _cd_eps_ndpd;
	}
};

struct ConfigParam
{				
	int		mt;					// Modulatation type
	int		num_steps;			// Number of integration steps
	int		npt;				// Number of periods in trans process
	int		np;					// Number of periods
	double	E;					// Bias in local potential
	double	A;					// Amplitude of modulation
	double	omega;				// Frequency of modulation
	double	phase;				// Phase of modulation
	double	gamma;				// Dissipation parameter
	double	J;					// Hopping

	double	cd_eps;				// Epsilon for correlated dimension
	int		cd_dim;				// Dimension of space for correlated dimension
	int		cd_nppp;				// Number of points per period

	double	U;					// Interaction
	int		seed;				// Seed

	double	T;					// Period

	ConfigParam(
		int	_task = 0,
		int	_mt = 0,
		int	_num_steps = 1000,
		int	_npt = 1000,
		int	_np = 1000,
		double _E = 1,
		double _A = 1.5,
		double _omega = 1.0,
		double _phase = 0.0,
		double _gamma = 0.1,
		double _J = 1.0
	)
	{
		mt = _mt;
		num_steps = _num_steps;
		npt = _npt;
		np = _np;
		E = _E;
		A = _A;
		omega = _omega;
		phase = _phase;
		gamma = _gamma;
		J = _J;

		cd_eps = 0.1;
		cd_dim = 2;
		cd_nppp = 100;

		U = 0.0;
		seed = 0.0;

		T = 2.0 * PI / omega;
	}
};


void set_param(RunParam &rp, ConfigParam &cp, string str, string val);

void init_params(RunParam &rp, ConfigParam &cp, char * file_name);

void output_params(RunParam &rp, ConfigParam &cp);
