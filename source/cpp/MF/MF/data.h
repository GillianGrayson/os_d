#pragma once
#include "config.h"
#include "utils.h"

struct MainData
{
	double step;			// Integration step

	int size;				// Data size

	double time;			// Current intrgration time

	double * data;			// Integration data
	double * args;			// Integration routines data
	double * k1s;			// Integration routines data
	double * k2s;			// Integration routines data
	double * k3s;			// Integration routines data
	double * k4s;			// Integration routines data

	int num_lpn;			// Number of Lyapunov exponents

	double ** data_lpn;		// Data for Lyapunov exps
	double ** args_lpn;		// Integration Lyapunov routines data
	double ** k1s_lpn;		// Integration Lyapunov routines data
	double ** k2s_lpn;		// Integration Lyapunov routines data
	double ** k3s_lpn;		// Integration Lyapunov routines data
	double ** k4s_lpn;		// Integration Lyapunov routines data
	double * norms_lpn;		// Norms for Lyapunov vectors
	double * exps_lpn;		// Lyapunov exponents
	double * rvm_lpn;		// Aux data for Lyapunov exponents

	int cd_ps;				// Points shift for correlated dimension
	int cd_M;				// Num points for correlated dimension
	int cd_size;			// Dimenson of reconstructed data
	double cd_obs;			// Observable
	int * cd_ti;			// Time iterations
	double ** cd_rd;		// Reconstruction data
	double cd_ci;			// Correlated integral

	int num_dumps;
	double ** data_evo;
	int dump_id;
};


void init_main_data(ConfigParam &cp, MainData &md);

void init_evo_data(ConfigParam &cp, MainData &md);

void init_lpn_data(ConfigParam &cp, MainData &md);

void init_cd_data(ConfigParam &cp, MainData &md);

void init_cd_d_data(ConfigParam &cp, MainData &md);

void init_cd_sd_data(ConfigParam &cp, MainData &md);

void print_cd_info(MainData &md);

void delete_main_data(MainData &dt);

void delete_evo_data(MainData &md);

void delete_lpn_data(MainData &md);

void delete_cd_data(MainData &md);

void delete_cd_d_data(MainData &md);

void delete_cd_sd_data(MainData &md);

void init_cond(RunParam &rp, ConfigParam &cp, MainData &md);

void init_cond_lpn(ConfigParam &cp, MainData &md);

void calc_norm_lpn(MainData &md, int lpn_id);

void normalization_lpn(MainData &md, int lpn_id);

void scalar_mult_lpn(MainData &md, double * mults, int lpn_id, int lpn_id_tmp);

void sub_lpn(MainData &md, double* mults, int lpn_id, int lpn_id_tmp);

void gsorth_lpn(MainData &md);

void calc_ci(ConfigParam &cp, MainData &md);

