#pragma once
#include "Config.h"
#include "utils.h"

struct MainData
{
	int size;
	double step;
	double T;
};

void init_main_data(RunParam &rp, ConfigParam &cp, MainData &md);
void delete_main_data(MainData &md);

struct PropData
{
	int total_num_dumps;

	int * dump_periods;

	double * deep_dump_times;
	double * deep_evals;
	MKL_Complex16 * deep_avg_rho;

	MKL_Complex16 * floquet;

	double* k1;
	double* k2;
	double* k3;
	double* k4;
	double* val;
	double* tmp;
	double* tmp_drv;
};

void f_basis_prop_std(RunParam &rp, ConfigParam &cp, MainData &md);
void f_basis_prop_floquet(RunParam &rp, ConfigParam &cp, MainData &md);

void init_prop_data_std(RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd);
void init_prop_data_floquet(RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd);

void free_prop_data_std(RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd);
void free_prop_data_floquet(RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd);

void dump_prop_data_std(RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd);
void dump_prop_data_floquet(RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd);

