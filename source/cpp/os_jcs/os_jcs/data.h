#pragma once
#include "config.h"
#include "utils.h"

struct MainData
{
	int size;

	double step_t_0;					// Integration step with 0-Hamiltonian
	double step_t_1;					// Integration step with 1-Hamiltonian

	double T;							// Full period
};

void init_main_data(RunParam &rp, ConfigParam &cp, MainData &md);
void delete_main_data(MainData &md);

struct PropData
{
	int total_num_dumps;

	int * dump_periods;

	double * dump_times;

	double * k1;
	double * k2;
	double * k3;
	double * k4;
	double * val;
	double * tmp;
};

void f_basis_prop_std(RunParam &rp, ConfigParam &cp, MainData &md);
void f_basis_prop_deep(RunParam &rp, ConfigParam &cp, MainData &md);

void init_prop_data_std(RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd);
void init_prop_data_deep(RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd);

void free_prop_data_std(RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd);
void free_prop_data_deep(RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd);

void dump_prop_data_std(RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd);
void dump_prop_data_deep(RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd);

