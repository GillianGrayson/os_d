#pragma once
#include "config.h"
#include "utils.h"

struct Data
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

	double * data_lpn;		// Data for Lyapunov exps
	double * args_lpn;		// Integration Lyapunov routines data
	double * k1s_lpn;		// Integration Lyapunov routines data
	double * k2s_lpn;		// Integration Lyapunov routines data
	double * k3s_lpn;		// Integration Lyapunov routines data
	double * k4s_lpn;		// Integration Lyapunov routines data
};


void init_main_data(ConfigParam &cp, Data &dt);

void delete_main_data(Data &dt);

void init_cond(ConfigParam &cp, Data &dt);
