#pragma once
#include "config.h"
#include "split.h"
#include "utils.h"

struct QJData
{
	VSLStreamStatePtr * streams;				// Random streams
	VSLStreamStatePtr * streams_var;			// Random streams for var trajectories

	MKL_Complex16 * phi_all;					// Phi for all trajectories
	double * abs_diag_rho_all;					// Abs diag rho for all trajectories

	double * etas;								// Etas

	int * dump_periods;							// Dump periods

	double * mean_start;						
	double * mean;
	double * dispersion;
	double * m2;
	double * energy;
	double * lambda;
	double * delta_s;
};
