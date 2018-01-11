#pragma once
#include "config.h"
#include "split.h"
#include "utils.h"

struct QJData
{
	VSLStreamStatePtr * streams;				// Random streams
	VSLStreamStatePtr * streams_var;			// Random streams for var trajectories

	MKL_Complex16 * phi_all;					// Phi for all trajectories
	MKL_Complex16 * phi_all_aux;				// Aux phi for all trajectories
	double * abs_diag_rho_all;					// Abs diag rho for all trajectories

	double * times_all;							// Curr times for each trajectory
	double * etas_all;							// Etas
	
	int dump_type;								// Type of dumping
	int num_dumps_total;						// Total number of dumps
	int * dump_periods;							// Dump periods

	// ======== Standart observables ========
	double * mean_start;

	double * mean;
	double * dispersion;
	double * m2;

	double * mean_evo;
	double * dispersion_evo;
	double * m2_evo;
	// ======================================

	// ======== Lyapunov observables ========
	double * delta_s;

	double * energy;
	double * lambda;
	double * mean_lpn;
	double * energy_lpn;
	double * lambda_lpn;

	double * energy_evo;
	double * lambda_evo;
	double * mean_lpn_evo;
	double * energy_lpn_evo;
	double * lambda_lpn_evo;
	// ======================================
};
