#pragma once
#include "config.h"
#include "split.h"
#include "utils.h"

struct MainData
{
	int sys_size;						// Number of states in system
	int num_diss;						// Number of dissipators
	int num_ham_qj;						// Number of qj hamiltonians

	double T;							// Period

	double * hamiltonian;				// Hamiltonian
	double * hamiltonian_drv;			// Hamiltonian driving

	MKL_Complex16 * non_drv_part;		// Non-drv part for rk
	MKL_Complex16 * drv_part;		// Non-drv part for rk

	MKL_Complex16 ** dissipators;		// Dissipators

	MKL_Complex16 ** hamiltonians_qj;	// Non-Hermitian Hamiltonians for qj

	Split * structure;					// Split structure
	Split * splits;						// Splits
};

struct ExpData
{
	VSLStreamStatePtr * streams;				// Random streams
	VSLStreamStatePtr * streams_var;			// Random streams for var trajectories

	MKL_Complex16 * phi_all;					// Phi for all trajectories
	MKL_Complex16 * phi_all_aux;				// Aux phi for all trajectories
	double * abs_diag_rho_all;					// Abs diag rho for all trajectories

	double * times_all;							// Curr times for each trajectory
	double * etas_all;							// Etas

	vector<double> * jump_times;				// Times of jumps
	vector<pair<double, double>> * jump_norms;	// Norms before and after during jumps

	int period_id;								// Curr period_id
	int dump_type;								// Type of dumping
	int num_dumps_total;						// Total number of dumps
	int * dump_periods;							// Dump periods

	// ======== RK Data ========
	double rk_step;

	MKL_Complex16 ** k1;
	MKL_Complex16 ** k2;
	MKL_Complex16 ** k3;
	MKL_Complex16 ** k4;

	MKL_Complex16 ** args;

	MKL_Complex16 ** non_drv_tmp;
	MKL_Complex16 ** drv_tmp;
	// =========================


	// ======== Standart observables ========
	double * mean_start;

	double * norm;
	double * mean;
	double * dispersion;
	double * m2;

	double * norm_evo;
	double * mean_evo;
	double * dispersion_evo;
	double * m2_evo;
	// ======================================

	// ======== Lyapunov observables ========
	double max_energy;

	double * delta_s;

	double * energy;
	double * lambda;
	double * lambda_now;
	double * mean_lpn;
	double * energy_lpn;

	double * energy_evo;
	double * lambda_evo;
	double * mean_lpn_evo;
	double * energy_lpn_evo;
	// ======================================

	// ===== Corr dimension observables =====
	int cd_shift_size;
	int cd_dim;
	int cd_num_points;
	double *** cd_rec_data;
	double * cd_i;
	// ======================================
};

struct AllData
{
	RunParam * rp;
	ConfigParam * cp;
	MainData * md;
	ExpData * ed;
};

