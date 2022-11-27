#pragma once
#include "config.h"
#include "split.h"
#include "utils.h"

struct MainData
{
	int sys_size;						// Number of states in system
	int num_diss;						// Number of dissipators
	int num_ham_qj;						// Number of qj hamiltonians

	int mbl_Nc;
	int mbl_Np;
	int mbl_Ns;

	int * mbl_adjacement;	// aux array for define adjacement states
	int * mbl_x_to_id;		// aux array for id definition from binary representation
	int * mbl_id_to_x;      // aux array for binary representation definition from id

	double T;							// Period

	MKL_Complex16* hamiltonian;				// Hamiltonian
	MKL_Complex16* hamiltonian_drv;			// Hamiltonian driving
	MKL_Complex16* hamiltonian_drv_2;			// Hamiltonian driving

	MKL_Complex16 * non_drv_part;		// Non-drv part for rk
	MKL_Complex16 * drv_part;			// Non-drv part for rk
	MKL_Complex16* drv_part_2;			// Non-drv part for rk

	MKL_Complex16 ** dissipators;		// Dissipators
	std::vector<sp_mtx> dissipators_eigen;


	MKL_Complex16 ** hamiltonians_qj;	// Non-Hermitian Hamiltonians for qj

	MKL_Complex16 * special;			// Matrix for special observable
	MKL_Complex16 * special_2;			// Matrix for special observable
	MKL_Complex16 * special_3;			// Matrix for special observable
	MKL_Complex16 * special_4;			// Matrix for special observable

	std::vector<MKL_Complex16*> random_obs_mtxs;

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
	int * jumps_counts;							// Jump count for each trajectory

	int is_obs;

	vector<double> * jump_times;				// Times of jumps
	vector<int> * diss_types;
	vector<double> * jump_norms;				// Norms before jumps
	vector<double> * jump_etas;					// Etas

	double curr_time;								// Curr period_id

	int dump_type;								// Type of dumping
	int dump_num_total;							// Total number of dumps
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
	MKL_Complex16** drv_tmp_2;
	// =========================

	// ======== Standart observables ========
	double * mean_start;

	double * norm;
	double * mean;
	double * dispersion;
	double * m2;
	MKL_Complex16 * spec;
	MKL_Complex16 * spec_2;
	MKL_Complex16 * spec_3;
	std::vector<std::vector<std::complex<double>>> random_obs; // 1 - trajectory, 2 - observable

	double * norm_evo;
	double * mean_evo;
	double * dispersion_evo;
	double * m2_evo;
	MKL_Complex16 * spec_evo;
	MKL_Complex16 * spec_2_evo;
	MKL_Complex16 * spec_3_evo;
	std::vector<std::vector<std::vector<std::complex<double>>>> random_obs_evo; // 1 - trajectory, 2 - observable, 3 - time
	// ======================================

	// ======== Lyapunov observables ========
	double max_energy;

	double * delta_s;

	double * energy;
	double * lambda;
	double * lambda_now;
	double * mean_lpn;
	double * energy_lpn;
	MKL_Complex16 * spec_lpn;
	MKL_Complex16 * spec_2_lpn;
	MKL_Complex16 * spec_3_lpn;
	std::vector<std::vector<std::complex<double>>> random_obs_lpn; // 1 - trajectory, 2 - observable

	double * energy_evo;
	double * lambda_evo;
	double * mean_lpn_evo;
	double * energy_lpn_evo;
	MKL_Complex16 * spec_lpn_evo;
	MKL_Complex16 * spec_2_lpn_evo;
	MKL_Complex16 * spec_3_lpn_evo;
	std::vector<std::vector<std::vector<std::complex<double>>>> random_obs_lpn_evo; // 1 - trajectory, 2 - observable, 3 - time

	int * num_renorms;
	vector<double>* lambdas;
	vector<double>* deltas_s;
	vector<double>* deltas_f;
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

