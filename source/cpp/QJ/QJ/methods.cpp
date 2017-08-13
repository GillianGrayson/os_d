#include "header.h"
#include "output.h"
#include "omp.h"
#include <iostream>
#include <mkl.h>
#include <complex>
#include <vector>
#include <math.h>
#include <algorithm>

#define EPS_ETA 1.0e-10

using namespace std;

#define IND(a,b)  (a) * N + (b)
#define IND2(a,b) (a) * N * N + (b)


MKL_Complex16 ZERO={0,0}, ONE={1,0}, I={0,1};
double * ts;

MKL_Complex16 * phi_global;
double * etas;
MKL_Complex16 * rho_omp;
double * abs_rho_diag_omp;
double * abs_rho_diag_all;
double * avg_rho_diag_all;

double * abs_rho_diag_all_stationary;
MKL_Complex16 * phi_normed;
MKL_Complex16 * phi_stationary;

MKL_Complex16 * rho_curr;
MKL_Complex16 * rho_curr_and;
MKL_Complex16 * hamiltonian;
MKL_Complex16 * eg_hamiltonian;
MKL_Complex16 * ev_hamiltonian;
MKL_Complex16 * ev_t_hamiltonian;

MKL_Complex16 * ev_hamiltonian_sorted;
MKL_Complex16 * ev_t_hamiltonian_sorted;

int * trans_process_end_period;
double * mean_start;
double * mean_start_stationary;

double * mean;
double * dispersion;
double * m2;

int num_att_peaks;
double * begins_att_peaks = NULL;
double * ends_att_peaks = NULL;
bool is_special_peak;
vector<int> ** sticking_times;
int * current_peaks_ids;

double * mean_stationary;
double * dispersion_stationary;
double * m2_stationary;
int * max_id_stationary;

int num_att_peaks_stationary;
double * begins_att_peaks_stationary = NULL;
double * ends_att_peaks_stationary = NULL;
bool is_special_peak_stationary;
vector<int> ** sticking_times_stationary;
int * current_peaks_ids_stationary;

vector<double> * jump_times;

int * histogramm;
double * energy_intervals;
double * mass_centers_intervals;
double e_min;
double e_max;
double e_shift;
double mc_min;
double mc_max;
double mc_shift;
MKL_Complex16 * sub_mult;

int period;

int periodic_index(
	int index,
	int max_index,
	int N
)
{
	int left = max_index - N / 2;
	int right = max_index + (N - N / 2);

	if (max_index < N / 2)
	{
		if (index >= right)
		{
			index -= N;
			return index;
		}
		else
		{
			return index;
		}
	}
	else
	{
		if (index < left)
		{
			index += N;
			return index;
		}
		else
		{
			return index;
		}
	}
}

double get_mean(
	int periodic,
	int N,
	int max_index,
	double * adr,
	int trajectory_id, 
	int mc_specific
)
{
	double current_mean = 0.0;
	double current_mean_on_period = 0.0;
	double previous_mean_on_period = 0.0;
	double distance_1 = 0.0;
	double distance_2 = 0.0;

	if (periodic == 0)
	{
		for (int state_id = 0; state_id < N; state_id++)
		{
			current_mean += double(state_id) * adr[state_id];
		}
	}
	else
	{
		if (mc_specific > 0)
		{
			std::vector<double> adr_vec(adr, adr + N);

			std::vector<int> order(adr_vec.size());

			std::size_t n(0);
			std::generate(std::begin(order), std::end(order), [&] { return n++; });

			std::sort(std::begin(order), std::end(order),
				[&](int i1, int i2) { return adr_vec[i1] > adr_vec[i2]; });

			int max_index_id = 0;
			double min_mc_dist = double(N);
			for (int index_id = 0; index_id < mc_specific; index_id++)
			{
				if (abs(order[index_id] - mean[trajectory_id]) > min_mc_dist)
				{
					min_mc_dist = abs(order[index_id] - mean[trajectory_id]);
					max_index_id = index_id;
				}
			}

			max_index = order[max_index_id];
		}

		for (int state_id = 0; state_id < N; state_id++)
		{
			current_mean += double(periodic_index(state_id, max_index, N)) * adr[state_id];
		}

		current_mean_on_period = fmod(current_mean, double(N));
		if (current_mean_on_period < 0.0)
		{
			current_mean_on_period += double(N);
		}

		previous_mean_on_period = fmod(mean[trajectory_id], double(N));
		if (previous_mean_on_period < 0.0)
		{
			previous_mean_on_period += double(N);
		}

		distance_1 = current_mean_on_period - previous_mean_on_period;
		if (distance_1 > 0)
		{
			distance_2 = -(double(N) - fabs(distance_1));
		}
		else
		{
			distance_2 = (double(N) - fabs(distance_1));
		}

		if (fabs(distance_1) < fabs(distance_2))
		{
			current_mean = mean[trajectory_id] + distance_1;
		}
		else
		{
			current_mean = mean[trajectory_id] + distance_2;
		}
	}

	return current_mean;
}

double get_mean_stationary(
	int periodic,
	int N,
	int max_index,
	double * adr,
	int trajectory_id,
	int mc_specific
)
{
	double current_mean = 0.0;
	double current_mean_on_period = 0.0;
	double previous_mean_on_period = 0.0;
	double distance_1 = 0.0;
	double distance_2 = 0.0;

	if (periodic == 0)
	{
		for (int state_id = 0; state_id < N; state_id++)
		{
			current_mean += double(state_id) * adr[state_id];
		}
	}
	else
	{
		current_mean = double(max_index);
		current_mean_on_period = current_mean;

		previous_mean_on_period = fmod(mean_stationary[trajectory_id], double(N));
		if (previous_mean_on_period < 0.0)
		{
			previous_mean_on_period += double(N);
		}

		distance_1 = current_mean_on_period - previous_mean_on_period;
		if (distance_1 > 0)
		{
			distance_2 = -(double(N) - fabs(distance_1));
		}
		else
		{
			distance_2 = (double(N) - fabs(distance_1));
		}

		if (fabs(distance_1) < fabs(distance_2))
		{
			current_mean = mean_stationary[trajectory_id] + distance_1;
		}
		else
		{
			current_mean = mean_stationary[trajectory_id] + distance_2;
		}
	}

	return current_mean;
}

double get_init_mean(
	int periodic,
	int N,
	int max_index,
	double * adr
)
{
	double current_mean = 0.0;

	int index = 0;
	for (int state_id = 0; state_id < N; state_id++)
	{
		if (periodic == 0)
		{
			index = state_id;
		}
		else
		{
			index = periodic_index(state_id, max_index, N);
		}

		current_mean += double(index) * adr[state_id];
	}
	

	return current_mean;
}

double get_init_mean_stationary(
	int periodic,
	int N,
	int max_index,
	double * adr
)
{
	double current_mean = double(max_index);

	return current_mean;
}

double get_mean_on_period(
	int N,
	double curr_mean
)
{
	double current_mean_on_period = 0.0;

	current_mean_on_period = fmod(curr_mean, double(N));
	if (current_mean_on_period < 0.0)
	{
		current_mean_on_period += double(N);
	}

	return current_mean_on_period;
}

double get_dispersion(
	double current_mean,
	int trajectory_id
)
{
	double curr_dispersion = 0.0;

	curr_dispersion = (current_mean - mean_start[trajectory_id]) * (current_mean - mean_start[trajectory_id]);

	return curr_dispersion;
}

double get_dispersion_stationary(
	double current_mean,
	int trajectory_id
)
{
	double curr_dispersion = 0.0;

	curr_dispersion = (current_mean - mean_start_stationary[trajectory_id]) * (current_mean - mean_start_stationary[trajectory_id]);

	return curr_dispersion;
}

double get_m2(
	int periodic,
	int N,
	double curr_mean_on_period,
	double * adr
)
{
	double curr_m2 = 0.0;

	int index = 0;
	for (int state_id = 0; state_id < N; state_id++)
	{
		if (periodic == 0)
		{
			index = state_id;
		}
		else
		{
			index = periodic_index(state_id, int(curr_mean_on_period), N);
		}

		curr_m2 += (double(index) - curr_mean_on_period) * (double(index) - curr_mean_on_period) * adr[state_id];
	}

	return curr_m2;
}

void st_processing_direct(
	int periodic,
	int N,
	double curr_mean_on_period,
	int trajectory_id
)
{

	double mean_in_period = curr_mean_on_period;

	if (periodic == 1)
	{
		if (num_att_peaks > 0)
		{
			if (mean_in_period + double(N) < ends_att_peaks[num_att_peaks - 1])
			{
				mean_in_period += double(N);
			}
		}
	}

	if (num_att_peaks > 0)
	{
		// searching attractor peak
		int attractor_peak_id = -1; // not belongs to attractor peak

		for (int peak_id = 0; peak_id < num_att_peaks; peak_id++)
		{
			if ((mean_in_period > begins_att_peaks[peak_id]) && (mean_in_period < ends_att_peaks[peak_id]))
			{
				attractor_peak_id = peak_id;
				break;
			}
		}

		if (attractor_peak_id == current_peaks_ids[trajectory_id])
		{
			if (attractor_peak_id >= 0)
			{
				sticking_times[trajectory_id][attractor_peak_id][sticking_times[trajectory_id][attractor_peak_id].size() - 1] ++;
			}
		}
		else
		{
			if (attractor_peak_id >= 0)
			{
				sticking_times[trajectory_id][attractor_peak_id].push_back(0);
			}

			current_peaks_ids[trajectory_id] = attractor_peak_id;
		}
	}
}

void st_processing_stationary(
	int periodic,
	int N,
	int max_index,
	int trajectory_id
)
{
	if (num_att_peaks_stationary > 0)
	{
		// searching attractor peak
		int attractor_peak_id_stationary = -1; // not belongs to attractor peak

		for (int peak_id = 0; peak_id < num_att_peaks_stationary; peak_id++)
		{
			if ((max_index > begins_att_peaks_stationary[peak_id]) && (max_index < ends_att_peaks_stationary[peak_id]))
			{
				attractor_peak_id_stationary = peak_id;
				break;
			}
		}

		if (attractor_peak_id_stationary == current_peaks_ids_stationary[trajectory_id])
		{
			if (attractor_peak_id_stationary >= 0)
			{
				sticking_times_stationary[trajectory_id][attractor_peak_id_stationary][sticking_times_stationary[trajectory_id][attractor_peak_id_stationary].size() - 1] ++;
			}
		}
		else
		{
			if (attractor_peak_id_stationary >= 0)
			{
				sticking_times_stationary[trajectory_id][attractor_peak_id_stationary].push_back(0);
			}

			current_peaks_ids_stationary[trajectory_id] = attractor_peak_id_stationary;
		}
	}
}

inline MKL_Complex16 Complex_mul(MKL_Complex16 a, double b)
{
	MKL_Complex16 res = {0.0, 0.0};
	res.real = a.real * b;
	res.imag = a.imag * b;
	return res;
}

inline MKL_Complex16 Complex_scalar_mul(MKL_Complex16 * a, MKL_Complex16 * b, int N)
{
	MKL_Complex16 res = {0.0, 0.0};
	for (int i = 0; i < N; i++)
	{
		res.real += a[i].real * b[i].real + a[i].imag * b[i].imag;
		res.imag += b[i].real * a[i].imag - a[i].real * b[i].imag;
	}
	return res;
}

int if_norm(MKL_Complex16 * phi, double * eta, int N)
{
	double norm=0.;
	for(int i=0; i<N; i++)
	{
		norm+=phi[i].real*phi[i].real+phi[i].imag*phi[i].imag;
	}
	if((norm)>eta[0])
		return 0;
	else
	{
		return 1;
	}
}

double norm_vector2(MKL_Complex16 * phi, int N)
{
	double norm=0.;
	for(int i=0; i<N; i++)
	{
		norm+=phi[i].real*phi[i].real+phi[i].imag*phi[i].imag;
	}
	return norm;
}

void print_complex_array(MKL_Complex16 * tmp, int N)
{
	printf("\n");
	for (int i = 0; i < N; i++)
	{	
		printf("%d: %0.2le %0.2le\n", i, tmp[i].real, tmp[i].imag);
	}
}

void print_double_array(double * tmp, int N)
{
	printf("\n");
	for (int i = 0; i < N; i++)
	{	
		printf("%d: %0.2le\n", i, tmp[i]);
	}
}


void init_data(int N,
			   int num_trajectories,
			   int num_omp_threads)
{
	phi_global = new MKL_Complex16[num_trajectories * N];
	abs_rho_diag_omp = new double[num_omp_threads * N];
	rho_omp = new MKL_Complex16[num_omp_threads * N * N];

	ts = new double[num_trajectories];

	etas = new double[num_trajectories];
	period = 0;

	for (int trajectoty_id = 0; trajectoty_id < num_trajectories; trajectoty_id++)
	{
		for (int state_id = 0; state_id < N; state_id++)
		{
			phi_global[trajectoty_id * N + state_id].real = 0.0;
			phi_global[trajectoty_id * N + state_id].imag = 0.0;
		}

		ts[trajectoty_id] = 0.0;

		etas[trajectoty_id] = 0.0;
	}

	for (int state_id_1 = 0; state_id_1 < N; state_id_1++)
	{
		for (int thread_id = 0; thread_id < num_omp_threads; thread_id ++)
		{
			abs_rho_diag_omp[thread_id * (N) + state_id_1] = 0.0;

			for (int state_id_2 = 0; state_id_2 < N; state_id_2++)
			{
				rho_omp[thread_id * (N * N) + state_id_1 * (N) + state_id_2].real = 0.0;
				rho_omp[thread_id * (N * N) + state_id_1 * (N) + state_id_2].imag = 0.0;
			}
		}
	}
}

void init_dump_statistics_data(int num_trajectories, int btw_jump_times, int stationary)
{
	trans_process_end_period = new int[num_trajectories];
	mean_start = new double[num_trajectories];

	mean = new double[num_trajectories];
	dispersion = new double[num_trajectories];
	m2 = new double[num_trajectories];

	if (stationary == 1)
	{
		mean_start_stationary = new double[num_trajectories];

		mean_stationary = new double[num_trajectories];
		dispersion_stationary = new double[num_trajectories];
		m2_stationary = new double[num_trajectories];
		max_id_stationary = new int[num_trajectories];
	}

	if (btw_jump_times > 0)
	{
		jump_times = new vector<double>[num_trajectories];
	}

	for (int trajectory_id = 0; trajectory_id < num_trajectories; trajectory_id++)
	{
		trans_process_end_period[trajectory_id] = 0;
		mean_start[trajectory_id] = 0.0;

		mean[trajectory_id] = 0.0;
		dispersion[trajectory_id] = 0.0;
		m2[trajectory_id] = 0.0;

		if (stationary == 1)
		{
			mean_start_stationary[trajectory_id] = 0.0;

			mean_stationary[trajectory_id] = 0.0;
			dispersion_stationary[trajectory_id] = 0.0;
			m2_stationary[trajectory_id] = 0.0;
			max_id_stationary[trajectory_id] = 0;
		}
	}
}

void init_dump_hist_data(int N, int num_trajectories, int num_e_intervals, int num_mc_intervals, double energy_min, double energy_max)
{
	energy_intervals = new double[num_e_intervals];
	mass_centers_intervals = new double[num_mc_intervals];
	histogramm = new int[num_trajectories * num_e_intervals * num_mc_intervals];

	for (int trajectory_id = 0; trajectory_id < num_trajectories; trajectory_id++)
	{
		for (int e_interval_id = 0; e_interval_id < num_e_intervals; e_interval_id ++)
		{
			for (int mc_interval_id = 0; mc_interval_id < num_mc_intervals; mc_interval_id ++)
			{
				histogramm[trajectory_id * (num_e_intervals * num_mc_intervals) + e_interval_id * num_mc_intervals + mc_interval_id] = 0;
			}
		}
	}

	e_min = energy_min;
	e_max = energy_max;
	e_shift = (e_max - e_min) / double(num_e_intervals);
	for (int e_interval_id = 0; e_interval_id < num_e_intervals; e_interval_id ++)
	{
		energy_intervals[e_interval_id] = e_min + double(e_interval_id) * e_shift + 0.5 * e_shift;
	}

	mc_min = 0.5;
	mc_max = double(N) + 0.5;
	mc_shift = (mc_max - mc_min) / double(num_mc_intervals);
	for (int mc_interval_id = 0; mc_interval_id < num_mc_intervals; mc_interval_id ++)
	{
		mass_centers_intervals[mc_interval_id] = mc_min + double(mc_interval_id) * mc_shift + 0.5 * mc_shift;
	}
}

void delete_dump_hist_data()
{
	delete[] energy_intervals;
	delete[] mass_centers_intervals;
	delete[] histogramm;
}

void refresh_dump_statistics_data(int num_trajectories)
{
	for (int trajectory_id = 0; trajectory_id < num_trajectories; trajectory_id++)
	{
		mean[trajectory_id] = 0.0;
		dispersion[trajectory_id] = 0.0;
		m2[trajectory_id] = 0.0;
	}
}

void delete_dump_statistics_data(int dump_characteristics, int btw_jump_times, int stationary)
{
	delete[] trans_process_end_period;
	delete[] mean_start;

	delete[] mean;
	delete[] dispersion;
	delete[] m2;

	if (stationary == 1)
	{
		delete[] mean_start_stationary;

		delete[] mean_stationary;
		delete[] dispersion_stationary;
		delete[] m2_stationary;
		delete[] max_id_stationary;
	}

	if (btw_jump_times > 0)
	{
		delete[] jump_times;
	}
}

void init_main_statistics_data(int N,
							   int num_trajectories,
							   int stationary)
{
	phi_global = new MKL_Complex16[num_trajectories * N];
	abs_rho_diag_all = new double[num_trajectories * N];
	avg_rho_diag_all = new double[num_trajectories * N];

	ts = new double[num_trajectories];

	if(stationary == 1)
	{
		abs_rho_diag_all_stationary = new double[num_trajectories * N];
		phi_normed = new MKL_Complex16[num_trajectories * N];
		phi_stationary = new MKL_Complex16[num_trajectories * N];
	}

	etas = new double[num_trajectories];
	period = 0;

	for (int trajectory_id = 0; trajectory_id < num_trajectories; trajectory_id++)
	{
		for (int state_id = 0; state_id < N; state_id++)
		{
			phi_global[trajectory_id * N + state_id].real = 0.0;
			phi_global[trajectory_id * N + state_id].imag = 0.0;

			abs_rho_diag_all[trajectory_id * N + state_id] = 0.0;
			avg_rho_diag_all[trajectory_id * N + state_id] = 0.0;

			if (stationary == 1)
			{
				abs_rho_diag_all_stationary[trajectory_id * N + state_id] = 0.0;

				phi_normed[trajectory_id * N + state_id].real = 0.0;
				phi_normed[trajectory_id * N + state_id].imag = 0.0;

				phi_stationary[trajectory_id * N + state_id].real = 0.0;
				phi_stationary[trajectory_id * N + state_id].imag = 0.0;
			}
		}

		etas[trajectory_id] = 0.0;

		ts[trajectory_id] = 0.0;
	}
}

void init_main_hist_data(int N,
						 int num_trajectories)
{
	phi_global = new MKL_Complex16[num_trajectories * N];
	sub_mult = new MKL_Complex16[num_trajectories * N];
	abs_rho_diag_all = new double[num_trajectories * N];

	ts = new double[num_trajectories];

	etas = new double[num_trajectories];
	period = 0;

	for (int trajectory_id = 0; trajectory_id < num_trajectories; trajectory_id++)
	{
		for (int state_id = 0; state_id < N; state_id++)
		{
			phi_global[trajectory_id * N + state_id].real = 0.0;
			phi_global[trajectory_id * N + state_id].imag = 0.0;

			sub_mult[trajectory_id * N + state_id].real = 0.0;
			sub_mult[trajectory_id * N + state_id].imag = 0.0;

			abs_rho_diag_all[trajectory_id * N + state_id] = 0.0;
		}

		ts[trajectory_id] = 0.0;

		etas[trajectory_id] = 0.0;
	}
}

void refresh_main_statistics_data(int N,
								  int num_trajectories,
								  int dump_characteristics,
								  int stationary)
{
	for (int trajectoty_id = 0; trajectoty_id < num_trajectories; trajectoty_id++)
	{
		for (int state_id = 0; state_id < N; state_id++)
		{
			abs_rho_diag_all[trajectoty_id * N + state_id] = 0.0;

			if (stationary == 1)
			{
				abs_rho_diag_all_stationary[trajectoty_id * N + state_id] = 0.0;
			}
		}
	}
}

void refresh_main_hist_data(int N,
							int num_trajectories)
{
	for (int trajectoty_id = 0; trajectoty_id < num_trajectories; trajectoty_id++)
	{
		for (int state_id = 0; state_id < N; state_id++)
		{
			abs_rho_diag_all[trajectoty_id * N + state_id] = 0.0;
		}
	}
}

void delete_main_statistics_data(int dump_characteristics, int stationary)
{
	delete[] phi_global;
	delete[] abs_rho_diag_all;
	delete[] avg_rho_diag_all;

	delete[] ts;

	if (stationary == 1)
	{
		delete[] abs_rho_diag_all_stationary;
		delete[] phi_normed;
		delete[] phi_stationary;
	}
	delete[] etas;
}

void delete_main_hist_data()
{
	delete[] phi_global;
	delete[] abs_rho_diag_all;
	delete[] sub_mult;

	delete[] ts;

	delete[] etas;
}


void refresh_dump_data(int N,
					   int num_omp_threads)
{
	for (int state_id_1 = 0; state_id_1 < N; state_id_1++)
	{
		for (int thread_id = 0; thread_id < num_omp_threads; thread_id ++)
		{
			abs_rho_diag_omp[thread_id * (N) + state_id_1] = 0.0;

			for (int state_id_2 = 0; state_id_2 < N; state_id_2++)
			{
				rho_omp[thread_id * (N * N) + state_id_1 * (N) + state_id_2].real = 0.0;
				rho_omp[thread_id * (N * N) + state_id_1 * (N) + state_id_2].imag = 0.0;
			}
		}
	}
}

inline void recovery_phi_full(MKL_Complex16 * phi, double eta, double * g, void * A, int N, unsigned int k, VSLStreamStatePtr streamRand)  // MKL_complex16 * A
{
	int index = 0;
	MKL_Complex16 * res = new MKL_Complex16[N];
	double norm = sqrt(norm_vector2(phi, N));
	for(int i=0; i<N; i++)
	{
		phi[i].real /= (norm);
		phi[i].imag /= (norm);
	}
	norm=0.;

	double * gnorms = new double[k];
	double tmp = 0;
	double ran;
	vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, streamRand, 1, &ran, 0, 1);

	for(unsigned int i = 0; i < k; i++)
	{
		cblas_zgemv (CblasRowMajor, CblasNoTrans, N, N, &ONE, &((MKL_Complex16*)A)[IND2(i, 0)], N, phi, 1, &ZERO, res, 1);
		gnorms[i]=(norm_vector2(res, N));
		gnorms[i] *= g[i];
		tmp += gnorms[i];
	}

	ran *= tmp;

	while(ran - gnorms[index] > 0.)
	{
		ran -= gnorms[index];
		index++;
		if(index == k - 1)
			break;
	}

	while(gnorms[index] == 0)
	{
		if (index == 0)
		{
			index++;
		}
		else
		{
			index--;
		}
	}

	memset(res, 0, N*sizeof(MKL_Complex16));
	cblas_zgemv (CblasRowMajor, CblasNoTrans,N, N, &ONE, &(((MKL_Complex16*)A)[IND2(index,0)]), N, phi, 1, &ZERO, res, 1);

	for(int i=0; i<N; i++)
	{
		phi[i].real=res[i].real/sqrt(gnorms[index]/g[index]);
		phi[i].imag=res[i].imag/sqrt(gnorms[index]/g[index]);
	}

	delete (res);
	delete (gnorms);
}

inline void QJ_step(MKL_Complex16 * phi, MKL_Complex16 * matrix, MKL_Complex16 * res, unsigned int N)
{
	cblas_zgemv (CblasRowMajor, CblasNoTrans, N, N, &ONE, matrix, N, phi, 1, &ZERO, res, 1);
}

void QJ_EXP_one_branch(
	MKL_Complex16 * phi,
	MKL_Complex16 * tmp_vec,
	double * eta,
	double * g,
	void * A,
	unsigned int k,
	split * branch,
	VSLStreamStatePtr streamRand,
	int N,
	int trajectory_id,
	int btw_jump_times,
	int deep_characteristics,
	int mc_specific,
	int periodic,
	int stationary
)
{
	if (branch->next == 0)
	{
		while (branch->counter != branch->steps)
		{
			QJ_step(phi, branch->matrix, tmp_vec, branch->N);
			if(if_norm(tmp_vec, eta, branch->N))
			{
				recovery_phi_full(tmp_vec, eta[0], g, A, branch->N, k, streamRand);
				vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, streamRand, 1, eta, 0, 1);
				while (eta[0] == 0.0)
				{
					vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, streamRand, 1, eta, 0, 1);
				}

				if (btw_jump_times > 0)
				{
					jump_times[trajectory_id].push_back(ts[trajectory_id] + branch->dt);
				}

				if (deep_characteristics > 0)
				{
					double * adr = &abs_rho_diag_all[trajectory_id * N];
					double * adr_stationary = &abs_rho_diag_all_stationary[trajectory_id * N];

					double curr_norm = norm_vector2(tmp_vec, N);

					MKL_Complex16 * phi_n = NULL;
					MKL_Complex16 * phi_st = NULL;
					if (stationary == 1)
					{
						phi_n = &phi_normed[trajectory_id * N];
						phi_st = &phi_stationary[trajectory_id * N];
					}

					double max_val_direct = 0.0;
					double max_index_direct = 0;
					for (int state_id = 0; state_id < N; state_id++)
					{
						adr[state_id] = Complex_mul(Complex_scalar_mul(&tmp_vec[state_id], &tmp_vec[state_id], 1), 1.0 / curr_norm).real;

						if (adr[state_id] > max_val_direct)
						{
							max_val_direct = adr[state_id];
							max_index_direct = state_id;
						}

						if (stationary == 1)
						{
							phi_n[state_id].real = tmp_vec[state_id].real / sqrt(curr_norm);
							phi_n[state_id].imag = tmp_vec[state_id].imag / sqrt(curr_norm);

							phi_st[state_id].real = 0.0;
							phi_st[state_id].imag = 0.0;
						}
					}

					double curr_mean = get_mean(periodic, N, max_index_direct, adr, trajectory_id, mc_specific);
					mean[trajectory_id] = curr_mean;

					if (stationary == 1)
					{
						cblas_zgemv(CblasRowMajor, CblasNoTrans, N, N, &ONE, ev_t_hamiltonian_sorted, N, phi_n, 1, &ZERO, phi_st, 1);
						double max_val_stationary = 0.0;
						int max_index_stationary = 0;
						for (int state_id = 0; state_id < N; state_id++)
						{
							adr_stationary[state_id] = Complex_mul(Complex_scalar_mul(&phi_st[state_id], &phi_st[state_id], 1), 1.0).real;

							if (adr_stationary[state_id] > max_val_stationary)
							{
								max_val_stationary = adr_stationary[state_id];
								max_index_stationary = state_id;
							}
						}

						double curr_mean_stationary = get_mean_stationary(periodic, N, max_index_direct, adr_stationary, trajectory_id, mc_specific);
						mean_stationary[trajectory_id] = curr_mean_stationary;
					}
				}
			}

			memcpy(phi, tmp_vec, sizeof(MKL_Complex16)*branch->N);
			branch->counter++;

			ts[trajectory_id] = ts[trajectory_id] + branch->dt;
		}
	}
	else
	{
		while (branch->counter != branch->steps)
		{
			QJ_step(phi, branch->matrix, tmp_vec, branch->N);
			if(if_norm(tmp_vec, eta, branch->N))
			{
				QJ_EXP_one_branch(phi, tmp_vec, eta, g, A, k, branch->next, streamRand, N, trajectory_id, btw_jump_times, deep_characteristics, mc_specific, periodic, stationary);
				ts[trajectory_id] = ts[trajectory_id] - branch->dt;
			}
			else
			{
				memcpy(phi, tmp_vec, sizeof(MKL_Complex16)*branch->N);
			}
			branch->counter++;
			ts[trajectory_id] = ts[trajectory_id] + branch->dt;
		}
	}

	branch->counter = 0;
}

void qj_propagate_one_period(
	MKL_Complex16 * phi,
	double * eta,
	split * head,
	VSLStreamStatePtr streamRand,
	int N,
	int trajectory_id,
	int btw_jump_times,
	int deep_characteristics,
	int mc_specific,
	int periodic,
	int stationary
)
{
	MKL_Complex16 * tmp_vec = new MKL_Complex16[head->N];
	for (unsigned int i = 0; i < head->counter; i++)
	{
		QJ_EXP_one_branch(phi, tmp_vec, eta, head->g, head->matrix, head->steps, &(head->next)[i], streamRand, N, trajectory_id, btw_jump_times, deep_characteristics, mc_specific, periodic, stationary);
	}
	delete(tmp_vec);
}

void set_init_conditions(MKL_Complex16 * phi, int N, int state_id)
{
	if (state_id >= 0 && state_id < N)
	{
		phi[state_id].real = 1.0;
	}
	else
	{
		for (int i = 0; i < N; i++)
		{
			phi[i].real = sqrt(1.0/double(N));
		}
	}
}

void calc_and_dump_characteristics(int N, char* write_type)
{
	char characteristics_file_name[] = "characteristics.txt";

	double start_time = omp_get_wtime();

	FILE * characteristics_file = fopen(characteristics_file_name, write_type);
	fprintf(characteristics_file, "%0.18le\n", 0.0);
	fclose(characteristics_file);

	double time = omp_get_wtime() - start_time;
	printf("routines: %0.3le\n", time);
}

void single_trajectory_init_prop(split * head,
											   VSLStreamStatePtr rnd_stream,
											   MKL_Complex16 * phi,
											   double * abs_rho_diag,
											   MKL_Complex16 * rho,
											   int num_periods_in_trans_proc,
											   int num_trajectories,
											   int trajectory_id,
											   int init_state_id,
											   int thread_id)
{
	int N = head->N;

	set_init_conditions(phi, N, init_state_id);

	double eta = 0.0;
	while (eta == 0.0)
	{
		vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, rnd_stream, 1, &eta, 0, 1);
	}

	if (num_periods_in_trans_proc > 0)
	{
		for(int period_id = 0; period_id < num_periods_in_trans_proc; period_id++)
		{
			qj_propagate_one_period(phi, &eta, head, rnd_stream, N, trajectory_id, 0, 0, 0, 0, 0);
		}
	}
	else
	{
		qj_propagate_one_period(phi, &eta, head, rnd_stream, N, trajectory_id, 0, 0, 0, 0, 0);
	}

	etas[trajectory_id] = eta;

	double norm = norm_vector2(phi, N);
	for(int state_id_1 = 0; state_id_1 < N; state_id_1++)
	{
		abs_rho_diag[state_id_1] += Complex_mul(Complex_scalar_mul(&phi[state_id_1], &phi[state_id_1], 1), 1.0 / norm).real / num_trajectories;

		for (int state_id_2 = 0; state_id_2 < N; state_id_2++)
		{
			rho[state_id_1 * (N) + state_id_2].real += Complex_mul(Complex_scalar_mul(&phi[state_id_1], &phi[state_id_2], 1), 1.0 / norm).real / num_trajectories;
			rho[state_id_1 * (N) + state_id_2].imag += Complex_mul(Complex_scalar_mul(&phi[state_id_1], &phi[state_id_2], 1), 1.0 / norm).imag / num_trajectories;
		}
	}
}

void single_trajectory_statistics_init_prop(
	split * head,
	VSLStreamStatePtr rnd_stream,
	MKL_Complex16 * phi,
	double * adr,
	int num_periods_in_trans_proc,
	int trajectory_id,
	int init_state_id,
	double mean_low_limit,
	double mean_high_limit,
	int dump_characteristics,
	double * adr_stationary,
	int stationary,
	int periodic
)
{
	int N = head->N;

	set_init_conditions(phi, N, init_state_id);

	double curr_norm = 0.0;

	double max_val_direct = 0.0;
	int max_index_direct = 0;

	double max_val_stationary = 0.0;
	int max_index_stationary = 0;

	double curr_mean = 0.0;
	double curr_mean_stationary = 0.0;

	MKL_Complex16 * phi_n = NULL;
	MKL_Complex16 * phi_st = NULL;
	if (stationary == 1)
	{
		phi_n = &phi_normed[trajectory_id * N];
		phi_st = &phi_stationary[trajectory_id * N];
	}

	double eta = 0.0;
	while (eta == 0.0)
	{
		vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, rnd_stream, 1, &eta, 0, 1);
	}


	for(int period_id = 0; period_id < num_periods_in_trans_proc; period_id++)
	{
		qj_propagate_one_period(phi, &eta, head, rnd_stream, N, trajectory_id, 0, 0, 0, periodic, stationary);
	}

	trans_process_end_period[trajectory_id] = num_periods_in_trans_proc;

	curr_norm = norm_vector2(phi, N);

	max_val_direct = 0.0;
	max_index_direct = 0;
	for(int state_id = 0; state_id < N; state_id++)
	{
		adr[state_id] = Complex_mul(Complex_scalar_mul(&phi[state_id], &phi[state_id], 1), 1.0 / curr_norm).real;

		if (adr[state_id] > max_val_direct)
		{
			max_val_direct = adr[state_id];
			max_index_direct = state_id;
		}

		if (stationary == 1)
		{
			phi_n[state_id].real = phi[state_id].real / sqrt(curr_norm);
			phi_n[state_id].imag = phi[state_id].imag / sqrt(curr_norm);

			phi_st[state_id].real = 0.0;
			phi_st[state_id].imag = 0.0;
		}
	}

	curr_mean = get_init_mean(periodic, N, max_index_direct, adr);

	if (stationary == 1)
	{
		cblas_zgemv (CblasRowMajor, CblasNoTrans, N, N, &ONE, ev_t_hamiltonian_sorted, N, phi_n, 1, &ZERO, phi_st, 1);
		max_val_stationary = 0.0;
		max_index_stationary = 0;
		for(int state_id = 0; state_id < N; state_id++)
		{
			adr_stationary[state_id] = Complex_mul(Complex_scalar_mul(&phi_st[state_id], &phi_st[state_id], 1), 1.0).real;
			
			if (adr_stationary[state_id] > max_val_stationary)
			{
				max_val_stationary = adr_stationary[state_id];
				max_index_stationary = state_id;
			}
		}

		curr_mean_stationary = get_init_mean_stationary(periodic, N, max_index_stationary, adr_stationary);
	}

	while (curr_mean < double(N) * mean_low_limit || curr_mean > double(N) * mean_high_limit)
	{
		qj_propagate_one_period(phi, &eta, head, rnd_stream, N, trajectory_id, 0, 0, 0, periodic, stationary);
		trans_process_end_period[trajectory_id]++;

		curr_norm = norm_vector2(phi, N);
		max_val_direct = 0.0;
		max_index_direct = 0;
		for (int state_id = 0; state_id < N; state_id++)
		{
			adr[state_id] = Complex_mul(Complex_scalar_mul(&phi[state_id], &phi[state_id], 1), 1.0 / curr_norm).real;

			if (adr[state_id] > max_val_direct)
			{
				max_val_direct = adr[state_id];
				max_index_direct = state_id;
			}

			if (stationary == 1)
			{
				phi_n[state_id].real = phi[state_id].real / sqrt(curr_norm);
				phi_n[state_id].imag = phi[state_id].imag / sqrt(curr_norm);

				phi_st[state_id].real = 0.0;
				phi_st[state_id].imag = 0.0;
			}
		}

		curr_mean = get_init_mean(periodic, N, max_index_direct, adr);

		if (stationary == 1)
		{
			cblas_zgemv(CblasRowMajor, CblasNoTrans, N, N, &ONE, ev_t_hamiltonian_sorted, N, phi_n, 1, &ZERO, phi_st, 1);
			max_val_stationary = 0.0;
			max_index_stationary = 0;
			for (int state_id = 0; state_id < N; state_id++)
			{
				adr_stationary[state_id] = Complex_mul(Complex_scalar_mul(&phi_st[state_id], &phi_st[state_id], 1), 1.0).real;

				if (adr_stationary[state_id] > max_val_stationary)
				{
					max_val_stationary = adr_stationary[state_id];
					max_index_stationary = state_id;
				}
			}

			curr_mean_stationary = get_init_mean_stationary(periodic, N, max_index_stationary, adr_stationary);
		}
	}

	mean_start[trajectory_id] = curr_mean;
	mean[trajectory_id] = curr_mean;
	dispersion[trajectory_id] = 0.0;
	m2[trajectory_id] = 0.0;

	if (stationary == 1)
	{
		mean_start_stationary[trajectory_id] = curr_mean_stationary;
		mean_stationary[trajectory_id] = curr_mean_stationary;
		dispersion_stationary[trajectory_id] = 0.0;
		m2_stationary[trajectory_id] = 0.0;
		max_id_stationary[trajectory_id] = max_index_stationary;
	}

	etas[trajectory_id] = eta;

	if (stationary == 1)
	{
		phi_n = NULL;
		phi_st = NULL;
	}
}


void single_trajectory_hist_init_prop(
	split * head,
	VSLStreamStatePtr rnd_stream,
	MKL_Complex16 * phi,
	double * adr,
	int num_periods_in_trans_proc,
	int trajectory_id,
	int init_state_id
)
{
	int N = head->N;
	set_init_conditions(phi, N, init_state_id);
	double current_norm = 0.0;

	double eta = 0.0;
	while (eta == 0.0)
	{
		vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, rnd_stream, 1, &eta, 0, 1);
	}


	for(int period_id = 0; period_id < num_periods_in_trans_proc; period_id++)
	{
		qj_propagate_one_period(phi, &eta, head, rnd_stream, N, trajectory_id, 0, 0, 0, 0, 0);
	}

	etas[trajectory_id] = eta;

	current_norm = norm_vector2(phi, N);
	for(int state_id = 0; state_id < N; state_id++)
	{
		adr[state_id] = Complex_mul(Complex_scalar_mul(&phi[state_id], &phi[state_id], 1), 1.0 / current_norm).real;
	}
}


void dump_propagation_data(int N,
						   int num_omp_threads,
						   int num_trajectories,
						   char* write_type,
						   int calc_characteristics,
						   int dump_rho)
{
	char abs_rho_diag_file_name[512];
	sprintf(abs_rho_diag_file_name, "abs_diag_rho_num_trajectories_%d.txt", num_trajectories);

	char rho_file_name[512];
	sprintf(rho_file_name, "rho_period_%d.txt", period);

	char periods_file_name[] = "periods.txt";

	FILE * abs_rho_diag_file = fopen(abs_rho_diag_file_name, write_type);
	double norm = 0.0;
	for (int state_id = 0; state_id < N; state_id++)
	{
		double abs_rho_diag_val = 0.0;
		for (int thread_id = 0; thread_id < num_omp_threads; thread_id ++)
		{
			abs_rho_diag_val += abs_rho_diag_omp[thread_id * (N) + state_id];
		}
		fprintf(abs_rho_diag_file, "%0.18le\n", abs_rho_diag_val);
		norm += abs_rho_diag_val;
	}
	printf("diff norm: %0.18le\n", 1.0 - norm);
	fclose(abs_rho_diag_file);

	FILE * rho_file;
	if (dump_rho)
	{
		rho_file = fopen(rho_file_name, "w");
	}
	for (int state_id_1 = 0; state_id_1 < N; state_id_1++)
	{
		for (int state_id_2 = 0; state_id_2 < N; state_id_2++)
		{
			double real_part = 0.0;
			double imag_part = 0.0;
			for (int thread_id = 0; thread_id < num_omp_threads; thread_id ++)
			{
				real_part += rho_omp[thread_id * (N * N) + state_id_1 * (N) + state_id_2].real;
				imag_part += rho_omp[thread_id * (N * N) + state_id_1 * (N) + state_id_2].imag;
			}

			rho_curr[state_id_1 * (N) + state_id_2].real = real_part;
			rho_curr[state_id_1 * (N) + state_id_2].imag = imag_part;

			if (dump_rho)
			{
				fprintf(rho_file, "%0.18le %0.18le\n", real_part, imag_part);
			}
		}
	}
	if (dump_rho)
	{
		fclose(rho_file);
	}

	FILE * periods_file = fopen(periods_file_name, write_type);
	fprintf(periods_file, "%d\n", period);
	fclose(periods_file);

	if (calc_characteristics)
	{
		calc_and_dump_characteristics(N, write_type);
	}

}

void dump_after_data(
	int N,
	int num_trajectories,
	char* write_type,
	int stationary
)
{
	for (int trajectory_id = 0; trajectory_id < num_trajectories; trajectory_id++)
	{
		char abs_rho_diag_file_name[512];
		sprintf(abs_rho_diag_file_name, "after_abs_diag_rho_trajectory_%d.txt", trajectory_id);
		FILE * abs_rho_diag_file = fopen(abs_rho_diag_file_name, write_type);
		for (int state_id = 0; state_id < N; state_id++)
		{
			fprintf(abs_rho_diag_file, "%0.18le\n", abs_rho_diag_all[trajectory_id * N + state_id]);
		}
		fclose(abs_rho_diag_file);

		if (stationary == 1)
		{
			char abs_rho_diag_stationary_file_name[512];
			sprintf(abs_rho_diag_stationary_file_name, "after_abs_diag_rho_stationary_trajectory_%d.txt", trajectory_id);
			FILE * abs_rho_diag_stationary_file = fopen(abs_rho_diag_stationary_file_name, write_type);
			for (int state_id = 0; state_id < N; state_id++)
			{
				fprintf(abs_rho_diag_stationary_file, "%0.18le\n", abs_rho_diag_all_stationary[trajectory_id * N + state_id]);
			}
			fclose(abs_rho_diag_stationary_file);
		}
	}
}


void dump_statistics_data(
	int N,
	int num_trajectories,
	char* write_type,
	int dump_rho,
	int avg_dump,
	int dump_characteristics,
	int stationary
)
{
	if (avg_dump == 0)
	{
		for (int trajectory_id = 0; trajectory_id < num_trajectories; trajectory_id++)
		{
			if (dump_rho)
			{
				char abs_rho_diag_file_name[512];
				sprintf(abs_rho_diag_file_name, "abs_diag_rho_trajectory_%d.txt", trajectory_id);
				FILE * abs_rho_diag_file = fopen(abs_rho_diag_file_name, write_type);
				for (int state_id = 0; state_id < N; state_id++)
				{
					fprintf(abs_rho_diag_file, "%0.18le\n", abs_rho_diag_all[trajectory_id * N + state_id]);
				}
				fclose(abs_rho_diag_file);

				if (stationary == 1)
				{
					char abs_rho_diag_stationary_file_name[512];
					sprintf(abs_rho_diag_stationary_file_name, "abs_diag_rho_stationary_trajectory_%d.txt", trajectory_id);
					FILE * abs_rho_diag_stationary_file = fopen(abs_rho_diag_stationary_file_name, write_type);
					for (int state_id = 0; state_id < N; state_id++)
					{
						fprintf(abs_rho_diag_stationary_file, "%0.18le\n", abs_rho_diag_all_stationary[trajectory_id * N + state_id]);
					}
					fclose(abs_rho_diag_stationary_file);
				}
			}

			char characteristics_file_name[512];
			sprintf(characteristics_file_name, "characteristics_trajectory_%d.txt", trajectory_id);
			FILE * characteristics_file = fopen(characteristics_file_name, write_type);
			if (stationary == 1)
			{
				fprintf(characteristics_file, "%0.18le %0.18le %0.18le %0.18le %0.18le %0.18le %d\n", mean[trajectory_id], dispersion[trajectory_id], m2[trajectory_id], mean_stationary[trajectory_id], dispersion_stationary[trajectory_id], m2_stationary[trajectory_id], max_id_stationary[trajectory_id]);
			}
			else
			{
				fprintf(characteristics_file, "%0.18le %0.18le %0.18le\n", mean[trajectory_id], dispersion[trajectory_id], m2[trajectory_id]);
			}
			fclose(characteristics_file);
		}
	}
	else if (avg_dump == 1)
	{
		double * abs_rho_diag_avg;
		double * abs_rho_diag_avg_stationary;

		double mean_avg = 0.0;
		double dispersion_avg = 0.0;
		double m2_avg = 0.0;

		double mean_stationary_avg = 0.0;
		double dispersion_stationary_avg = 0.0;
		double m2_stationary_avg = 0.0;
		double max_id_stationary_avg = 0.0;

		if (dump_rho)
		{
			abs_rho_diag_avg = new double[N];
			if (stationary == 1)
			{
				abs_rho_diag_avg_stationary = new double[N];
			}
			for (int state_id = 0; state_id < N; state_id++)
			{
				abs_rho_diag_avg[state_id] = 0;
				if (stationary == 1)
				{
					abs_rho_diag_avg_stationary[state_id] = 0;
				}
			}
		}

		for (int trajectory_id = 0; trajectory_id < num_trajectories; trajectory_id++)
		{
			if (dump_rho)
			{
				for (int state_id = 0; state_id < N; state_id++)
				{
					abs_rho_diag_avg[state_id] += abs_rho_diag_all[trajectory_id * N + state_id] / double(num_trajectories);

					if (stationary == 1)
					{
						abs_rho_diag_avg_stationary[state_id] += abs_rho_diag_all_stationary[trajectory_id * N + state_id] / double(num_trajectories);
					}
				}
			}

			mean_avg += mean[trajectory_id] / double(num_trajectories);
			dispersion_avg += dispersion[trajectory_id] / double(num_trajectories);
			m2_avg += m2[trajectory_id] / double(num_trajectories);

			if (stationary == 1)
			{
				mean_stationary_avg += mean_stationary[trajectory_id] / double(num_trajectories);
				dispersion_stationary_avg += dispersion_stationary[trajectory_id] / double(num_trajectories);
				m2_stationary_avg += m2_stationary[trajectory_id] / double(num_trajectories);
				max_id_stationary_avg += double(max_id_stationary[trajectory_id]) / double(num_trajectories);
			}
		}

		if (dump_rho)
		{
			char abs_rho_diag_file_name[512];
			sprintf(abs_rho_diag_file_name, "abs_diag_rho_avg_%d.txt", num_trajectories);
			FILE * abs_rho_diag_file = fopen(abs_rho_diag_file_name, write_type);
			for (int state_id = 0; state_id < N; state_id++)
			{
				fprintf(abs_rho_diag_file, "%0.18le\n", abs_rho_diag_avg[state_id]);
			}
			fclose(abs_rho_diag_file);

			if (stationary == 1)
			{
				char abs_rho_diag_stationary_file_name[512];
				sprintf(abs_rho_diag_stationary_file_name, "abs_diag_rho_stationary_avg_%d.txt", num_trajectories);
				FILE * abs_rho_diag_stationary_file = fopen(abs_rho_diag_stationary_file_name, write_type);
				for (int state_id = 0; state_id < N; state_id++)
				{
					fprintf(abs_rho_diag_stationary_file, "%0.18le\n", abs_rho_diag_avg_stationary[state_id]);
				}
				fclose(abs_rho_diag_stationary_file);
			}
		}

		char characteristics_file_name[512];
		sprintf(characteristics_file_name, "characteristics_avg_%d.txt", num_trajectories);
		FILE * characteristics_file = fopen(characteristics_file_name, write_type);
		if (stationary == 1)
		{
			fprintf(characteristics_file, "%0.18le %0.18le %0.18le %0.18le %0.18le %0.18le %0.18le\n", mean_avg, dispersion_avg, m2_avg, mean_stationary_avg, dispersion_stationary_avg, m2_stationary_avg, max_id_stationary_avg);
		}
		else
		{
			fprintf(characteristics_file, "%0.18le %0.18le %0.18le\n", mean_avg, dispersion_avg, m2_avg);
		}
		fclose(characteristics_file);

		if (dump_rho)
		{
			delete[] abs_rho_diag_avg;
			if (stationary == 1)
			{
				delete[] abs_rho_diag_avg_stationary;
			}
		}
	}

	char periods_file_name[] = "periods.txt";
	FILE * periods_file = fopen(periods_file_name, write_type);
	fprintf(periods_file, "%d\n", period);
	fclose(periods_file);
}

void dump_aux_statistics_data(
	int N,
	int num_periods,
	int num_trajectories,
	char* write_type,
	int avg_dump,
	int dump_characteristics,
	int btw_jump_times,
	int stationary
)
{
	if (avg_dump == 0)
	{
		for (int trajectory_id = 0; trajectory_id < num_trajectories; trajectory_id++)
		{
			char characteristics_file_name[512];
			sprintf(characteristics_file_name, "aux_characteristics_trajectory_%d.txt", trajectory_id);
			FILE * characteristics_file = fopen(characteristics_file_name, write_type);
			fprintf(characteristics_file, "%d %0.18le\n", trans_process_end_period[trajectory_id], mean_start[trajectory_id]);
			fclose(characteristics_file);

			char abs_rho_diag_file_name[512];
			sprintf(abs_rho_diag_file_name, "avg_diag_rho_trajectory_%d.txt", trajectory_id);
			FILE * abs_rho_diag_file = fopen(abs_rho_diag_file_name, write_type);
			for (int state_id = 0; state_id < N; state_id++)
			{
				fprintf(abs_rho_diag_file, "%0.18le\n", avg_rho_diag_all[trajectory_id * N + state_id] / (double(num_periods-1)));
			}
			fclose(abs_rho_diag_file);
		}
	}

	if (dump_characteristics >= 1)
	{
		if (num_att_peaks > 0)
		{
			char file_name[512];
			FILE * file;

			sprintf(file_name, "peaks_limits.txt");
			file = fopen(file_name, write_type);
			for (int peak_id = 0; peak_id < num_att_peaks; peak_id++)
			{
				fprintf(file, "%0.18le %0.18le\n", begins_att_peaks[peak_id], ends_att_peaks[peak_id]);
			}
			fclose(file);

			if (avg_dump == 0)
			{
				for (int trajectory_id = 0; trajectory_id < num_trajectories; trajectory_id++)
				{
					for (int peak_id = 0; peak_id < num_att_peaks; peak_id++)
					{
						sprintf(file_name, "sticking_times_trajectory_%d_peak_id_%d.txt", trajectory_id, peak_id);
						file = fopen(file_name, write_type);
						if (!sticking_times[trajectory_id][peak_id].empty())
						{
							for (int i = 0; i < sticking_times[trajectory_id][peak_id].size(); i ++)
							{
								fprintf(file, "%d\n",sticking_times[trajectory_id][peak_id][i]);
							}
						}
						fclose(file);
					}
				}
			}
			else
			{
				for (int peak_id = 0; peak_id < num_att_peaks; peak_id++)
				{
					sprintf(file_name, "sticking_times_all_trajectories_%d_peak_id_%d.txt", num_trajectories, peak_id);
					file = fopen(file_name, write_type);

					for (int trajectory_id = 0; trajectory_id < num_trajectories; trajectory_id++)
					{
						if (!sticking_times[trajectory_id][peak_id].empty())
						{
							for (int i = 0; i < sticking_times[trajectory_id][peak_id].size(); i ++)
							{
								fprintf(file, "%d\n",sticking_times[trajectory_id][peak_id][i]);
							}
						}
					}

					fclose(file);
				}
			}
		}
	}

	if (stationary == 1 && 	dump_characteristics >= 1)
	{
		if (num_att_peaks_stationary > 0)
		{
			char file_name[512];
			FILE * file;

			sprintf(file_name, "peaks_limits_stationary.txt");
			file = fopen(file_name, write_type);
			for (int peak_id = 0; peak_id < num_att_peaks_stationary; peak_id++)
			{
				fprintf(file, "%0.18le %0.18le\n", begins_att_peaks_stationary[peak_id], ends_att_peaks_stationary[peak_id]);
			}
			fclose(file);

			if (avg_dump == 0)
			{
				for (int trajectory_id = 0; trajectory_id < num_trajectories; trajectory_id++)
				{
					for (int peak_id = 0; peak_id < num_att_peaks_stationary; peak_id++)
					{
						sprintf(file_name, "sticking_times_trajectory_%d_peak_id_%d_stationary.txt", trajectory_id, peak_id);
						file = fopen(file_name, write_type);
						if (!sticking_times_stationary[trajectory_id][peak_id].empty())
						{
							for (int i = 0; i < sticking_times_stationary[trajectory_id][peak_id].size(); i ++)
							{
								fprintf(file, "%d\n",sticking_times_stationary[trajectory_id][peak_id][i]);
							}
						}
						fclose(file);
					}
				}
			}
			else
			{
				for (int peak_id = 0; peak_id < num_att_peaks_stationary; peak_id++)
				{
					sprintf(file_name, "sticking_times_all_trajectories_%d_peak_id_%d_stationary.txt", num_trajectories, peak_id);
					file = fopen(file_name, write_type);

					for (int trajectory_id = 0; trajectory_id < num_trajectories; trajectory_id++)
					{
						if (!sticking_times_stationary[trajectory_id][peak_id].empty())
						{
							for (int i = 0; i < sticking_times_stationary[trajectory_id][peak_id].size(); i ++)
							{
								fprintf(file, "%d\n",sticking_times_stationary[trajectory_id][peak_id][i]);
							}
						}
					}

					fclose(file);
				}
			}
		}
	}

	if (btw_jump_times > 0)
	{
		char file_name[512];
		FILE * file;

		if (avg_dump == 0)
		{
			for (int trajectory_id = 0; trajectory_id < num_trajectories; trajectory_id++)
			{
				sprintf(file_name, "jump_times_trajectory_%d.txt", trajectory_id);
				file = fopen(file_name, write_type);
				if (!jump_times[trajectory_id].empty())
				{
					for (int i = 0; i < jump_times[trajectory_id].size(); i ++)
					{
						fprintf(file, "%0.18le\n", jump_times[trajectory_id][i]);
					}
				}
				fclose(file);
			}
		}
		else
		{
			sprintf(file_name, "jump_times_all_trajectories_%d.txt", num_trajectories);
			file = fopen(file_name, write_type);
			for (int trajectory_id = 0; trajectory_id < num_trajectories; trajectory_id++)
			{
				
				if (!jump_times[trajectory_id].empty())
				{
					for (int i = 0; i < jump_times[trajectory_id].size(); i ++)
					{
						fprintf(file, "%0.18le\n", jump_times[trajectory_id][i] - double(trans_process_end_period[trajectory_id]));
					}
				}
			}
			fclose(file);
		}
	}
}

void dump_hist_data(
	int num_trajectories,
	char* write_type,
	int num_e_intervals,
	int num_mc_intervals
)
{
	
	char file_name[512];
	FILE * file;

	sprintf(file_name, "energy_intervals.txt");
	file = fopen(file_name, write_type);
	for (int e_interval_id = 0; e_interval_id < num_e_intervals; e_interval_id ++)
	{
		fprintf(file, "%0.18le \n", energy_intervals[e_interval_id]);
	}
	fclose(file);

	sprintf(file_name, "mass_centers_intervals.txt");
	file = fopen(file_name, write_type);
	for (int mc_interval_id = 0; mc_interval_id < num_mc_intervals; mc_interval_id ++)
	{
		fprintf(file, "%0.18le \n", mass_centers_intervals[mc_interval_id]);
	}
	fclose(file);



	for (int trajectory_id = 0; trajectory_id < num_trajectories; trajectory_id++)
	{
		sprintf(file_name, "histogramm_trajectory_%d.txt", trajectory_id);
		file = fopen(file_name, write_type);

		for (int e_interval_id = 0; e_interval_id < num_e_intervals; e_interval_id ++)
		{
			for (int mc_interval_id = 0; mc_interval_id < num_mc_intervals; mc_interval_id ++)
			{
				fprintf(file, "%d \n", histogramm[trajectory_id * (num_e_intervals * num_mc_intervals) + e_interval_id * num_mc_intervals + mc_interval_id]);
			}
		}

		fclose(file);
	}
}


void single_trajectory_dump_prop(
	split * head,
	VSLStreamStatePtr rnd_stream,
	MKL_Complex16 * phi,
	double * abs_rho_diag,
	MKL_Complex16 * rho,
	int num_trajectories,
	int trajectory_id,
	int begin_propagation_loop,
	int end_propagation_loop,
	int thread_id
)
{
	int N = head->N;

	double eta = etas[trajectory_id];

	for(int period_id = begin_propagation_loop; period_id < end_propagation_loop; period_id++)
	{
		qj_propagate_one_period(phi, &eta, head, rnd_stream, N, trajectory_id, 0, 0, 0, 0, 0);
	}

	etas[trajectory_id] = eta;

	double norm = norm_vector2(phi, N);
	for(int state_id_1 = 0; state_id_1 < N; state_id_1++)
	{
		abs_rho_diag[state_id_1] += Complex_mul(Complex_scalar_mul(&phi[state_id_1], &phi[state_id_1], 1), 1.0 / norm).real / num_trajectories;

		for (int state_id_2 = 0; state_id_2 < N; state_id_2++)
		{
			rho[state_id_1 * (N) + state_id_2].real += Complex_mul(Complex_scalar_mul(&phi[state_id_1], &phi[state_id_2], 1), 1.0 / norm).real / num_trajectories;
			rho[state_id_1 * (N) + state_id_2].imag += Complex_mul(Complex_scalar_mul(&phi[state_id_1], &phi[state_id_2], 1), 1.0 / norm).imag / num_trajectories;
		}
	}
}

void single_trajectory_statistics_dump_prop(
	split * head,
	VSLStreamStatePtr rnd_stream,
	MKL_Complex16 * phi,
	double * adr,
	int trajectory_id,
	int begin_propagation_loop,
	int end_propagation_loop,
	int dump_characteristics,
	int btw_jump_times,
	double * adr_stationary,
	int stationary,
	int periodic, 
	int deep_characteristics,
	int mc_specific,
	int mc_type
)
{
	if (end_propagation_loop > begin_propagation_loop)
	{
		int N = head->N;
		double curr_norm = 0.0;

		double max_val_direct = 0.0;
		int max_index_direct = 0;

		double max_val_stationary = 0.0;
		int max_index_stationary = 0;

		double curr_eta = etas[trajectory_id];

		double * avg_adr = &avg_rho_diag_all[trajectory_id * N];

		MKL_Complex16 * phi_n = NULL;
		MKL_Complex16 * phi_st = NULL;
		if (stationary == 1)
		{
			phi_n = &phi_normed[trajectory_id * N];
			phi_st = &phi_stationary[trajectory_id * N];
		}

		double eta = etas[trajectory_id];

		for(int period_id = begin_propagation_loop; period_id < end_propagation_loop; period_id++)
		{
			qj_propagate_one_period(phi, &eta, head, rnd_stream, N, trajectory_id, btw_jump_times, deep_characteristics, mc_specific, periodic, stationary);

			curr_norm = norm_vector2(phi, N);
			
			max_val_direct = 0.0;
			max_index_direct = 0;
			for(int state_id = 0; state_id < N; state_id++)
			{
				adr[state_id] = Complex_mul(Complex_scalar_mul(&phi[state_id], &phi[state_id], 1), 1.0 / curr_norm).real;

				if (adr[state_id] > max_val_direct)
				{
					max_val_direct = adr[state_id];
					max_index_direct = state_id;
				}

				if (stationary == 1)
				{
					phi_n[state_id].real = phi[state_id].real / sqrt(curr_norm);
					phi_n[state_id].imag = phi[state_id].imag / sqrt(curr_norm);

					phi_st[state_id].real = 0.0;
					phi_st[state_id].imag = 0.0;
				}
			}


			double curr_mean = 0.0;
			if (mc_type == 0)
			{
				curr_mean = mean[trajectory_id];
			}
			else
			{
				curr_mean = get_mean(periodic, N, max_index_direct, adr, trajectory_id, mc_specific);
			}

			double curr_mean_on_period = get_mean_on_period(N, curr_mean);
			double curr_dispersion = get_dispersion(curr_mean, trajectory_id);
			double curr_m2 = get_m2(periodic, N, curr_mean_on_period, adr);

			mean[trajectory_id] = curr_mean;
			dispersion[trajectory_id] = curr_dispersion;
			m2[trajectory_id] = curr_m2;

			if (dump_characteristics >= 1)
			{
				st_processing_direct(periodic, N, curr_mean_on_period, trajectory_id);
			}

			if (stationary == 1)
			{
				cblas_zgemv (CblasRowMajor, CblasNoTrans, N, N, &ONE, ev_t_hamiltonian_sorted, N, phi_n, 1, &ZERO, phi_st, 1);
				max_val_stationary = 0.0;
				max_index_stationary = 0;
				for(int state_id = 0; state_id < N; state_id++)
				{
					adr_stationary[state_id] = Complex_mul(Complex_scalar_mul(&phi_st[state_id], &phi_st[state_id], 1), 1.0).real;;

					if (adr_stationary[state_id] > max_val_stationary)
					{
						max_val_stationary = adr_stationary[state_id];
						max_index_stationary = state_id;
					}
				}

				double current_mean_stationary = get_mean_stationary(periodic, N, max_index_stationary, adr_stationary, trajectory_id, mc_specific);
				double current_mean_on_period_stationary = get_mean_on_period(N, current_mean_stationary);
				double current_dispersion_stationary = get_dispersion_stationary(current_mean_stationary, trajectory_id);
				double curr_m2_stationary = get_m2(periodic, N, current_mean_on_period_stationary, adr_stationary);

				mean_stationary[trajectory_id] = current_mean_stationary;
				dispersion_stationary[trajectory_id] = current_dispersion_stationary;
				m2_stationary[trajectory_id] = curr_m2_stationary;

				max_id_stationary[trajectory_id] = max_index_stationary;

				if (dump_characteristics >= 1)
				{
					st_processing_stationary(periodic, N, max_index_stationary, trajectory_id);
				}
			}
		}

		for (int state_id = 0; state_id < N; state_id++)
		{
			avg_adr[state_id] += adr[state_id];
		}

		etas[trajectory_id] = eta;

		if (stationary == 1)
		{
			phi_n = NULL;
			phi_st = NULL;
		}
	}
}

void single_trajectory_hist_dump_prop(
	split * head,
	VSLStreamStatePtr rnd_stream,
	MKL_Complex16 * phi,
	double * adr,
	MKL_Complex16 * sm,
	int * hist,
	int trajectory_id,
	int begin_propagation_loop,
	int end_propagation_loop,
	int num_e_intervals,
	int num_mc_intervals,
	int periodic
)
{
	if (end_propagation_loop > begin_propagation_loop)
	{
		int N = head->N;
		double current_norm = 0.0;

		double max_val = 0.0;
		int max_index = 0;

		double curr_mean = 0.0;
		double curr_mean_on_period = 0;

		double eta = etas[trajectory_id];

		for(int period_id = begin_propagation_loop; period_id < end_propagation_loop; period_id++)
		{
			qj_propagate_one_period(phi, &eta, head, rnd_stream, N, trajectory_id, 0, 0, 0, periodic, 0);
		}

		etas[trajectory_id] = eta;

		current_norm = norm_vector2(phi, N);

		max_val = 0.0;
		max_index = 0;
		for(int state_id = 0; state_id < N; state_id++)
		{
			adr[state_id] = Complex_mul(Complex_scalar_mul(&phi[state_id], &phi[state_id], 1), 1.0 / current_norm).real;

			if (adr[state_id] > max_val)
			{
				max_val = adr[state_id];
				max_index = state_id;
			}
		}

		curr_mean = get_init_mean(periodic, N, max_index, adr);
		curr_mean_on_period = get_mean_on_period(N, curr_mean);

		for (int state_id_1 = 0; state_id_1 < N; state_id_1++)
		{
			MKL_Complex16 tmp;
			tmp.real = 0;
			tmp.imag = 0;
			for (int state_id_2 = 0; state_id_2 < N; state_id_2++)
			{
				tmp.real += (hamiltonian[state_id_1 * N + state_id_2].real * phi[state_id_2].real / sqrt(current_norm) - hamiltonian[state_id_1 * N + state_id_2].imag * phi[state_id_2].imag / sqrt(current_norm));
				tmp.imag += (hamiltonian[state_id_1 * N + state_id_2].real * phi[state_id_2].imag / sqrt(current_norm) + hamiltonian[state_id_1 * N + state_id_2].imag * phi[state_id_2].real / sqrt(current_norm));
			}

			sm[state_id_1].real = tmp.real;
			sm[state_id_1].imag = tmp.imag;
		}

		double current_e = 0.0;
		for (int state_id_1 = 0; state_id_1 < N; state_id_1++)
		{
			current_e += (phi[state_id_1].real / sqrt(current_norm) * sm[state_id_1].real + phi[state_id_1].imag / sqrt(current_norm) * sm[state_id_1].imag);
		}

		if ((current_e < e_max) && (current_e > e_min) && (curr_mean_on_period > mc_min) && (curr_mean_on_period < mc_max))
		{
			int energy_id = floor((current_e - e_min) * num_e_intervals / (e_max - e_min + 0.000001));
			int mass_center_id = floor((curr_mean_on_period - mc_min) * num_mc_intervals / (mc_max - mc_min + 0.000001));

			hist[energy_id * num_mc_intervals + mass_center_id] ++;
		}
	}
}

void init_aux_data(char aux_file_name[], int N)
{
	rho_curr = new MKL_Complex16[N*N];
	rho_curr_and = new MKL_Complex16[N*N];
	hamiltonian = new MKL_Complex16[N*N];
	eg_hamiltonian = new MKL_Complex16[N];
	ev_hamiltonian = new MKL_Complex16[N*N];
	ev_t_hamiltonian = new MKL_Complex16[N*N];

	for (int state_id_1 = 0; state_id_1 < N; state_id_1++)
	{
		for (int state_id_2 = 0; state_id_2 < N; state_id_2++)
		{
			rho_curr[state_id_1 * (N) + state_id_2].real = 0.0;
			rho_curr[state_id_1 * (N) + state_id_2].imag = 0.0;

			rho_curr_and[state_id_1 * (N) + state_id_2].real = 0.0;
			rho_curr_and[state_id_1 * (N) + state_id_2].imag = 0.0;

			hamiltonian[state_id_1 * (N) + state_id_2].real = 0.0;
			hamiltonian[state_id_1 * (N) + state_id_2].imag = 0.0;

			ev_hamiltonian[state_id_1 * (N) + state_id_2].real = 0.0;
			ev_hamiltonian[state_id_1 * (N) + state_id_2].imag = 0.0;

			ev_t_hamiltonian[state_id_1 * (N) + state_id_2].real = 0.0;
			ev_t_hamiltonian[state_id_1 * (N) + state_id_2].imag = 0.0;
		}

		eg_hamiltonian[state_id_1].real = 0.0;
		eg_hamiltonian[state_id_1].imag = 0.0;
	}

	FILE * aux_file = fopen(aux_file_name, "rb");
	fread(hamiltonian, sizeof(MKL_Complex16), N * N, aux_file);
	fread(eg_hamiltonian, sizeof(MKL_Complex16), N, aux_file);
	fread(ev_hamiltonian, sizeof(MKL_Complex16), N * N, aux_file);
	fread(ev_t_hamiltonian, sizeof(MKL_Complex16), N * N, aux_file);
	fclose(aux_file);
}

void init_aux_statistics_data(char aux_file_name[],
							  int N,
							  int num_trajectories,
							  int borders_type,
							  int dump_characteristics,
							  int stationary)
{
	hamiltonian = new MKL_Complex16[N*N];
	eg_hamiltonian = new MKL_Complex16[N];
	ev_hamiltonian = new MKL_Complex16[N*N];
	ev_t_hamiltonian = new MKL_Complex16[N*N];

	ev_hamiltonian_sorted = new MKL_Complex16[N*N];
	ev_t_hamiltonian_sorted = new MKL_Complex16[N*N];

	for (int state_id_1 = 0; state_id_1 < N; state_id_1++)
	{
		for (int state_id_2 = 0; state_id_2 < N; state_id_2++)
		{
			hamiltonian[state_id_1 * (N) + state_id_2].real = 0.0;
			hamiltonian[state_id_1 * (N) + state_id_2].imag = 0.0;

			ev_hamiltonian[state_id_1 * (N) + state_id_2].real = 0.0;
			ev_hamiltonian[state_id_1 * (N) + state_id_2].imag = 0.0;

			ev_t_hamiltonian[state_id_1 * (N) + state_id_2].real = 0.0;
			ev_t_hamiltonian[state_id_1 * (N) + state_id_2].imag = 0.0;


			ev_hamiltonian_sorted[state_id_1 * (N) + state_id_2].real = 0.0;
			ev_hamiltonian_sorted[state_id_1 * (N) + state_id_2].imag = 0.0;

			ev_t_hamiltonian_sorted[state_id_1 * (N) + state_id_2].real = 0.0;
			ev_t_hamiltonian_sorted[state_id_1 * (N) + state_id_2].imag = 0.0;
		}

		eg_hamiltonian[state_id_1].real = 0.0;
		eg_hamiltonian[state_id_1].imag = 0.0;
	}

	FILE * aux_file = fopen(aux_file_name, "rb");
	fread(hamiltonian, sizeof(MKL_Complex16), N * N, aux_file);
	fread(eg_hamiltonian, sizeof(MKL_Complex16), N, aux_file);
	fread(ev_hamiltonian, sizeof(MKL_Complex16), N * N, aux_file);
	fread(ev_t_hamiltonian, sizeof(MKL_Complex16), N * N, aux_file);

	fread(ev_hamiltonian_sorted, sizeof(MKL_Complex16), N * N, aux_file);
	fread(ev_t_hamiltonian_sorted, sizeof(MKL_Complex16), N * N, aux_file);

	if (dump_characteristics > 0)
	{

		fread(&num_att_peaks, sizeof(int), 1, aux_file);
		if (num_att_peaks > 0)
		{
			begins_att_peaks = new double[num_att_peaks];
			ends_att_peaks = new double[num_att_peaks];
			fread(begins_att_peaks, sizeof(double), num_att_peaks, aux_file);
			fread(ends_att_peaks, sizeof(double), num_att_peaks, aux_file);
		}

		fread(&num_att_peaks_stationary, sizeof(int), 1, aux_file);
		if (num_att_peaks_stationary > 0)
		{
			begins_att_peaks_stationary = new double[num_att_peaks_stationary];
			ends_att_peaks_stationary = new double[num_att_peaks_stationary];
			fread(begins_att_peaks_stationary, sizeof(double), num_att_peaks_stationary, aux_file);
			fread(ends_att_peaks_stationary, sizeof(double), num_att_peaks_stationary, aux_file);
		}
	}

	fclose(aux_file);

	if (dump_characteristics > 0)
	{
		if (borders_type == 1)
		{
			if (num_att_peaks > 0)
			{
				if ((fabs(double(begins_att_peaks[0]) + 0.5) < 1.0e-10) && ((fabs(double(ends_att_peaks[num_att_peaks - 1]) - (double(N - 1) + 0.5)) < 1.0e-10)))
				{
					if (num_att_peaks > 1)
					{
						double begin_special_peak = begins_att_peaks[num_att_peaks - 1];
						double end_special_peak = double(N) + (ends_att_peaks[0] - begins_att_peaks[0]) + 0.5;

						for (int peak_id = 0; peak_id < num_att_peaks - 2; peak_id++)
						{
							begins_att_peaks[peak_id] = begins_att_peaks[peak_id + 1];
							ends_att_peaks[peak_id] = ends_att_peaks[peak_id + 1];
						}

						begins_att_peaks[num_att_peaks - 2] = begin_special_peak;
						ends_att_peaks[num_att_peaks - 2] = end_special_peak;
					}

					num_att_peaks = num_att_peaks - 1;
				}
			}

			if (num_att_peaks_stationary > 0)
			{
				if ((fabs(begins_att_peaks_stationary[0] + 0.5) < 1.0e-10) && ((fabs(ends_att_peaks_stationary[num_att_peaks_stationary - 1] - (double(N - 1) + 0.5)) < 1.0e-10)))
				{
					if (num_att_peaks_stationary > 1)
					{
						double begin_special_peak_stationary = begins_att_peaks_stationary[num_att_peaks_stationary - 1];
						double end_special_peak_stationary = double(N) + (ends_att_peaks_stationary[0] - begins_att_peaks_stationary[0]) + 0.5;

						for (int peak_id = 0; peak_id < num_att_peaks_stationary - 2; peak_id++)
						{
							begins_att_peaks_stationary[peak_id] = begins_att_peaks_stationary[peak_id + 1];
							ends_att_peaks_stationary[peak_id] = ends_att_peaks_stationary[peak_id + 1];
						}

						begins_att_peaks_stationary[num_att_peaks_stationary - 2] = begin_special_peak_stationary;
						ends_att_peaks_stationary[num_att_peaks_stationary - 2] = end_special_peak_stationary;
					}

					num_att_peaks_stationary = num_att_peaks_stationary - 1;
				}
			}
		}

		printf("num_att_peaks: %d\n", num_att_peaks);
		printf("num_att_peaks_stationary: %d\n", num_att_peaks_stationary);


		if (num_att_peaks > 0)
		{
			sticking_times = new vector<int> *[num_trajectories];
			current_peaks_ids = new int[num_trajectories];
			for (int trajectory_id = 0; trajectory_id < num_trajectories; trajectory_id++)
			{
				sticking_times[trajectory_id] = new vector<int>[num_att_peaks];
				current_peaks_ids[trajectory_id] = -1; // not belongs to any peak
			}
		}

		if (num_att_peaks_stationary > 0)
		{
			sticking_times_stationary = new vector<int> *[num_trajectories];
			current_peaks_ids_stationary = new int[num_trajectories];
			for (int trajectory_id = 0; trajectory_id < num_trajectories; trajectory_id++)
			{
				sticking_times_stationary[trajectory_id] = new vector<int>[num_att_peaks_stationary];
				current_peaks_ids_stationary[trajectory_id] = -1; // not belongs to any peak
			}
		}
	}
}

void init_aux_hist_data(
	char aux_file_name[],
	int N
)
{
	hamiltonian = new MKL_Complex16[N*N];
	eg_hamiltonian = new MKL_Complex16[N];
	ev_hamiltonian = new MKL_Complex16[N*N];
	ev_t_hamiltonian = new MKL_Complex16[N*N];

	for (int state_id_1 = 0; state_id_1 < N; state_id_1++)
	{
		for (int state_id_2 = 0; state_id_2 < N; state_id_2++)
		{
			hamiltonian[state_id_1 * (N) + state_id_2].real = 0.0;
			hamiltonian[state_id_1 * (N) + state_id_2].imag = 0.0;

			ev_hamiltonian[state_id_1 * (N) + state_id_2].real = 0.0;
			ev_hamiltonian[state_id_1 * (N) + state_id_2].imag = 0.0;

			ev_t_hamiltonian[state_id_1 * (N) + state_id_2].real = 0.0;
			ev_t_hamiltonian[state_id_1 * (N) + state_id_2].imag = 0.0;
		}

		eg_hamiltonian[state_id_1].real = 0.0;
		eg_hamiltonian[state_id_1].imag = 0.0;
	}

	FILE * aux_file = fopen(aux_file_name, "rb");
	fread(hamiltonian, sizeof(MKL_Complex16), N * N, aux_file);
	fread(eg_hamiltonian, sizeof(MKL_Complex16), N, aux_file);
	fread(ev_hamiltonian, sizeof(MKL_Complex16), N * N, aux_file);
	fread(ev_t_hamiltonian, sizeof(MKL_Complex16), N * N, aux_file);

	fclose(aux_file);
}

void delete_aux_hist_data()
{
	delete[] hamiltonian;
	delete[] eg_hamiltonian;
	delete[] ev_hamiltonian;
	delete[] ev_t_hamiltonian;
}

void delete_aux_statistics_data(int num_trajectories, int dump_characteristics, int stationary)
{
	delete[] hamiltonian;
	delete[] eg_hamiltonian;
	delete[] ev_hamiltonian;
	delete[] ev_t_hamiltonian;
	delete[] ev_hamiltonian_sorted;
	delete[] ev_t_hamiltonian_sorted;

	if (dump_characteristics > 0)
	{

		if (num_att_peaks > 0)
		{
			delete[] begins_att_peaks;
			delete[] ends_att_peaks;
		}

		if (num_att_peaks > 0)
		{
			for (int trajectory_id = 0; trajectory_id < num_trajectories; trajectory_id++)
			{
				delete[] sticking_times[trajectory_id];
			}
			delete[] sticking_times;
			delete[] current_peaks_ids;
		}

		if (stationary == 1)
		{

			if (num_att_peaks_stationary > 0)
			{
				delete[] begins_att_peaks_stationary;
				delete[] ends_att_peaks_stationary;
			}


			if (num_att_peaks_stationary > 0)
			{
				for (int trajectory_id = 0; trajectory_id < num_trajectories; trajectory_id++)
				{
					delete[] sticking_times_stationary[trajectory_id];
				}
				delete[] sticking_times_stationary;
				delete[] current_peaks_ids_stationary;
			}
		}
	}
}


void refresh_aux_data(int N)
{
	for (int state_id_1 = 0; state_id_1 < N; state_id_1++)
	{
		for (int state_id_2 = 0; state_id_2 < N; state_id_2++)
		{
			rho_curr[state_id_1 * (N) + state_id_2].real = 0.0;
			rho_curr[state_id_1 * (N) + state_id_2].imag = 0.0;

			rho_curr_and[state_id_1 * (N) + state_id_2].real = 0.0;
			rho_curr_and[state_id_1 * (N) + state_id_2].imag = 0.0;
		}
	}
}

void delete_aux_data()
{
	delete[] rho_curr;
	delete[] rho_curr_and;
	delete[] hamiltonian;
	delete[] eg_hamiltonian;
	delete[] ev_hamiltonian;
	delete[] ev_t_hamiltonian;
}

void omp_qj(
	char input_file_name[],
	char aux_file_name[],
	int num_periods,
	int num_dumps,
	int dump_type,
	int num_periods_in_trans_proc,
	int num_omp_threads,
	int num_trajectories,
	int rnd_max,
	int rnd_cur,
	int init_state_id,
	int calc_characteristics,
	int dump_rho
)
{
	printf("num_omp_threads: %d\n\n", num_omp_threads);
	omp_set_num_threads(num_omp_threads);

	FILE * input_file;
	input_file = fopen(input_file_name, "rb");
	split * head = create_struct_bin(input_file);
	fclose(input_file);

	int N = head->N;

	VSLStreamStatePtr * rnd_streams = new VSLStreamStatePtr[num_trajectories];
	vslNewStream (&rnd_streams[0], VSL_BRNG_MCG59, 777);
	for(int i = 1; i < num_trajectories; i++)
	{
		vslCopyStream(&rnd_streams[i], rnd_streams[0]);
	}
	for(int i = 0; i < num_trajectories; i++)
	{
		vslLeapfrogStream(rnd_streams[i], rnd_cur + i, rnd_max);
	}

	split * heads = new split[num_omp_threads];
	for(int i = 0; i < num_omp_threads; i++)
	{
		cmp_struct_not_member(head, &heads[i]);
	}

	init_data(N, num_trajectories, num_omp_threads);
	init_aux_data(aux_file_name, N);

	double step = double(num_periods) / double(num_dumps - 1);
	double start = 0.0;
	if (dump_type == 1)
	{
		start = 0.0;
		step = log10(double(num_periods)) / double(num_dumps - 1);
	}

#pragma omp parallel for
	for(int i = 0; i < num_trajectories; i++)
	{
		int thread_id = omp_get_thread_num();

		single_trajectory_init_prop(
			&heads[thread_id],
			rnd_streams[i],
			&phi_global[i * N],
			&abs_rho_diag_omp[thread_id * N],
			&rho_omp[thread_id * N * N],
			num_periods_in_trans_proc,
			num_trajectories,
			i,
			init_state_id,
			thread_id
		);
	}

	for (int trajectory_id = 0; trajectory_id < num_trajectories; trajectory_id++)
	{
		ts[trajectory_id] = 0.0;
	}

	period = 1;
	printf("Period = %d\n", period);

	dump_propagation_data(N, num_omp_threads, num_trajectories, "w", calc_characteristics, dump_rho);
	refresh_dump_data(N, num_omp_threads);
	refresh_aux_data(N);

	double current_dump_limit = start;
	int begin_propagation_loop = period;
	int end_propagation_loop = 0;

	while (begin_propagation_loop != num_periods)
	{
		current_dump_limit += step;
		end_propagation_loop = num_periods;

		if (dump_type == 0)
		{
			if (current_dump_limit < double(end_propagation_loop))
			{
				end_propagation_loop = int(current_dump_limit);
			}
		}
		else if (dump_type == 1)
		{
			int exp_limit = int(pow(10.0, current_dump_limit + 1.0e-10));
			if (exp_limit < end_propagation_loop)
			{
				if (exp_limit <= begin_propagation_loop)
				{
					end_propagation_loop = begin_propagation_loop + 1;
				}
				else
				{
					end_propagation_loop = exp_limit;
				}
			}
		}

		if (end_propagation_loop > begin_propagation_loop)
		{

#pragma omp parallel for
			for(int i = 0; i < num_trajectories; i++)
			{
				int thread_id = omp_get_thread_num();

				single_trajectory_dump_prop(
					&heads[thread_id],
					rnd_streams[i],
					&phi_global[i * N],
					&abs_rho_diag_omp[thread_id * N],
					&rho_omp[thread_id * N * N],
					num_trajectories,
					i,
					begin_propagation_loop,
					end_propagation_loop,
					thread_id
					);
			}

			printf("Period = %d\n", end_propagation_loop);
			begin_propagation_loop = end_propagation_loop;
			period = end_propagation_loop;

			dump_propagation_data(N, num_omp_threads, num_trajectories, "a", calc_characteristics, dump_rho);
			refresh_dump_data(N, num_omp_threads);
			refresh_aux_data(N);
		}
	}


	delete_aux_data();

	for(int i = 0; i < num_omp_threads; i++)
	{
		delete_split_struct_not_member (&heads[i]);
	}
	delete (heads);
	delete_split_struct (head);
	delete (rnd_streams);
	delete[] phi_global;
	delete[] abs_rho_diag_omp;
	delete[] rho_omp;
	delete[] ts;
	delete[] etas;
}


void omp_qj_statistic(char input_file_name[],
	char aux_file_name[],
	int num_periods,
	int num_dumps,
	int dump_type,
	int num_periods_in_trans_proc,
	int num_omp_threads,
	int num_trajectories,
	int rnd_max,
	int rnd_cur,
	int init_state_id,
	double mean_low_limit,
	double mean_high_limit,
	int dump_rho,
	int avg_dump,
	int dump_characteristics,
	int btw_jump_times,
	int borders_type,
	int stationary,
	int after_dump,
	int double_scale_dump,
	int deep_characteristic,
	int mc_specific,
	int mc_type
)
{
	printf("num_omp_threads: %d\n\n", num_omp_threads);
	omp_set_num_threads(num_omp_threads);

	FILE * input_file;
	input_file = fopen(input_file_name, "rb");
	split * head = create_struct_bin(input_file);
	fclose(input_file);

	int N = head->N;

	// #### Random initialization begin ####
	VSLStreamStatePtr * rnd_streams = new VSLStreamStatePtr[num_trajectories];
	vslNewStream (&rnd_streams[0], VSL_BRNG_MCG59, 777);
	for(int i = 1; i < num_trajectories; i++)
	{
		vslCopyStream(&rnd_streams[i], rnd_streams[0]);
	}
	for(int i = 0; i < num_trajectories; i++)
	{
		vslLeapfrogStream(rnd_streams[i], rnd_cur + i, rnd_max);
	}
	// #### Random initialization end ####

	split * heads = new split[num_omp_threads];
	for(int i = 0; i < num_omp_threads; i++)
	{
		cmp_struct_not_member(head, &heads[i]);
	}

	init_main_statistics_data(N, num_trajectories, stationary);
	init_dump_statistics_data(num_trajectories, btw_jump_times, stationary);
	//init_aux_statistics_data(aux_file_name, N, num_trajectories, borders_type, dump_characteristics, stationary);

	double step = double(num_periods) / double(num_dumps - 1);
	double start = 0.0;
	if (dump_type == 1)
	{
		start = 0.0;
		step = log10(double(num_periods)) / double(num_dumps - 1);
	}

#pragma omp parallel for
	for (int i = 0; i < num_trajectories; i++)
	{
		int thread_id = omp_get_thread_num();
		single_trajectory_statistics_init_prop(
			&heads[thread_id],
			rnd_streams[i],
			&phi_global[i * N],
			&abs_rho_diag_all[i * N],
			num_periods_in_trans_proc,
			i,
			init_state_id,
			mean_low_limit,
			mean_high_limit,
			dump_characteristics,
			&abs_rho_diag_all_stationary[i * N],
			stationary,
			borders_type
		);
	}

	for (int trajectory_id = 0; trajectory_id < num_trajectories; trajectory_id++)
	{
		ts[trajectory_id] = 0.0;
	}

	period = 1;
	printf("Period = %d\n", period);

	dump_statistics_data(N,	num_trajectories, "w", dump_rho, avg_dump, dump_characteristics, stationary);
	refresh_main_statistics_data(N, num_trajectories, dump_characteristics, stationary);

	double current_dump_limit = start;
	int begin_propagation_loop = period;
	int end_propagation_loop = 0;

	while (begin_propagation_loop != num_periods)
	{
		current_dump_limit += step;
		end_propagation_loop = num_periods;

		if (dump_type == 0)
		{
			if (current_dump_limit < double(end_propagation_loop))
			{
				end_propagation_loop = int(current_dump_limit);
			}
		}
		else if (dump_type == 1)
		{
			int exp_limit = int(pow(10.0, current_dump_limit + 1.0e-10));
			if (exp_limit < end_propagation_loop)
			{
				if (exp_limit <= begin_propagation_loop)
				{
					end_propagation_loop = begin_propagation_loop + 1;
				}
				else
				{
					end_propagation_loop = exp_limit;
				}
			}
		}

		if (end_propagation_loop > begin_propagation_loop)
		{

#pragma omp parallel for
			for(int i = 0; i < num_trajectories; i++)
			{
				int thread_id = omp_get_thread_num();

				single_trajectory_statistics_dump_prop(
					&heads[thread_id],
					rnd_streams[i],
					&phi_global[i * N],
					&abs_rho_diag_all[i * N],
					i,
					begin_propagation_loop,
					end_propagation_loop,
					dump_characteristics,
					btw_jump_times,
					&abs_rho_diag_all_stationary[i * N],
					stationary,
					borders_type,
					deep_characteristic,
					mc_specific,
					mc_type
				);
			}

			//printf("Period = %d\n", end_propagation_loop);
			begin_propagation_loop = end_propagation_loop;
			period = end_propagation_loop;

			dump_statistics_data(N,	num_trajectories, "a", dump_rho, avg_dump, dump_characteristics, stationary);
			refresh_main_statistics_data(N, num_trajectories, dump_characteristics, stationary);
		}
	}

	if (double_scale_dump > 0)
	{
		double current_dump_limit = start;
		int begin_propagation_loop = 0;
		int end_propagation_loop = 0;

		int period_shift = num_periods;

		while (begin_propagation_loop != num_periods)
		{
			current_dump_limit += step;
			end_propagation_loop = num_periods;

			if (dump_type == 0)
			{
				if (current_dump_limit < double(end_propagation_loop))
				{
					end_propagation_loop = int(current_dump_limit);
				}
			}
			else if (dump_type == 1)
			{
				int exp_limit = int(pow(10.0, current_dump_limit + 1.0e-10));
				if (exp_limit < end_propagation_loop)
				{
					if (exp_limit <= begin_propagation_loop)
					{
						end_propagation_loop = begin_propagation_loop + 1;
					}
					else
					{
						end_propagation_loop = exp_limit;
					}
				}
			}

			if (end_propagation_loop > begin_propagation_loop)
			{

#pragma omp parallel for
				for (int i = 0; i < num_trajectories; i++)
				{
					int thread_id = omp_get_thread_num();

					single_trajectory_statistics_dump_prop(
						&heads[thread_id],
						rnd_streams[i],
						&phi_global[i * N],
						&abs_rho_diag_all[i * N],
						i,
						begin_propagation_loop,
						end_propagation_loop,
						dump_characteristics,
						btw_jump_times,
						&abs_rho_diag_all_stationary[i * N],
						stationary,
						borders_type, 
						deep_characteristic,
						mc_specific,
						mc_type
					);
				}
			
				begin_propagation_loop = end_propagation_loop;
				period = end_propagation_loop + period_shift;
				//printf("Period = %d\n", period);

				dump_statistics_data(N, num_trajectories, "a", dump_rho, avg_dump, dump_characteristics, stationary);
				refresh_main_statistics_data(N, num_trajectories, dump_characteristics, stationary);
			}
		}
	}

	dump_aux_statistics_data(N, num_periods, num_trajectories, "w", avg_dump, dump_characteristics, btw_jump_times, stationary);

	if (after_dump > 0)
	{
		for (int trajectory_id = 0; trajectory_id < num_trajectories; trajectory_id++)
		{
			ts[trajectory_id] = 0.0;
		}

		double current_dump_limit = 1;
		int begin_propagation_loop = 1;
		int end_propagation_loop = after_dump;

		while (begin_propagation_loop != after_dump)
		{
			current_dump_limit += 1;
			end_propagation_loop = after_dump;

			if (current_dump_limit < double(end_propagation_loop))
			{
				end_propagation_loop = int(current_dump_limit);
			}

			if (end_propagation_loop > begin_propagation_loop)
			{

#pragma omp parallel for
				for (int i = 0; i < num_trajectories; i++)
				{
					int thread_id = omp_get_thread_num();

					single_trajectory_statistics_dump_prop(
						&heads[thread_id],
						rnd_streams[i],
						&phi_global[i * N],
						&abs_rho_diag_all[i * N],
						i,
						begin_propagation_loop,
						end_propagation_loop,
						dump_characteristics,
						btw_jump_times,
						&abs_rho_diag_all_stationary[i * N],
						stationary,
						borders_type,
						deep_characteristic,
						mc_specific,
						mc_type
					);
				}

				printf("aux dump period = %d\n", end_propagation_loop);
				begin_propagation_loop = end_propagation_loop;
				period = end_propagation_loop;

				dump_after_data(N, num_trajectories, "a", stationary);
				refresh_main_statistics_data(N, num_trajectories, dump_characteristics, stationary);
			}
		}

	}

	delete_main_statistics_data(dump_characteristics, stationary);
	delete_dump_statistics_data(dump_characteristics, btw_jump_times, stationary);

	//delete_aux_statistics_data(num_trajectories, dump_characteristics, stationary);

	for(int i = 0; i < num_omp_threads; i++)
	{
		delete_split_struct_not_member (&heads[i]);
	}
	delete (heads);
	delete_split_struct (head);
	delete (rnd_streams);
}


void seq_qj_hist(char input_file_name[],
				 char aux_file_name[],
				 int num_periods,
				 int num_dumps,
				 int num_periods_in_trans_proc,
				 int num_omp_threads,
				 int num_trajectories,
				 int rnd_max,
				 int rnd_cur,
				 int init_state_id,
				 int borders_type,
				 int num_e_intervals,
				 int num_mc_intervals,
				 double energy_min,
				 double energy_max
				 )
{
	printf("num_omp_threads: %d\n\n", num_omp_threads);
	omp_set_num_threads(num_omp_threads);

	FILE * input_file;
	input_file = fopen(input_file_name, "rb");
	split * head = create_struct_bin(input_file);
	fclose(input_file);

	int N = head->N;

	VSLStreamStatePtr * rnd_streams = new VSLStreamStatePtr[num_trajectories];
	vslNewStream (&rnd_streams[0], VSL_BRNG_MCG59, 777);
	for(int i = 1; i < num_trajectories; i++)
	{
		vslCopyStream(&rnd_streams[i], rnd_streams[0]);
	}
	for(int i = 0; i < num_trajectories; i++)
	{
		vslLeapfrogStream(rnd_streams[i], rnd_cur + i, rnd_max);
	}

	split * heads = new split[num_omp_threads];
	for(int i = 0; i < num_omp_threads; i++)
	{
		cmp_struct_not_member(head, &heads[i]);
	}

	init_main_hist_data(N, num_trajectories);
	init_aux_hist_data(aux_file_name, N);
	init_dump_hist_data(N, num_trajectories, num_e_intervals, num_mc_intervals, energy_min, energy_max);

	double step = double(num_periods) / double(num_dumps - 1);
	double start = 0.0;

#pragma omp parallel for
	for(int i = 0; i < num_trajectories; i++)
	{
		int thread_id = omp_get_thread_num();

		single_trajectory_hist_init_prop(
			&heads[thread_id],
			rnd_streams[i],
			&phi_global[i * N],
			&abs_rho_diag_all[i * N],
			num_periods_in_trans_proc,
			i,
			init_state_id);
	}

	period = 1;
	if (period % 1000 == 0)
	{
		printf("Period = %d\n", period);
	}
	
	refresh_main_hist_data(N, num_trajectories);

	double current_dump_limit = start;
	int begin_propagation_loop = period;
	int end_propagation_loop = 0;


	while (begin_propagation_loop != num_periods)
	{
		current_dump_limit += step;
		end_propagation_loop = num_periods;

		if (current_dump_limit < double(end_propagation_loop))
		{
			end_propagation_loop = int(current_dump_limit);
		}

		if (end_propagation_loop > begin_propagation_loop)
		{

#pragma omp parallel for
			for(int i = 0; i < num_trajectories; i++)
			{
				int thread_id = omp_get_thread_num();

				single_trajectory_hist_dump_prop(
					&heads[thread_id],
					rnd_streams[i],
					&phi_global[i * N],
					&abs_rho_diag_all[i * N],
					&sub_mult[i * N],
					&histogramm[i * num_e_intervals * num_mc_intervals],
					i,
					begin_propagation_loop,
					end_propagation_loop,
					num_e_intervals,
					num_mc_intervals,
					borders_type);
			}

			if (period % 1000 == 0)
			{
				printf("Period = %d\n", period);
			}
			begin_propagation_loop = end_propagation_loop;
			period = end_propagation_loop;

			refresh_main_hist_data(N, num_trajectories);
		}
	}

	dump_hist_data(num_trajectories, "w", num_e_intervals, num_mc_intervals);

	delete_main_hist_data();
	delete_aux_hist_data();
	delete_dump_hist_data();

	for(int i = 0; i < num_omp_threads; i++)
	{
		delete_split_struct_not_member (&heads[i]);
	}
	delete (heads);
	delete_split_struct (head);
	delete (rnd_streams);
}