#include "mkl.h"
#include "stdio.h"
#include "math.h"
#include "string.h"


struct split
{
	bool type;
	split * next;
	split * prev;
	double dt;
	unsigned int steps;
	unsigned int counter;
	unsigned int N;
	MKL_Complex16 * matrix;
	double * g;
};

split * create_struct_bin (FILE * file);
split * create_struct (FILE * file);

void delete_split_struct (split * head);

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
);


double norm_vector2(MKL_Complex16 * phi, int N);

void cmp_struct_not_member(split * head1, split * head2);
void delete_split_struct_not_member (split * head);

void set_init_conditions(MKL_Complex16 * phi, int N, int state_id);

void omp_qj(char input_file_name[],
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
							  int dump_rho);

void omp_qj_statistic(
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
	int mc_type,
	int num_att_trajectories,
	double var_eps
);


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
				 );

