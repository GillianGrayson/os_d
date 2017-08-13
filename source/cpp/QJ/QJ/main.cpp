#define _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_DEPRECATE
#include "header.h"
#include "omp.h"
#include <iostream>

int main(int argc, char **argv)
{
	if(argc < 3)
	{
		printf("Not enough arguments\n");
		return -2;
	}

	int num_periods;
	int init_state_id;
	int propagation_type;
	int num_dumps;
	int dump_type;
	int num_periods_in_trans_proc;
	int num_omp_threads;
	int num_trajectories;
	int rnd_max;
	int rnd_cur;
	int calc_characteristics;
	int dump_rho;
	double mean_low_limit;
	double mean_high_limit;
	int avg_dump;
	int dump_characteristics;
	int btw_jump_times;
	int borders_type;
	int stationary;
	int num_e_intervals;
	int num_mc_intervals;
	double energy_min;
	double energy_max;
	int after_dump;
	int double_scale_dump;
	int deep_characteristic;
	int mc_specific;
	int mc_type;
	
	FILE * config_file = fopen("config.txt", "r");
	fscanf(config_file, "num_periods = %d\n", &num_periods);
	fscanf(config_file, "init_state_id = %d\n", &init_state_id);
	fscanf(config_file, "propagation_type = %d\n", &propagation_type);
	fscanf(config_file, "num_dumps = %d\n", &num_dumps);
	fscanf(config_file, "dump_type = %d\n", &dump_type);
	fscanf(config_file, "num_periods_in_trans_proc = %d\n", &num_periods_in_trans_proc);
	fscanf(config_file, "num_omp_threads = %d\n", &num_omp_threads);
	fscanf(config_file, "num_trajectories = %d\n", &num_trajectories);
	fscanf(config_file, "rnd_max = %d\n", &rnd_max);
	fscanf(config_file, "rnd_cur = %d\n", &rnd_cur);
	fscanf(config_file, "calc_characteristics = %d\n", &calc_characteristics);
	fscanf(config_file, "dump_rho = %d\n", &dump_rho);
	fscanf(config_file, "mean_low_limit = %lf\n", &mean_low_limit);
	fscanf(config_file, "mean_high_limit = %lf\n", &mean_high_limit);
	fscanf(config_file, "avg_dump = %d\n", &avg_dump);
	fscanf(config_file, "dump_characteristics = %d\n", &dump_characteristics);
	fscanf(config_file, "btw_jump_times = %d\n", &btw_jump_times);
	fscanf(config_file, "borders_type = %d\n", &borders_type);
	fscanf(config_file, "stationary = %d\n", &stationary);
	fscanf(config_file, "num_e_intervals = %d\n", &num_e_intervals);
	fscanf(config_file, "num_mc_intervals = %d\n", &num_mc_intervals);
	fscanf(config_file, "energy_min = %lf\n", &energy_min);
	fscanf(config_file, "energy_max = %lf\n", &energy_max);
	fscanf(config_file, "after_dump = %d\n", &after_dump);
	fscanf(config_file, "double_scale_dump = %d\n", &double_scale_dump);
	fscanf(config_file, "deep_characteristic = %d\n", &deep_characteristic);
	fscanf(config_file, "mc_specific = %d\n", &mc_specific);
	fscanf(config_file, "mc_type = %d\n", &mc_type);
	fclose(config_file);

	std::cout << argv[1] << std::endl;
	std::cout << argv[2] << std::endl;

	double start_time = omp_get_wtime();

	if (propagation_type == 0)
	{
		omp_qj(
			argv[1],
			argv[2],
			num_periods,
			num_dumps,
			dump_type,
			num_periods_in_trans_proc,
			num_omp_threads,
			num_trajectories,
			rnd_max, rnd_cur,
			init_state_id,
			calc_characteristics,
			dump_rho);
	}
	else if (propagation_type == 1)
	{
		omp_qj_statistic(
			argv[1],
			argv[2],
			num_periods,
			num_dumps,
			dump_type,
			num_periods_in_trans_proc,
			num_omp_threads,
			num_trajectories,
			rnd_max,
			rnd_cur,
			init_state_id,
			mean_low_limit,
			mean_high_limit,
			dump_rho,
			avg_dump,
			dump_characteristics,
			btw_jump_times,
			borders_type,
			stationary,
			after_dump,
			double_scale_dump,
			deep_characteristic,
			mc_specific,
			mc_type);
	}
	else if (propagation_type == 2)
	{
		seq_qj_hist(
			argv[1],
			argv[2],
			num_periods,
			num_dumps,
			num_periods_in_trans_proc,
			num_omp_threads,
			num_trajectories,
			rnd_max,
			rnd_cur,
			init_state_id,
			borders_type,
			num_e_intervals,
			num_mc_intervals,
			energy_min,
			energy_max
			);
	}
	
	double time = omp_get_wtime() - start_time;

	printf("elapsed_time: %0.3le\n", time);

	return 0;
}