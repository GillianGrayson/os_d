#pragma once
#include "config.h"
#include "data.h"
#include "dump.h"

class ExperimentBehavior
{
public:
	virtual void trans_process(AllData * ad) const = 0;
	virtual void obser_process(AllData * ad) const = 0;
};

class LpnExperimentBehaviour : public ExperimentBehavior
{
public:
	virtual void trans_process(AllData * ad) const;
	virtual void obser_process(AllData * ad) const;
};

class StdExperimentBehaviour : public ExperimentBehavior
{
public:
	virtual void trans_process(AllData * ad) const;
	virtual void obser_process(AllData * ad) const;
};

class CorrDimExperimentBehaviour : public ExperimentBehavior
{
public:
	virtual void trans_process(AllData * ad) const;
	virtual void obser_process(AllData * ad) const;
};

class SigmaExperimentBehaviour : public ExperimentBehavior
{
public:
	virtual void trans_process(AllData * ad) const;
	virtual void obser_process(AllData * ad) const;
};

void prop_step(MKL_Complex16 * phi, MKL_Complex16 * matrix, MKL_Complex16 * res, int sys_size);

MKL_Complex16 mult_scalar_double(MKL_Complex16 a, double b);

MKL_Complex16 mult_scalar_complex(MKL_Complex16 * a, MKL_Complex16 * b, int N);

int is_norm_crossed(MKL_Complex16 * phi, double * eta, int sys_size);

double norm_square(MKL_Complex16 * phi, int sys_size);

void recovery(AllData * ad, Split * head, int tr_id);

void one_period_branch(AllData * ad, Split * head, int tr_id, Split * branch);

void one_sub_period_cd(AllData * ad, int tr_id, int part_id, int thread_id);

double get_norm_cd(double * vec, int size);

double get_mean_simple(double * adr, int sys_size);

double get_dispersion_simple(double mean_curr, double mean_start);

double get_m2(double * adr, int sys_size, double mean);

double get_energy(AllData * ad, int tr_id);

void calc_chars_start_std(AllData * ad, int tr_id);

void calc_chars_std(AllData * ad, int tr_id);

void evo_chars_std(AllData * ad, int tr_id, int dump_id);

void calc_chars_start_lpn(AllData * ad, int tr_id);

void calc_chars_lpn(AllData * ad, int tr_id);

void calc_ci(AllData * ad, int tr_id);

void evo_chars_lpn(AllData * ad, int tr_id, int dump_id);

void resresh_times(AllData * ad, int tr_id);

void copy_trajectory_lpn(AllData * ad, int tr_id);

void copy_trajectory_data(AllData * ad, int tr_id);

void var_trajectory_lpn(AllData * ad, int tr_id);

void lambda_lpn(AllData * ad, int tr_id);

void trans_process_single_std(AllData * ad, int tr_id, int thread_id);

void trans_process_single_cd(AllData * ad, int tr_id, int thread_id);