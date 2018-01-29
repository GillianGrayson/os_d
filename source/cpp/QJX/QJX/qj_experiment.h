#pragma once
#include "config.h"
#include "data.h"
#include "qj_data.h"
#include "dump.h"

class QJExperimentBehavior
{
public:
	virtual void trans_process(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd) const = 0;
	virtual void obs_process(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd) const = 0;
};

class LpnExperimentBehaviour : public QJExperimentBehavior
{
public:
	virtual void trans_process(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd) const;
	virtual void obs_process(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd) const;
};

class StdExperimentBehaviour : public QJExperimentBehavior
{
public:
	virtual void trans_process(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd) const;
	virtual void obs_process(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd) const;
};

class CorrDimExperimentBehaviour : public QJExperimentBehavior
{
public:
	virtual void trans_process(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd) const;
	virtual void obs_process(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd) const;
};

class SigmaExperimentBehaviour : public QJExperimentBehavior
{
public:
	virtual void trans_process(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd) const;
	virtual void obs_process(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd) const;
};

void prop_step(MKL_Complex16 * phi, MKL_Complex16 * matrix, MKL_Complex16 * res, int sys_size);

MKL_Complex16 mult_scalar_double(MKL_Complex16 a, double b);

MKL_Complex16 mult_scalar_complex(MKL_Complex16 * a, MKL_Complex16 * b, int N);

int is_norm_crossed(MKL_Complex16 * phi, double * eta, int sys_size);

double norm_square(MKL_Complex16 * phi, int sys_size);

void recovery(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, Split * head, int tr_id);

void one_period_branch(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, Split * head, int tr_id, Split * branch);

void one_sub_period_cd(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int tr_id, int part_id, int thread_id);

void one_period_cd_tp(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int tr_id, int thread_id);

void one_period_sigma_obs(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int tr_id, int thread_id, int period_id);

void one_period_cd_obs(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int tr_id, int thread_id, int period_id);

void one_period(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int tr_id, int thread_id);

double get_norm_cd(double * vec, int size);

double get_mean_simple(double * adr, int sys_size);

double get_dispersion_simple(double mean_curr, double mean_start);

double get_m2(double * adr, int sys_size, double mean);

double get_energy(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int tr_id);

void calc_chars_start_std(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int tr_id);

void calc_chars_std(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int tr_id);

void evo_chars_std(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int tr_id, int dump_id);

void calc_chars_start_lpn(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int tr_id);

void calc_chars_lpn(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int tr_id);

void calc_ci(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int tr_id);

void evo_chars_lpn(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int tr_id, int dump_id);

void resresh_times(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int tr_id);

void copy_trajectory_lpn(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int tr_id);

void copy_trajectory_data(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int tr_id);

void var_trajectory_lpn(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int tr_id);

void lambda_lpn(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int tr_id);

void trans_process_single_std(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int tr_id, int thread_id);

void trans_process_single_cd(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int tr_id, int thread_id);