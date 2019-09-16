#pragma once
#include "config.h"
#include "data.h"
#include "dump.h"
#include "propagator.h"
#include "core.h"

class ExperimentBehavior
{
public:
	virtual void trans_process(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const = 0;
	virtual void obser_process(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const = 0;
};

class LpnExperimentBehaviour : public ExperimentBehavior
{
public:
	virtual void trans_process(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const;
	virtual void obser_process(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const;
};

class LpnMultExperimentBehaviour: public ExperimentBehavior
{
public:
	virtual void trans_process(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const;
	virtual void obser_process(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const;
};

class StdExperimentBehaviour : public ExperimentBehavior
{
public:
	virtual void trans_process(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const;
	virtual void obser_process(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const;
};

class CorrDimExperimentBehaviour : public ExperimentBehavior
{
public:
	virtual void trans_process(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const;
	virtual void obser_process(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const;
};

class SigmaExperimentBehaviour : public ExperimentBehavior
{
public:
	virtual void trans_process(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const;
	virtual void obser_process(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const;
};

class StdDeepExperimentBehaviour : public ExperimentBehavior
{
public:
	virtual void trans_process(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const;
	virtual void obser_process(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const;
};

class LpnDeepExperimentBehaviour : public ExperimentBehavior
{
public:
	virtual void trans_process(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const;
	virtual void obser_process(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const;
};

class LpnDeepPer1TExperimentBehaviour : public ExperimentBehavior
{
public:
	virtual void trans_process(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const;
	virtual void obser_process(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const;
};

class LpnAllExperimentBehaviour : public ExperimentBehavior
{
public:
	virtual void trans_process(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const;
	virtual void obser_process(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const;
};

MKL_Complex16 mult_scalar_double(MKL_Complex16 a, double b);

MKL_Complex16 mult_scalar_complex(MKL_Complex16 * a, MKL_Complex16 * b, int N);

int is_norm_crossed(MKL_Complex16 * phi, double * eta, int sys_size);

double norm_square(MKL_Complex16 * phi, int sys_size);

void recovery(AllData * ad, Split * head, int tr_id);

double get_norm_cd(double * vec, int size);

double get_mean_simple(double * adr, int sys_size);

double get_dispersion_simple(double mean_curr, double mean_start);

double get_m2(double * adr, int sys_size, double mean);

double get_energy(AllData * ad, int tr_id);

MKL_Complex16 get_spec_jcs(AllData * ad, int tr_id);
MKL_Complex16 get_spec_ps(AllData * ad, int tr_id);
MKL_Complex16 get_spec_2_ps(AllData * ad, int tr_id);
MKL_Complex16 get_spec_3_ps(AllData * ad, int tr_id);

MKL_Complex16 get_spec_mbl(AllData * ad, int tr_id);

MKL_Complex16 get_num_photons_jcs(AllData * ad, int tr_id);
MKL_Complex16 get_num_photons_ps(AllData * ad, int tr_id);

void resresh_times(AllData * ad, int tr_id);

void copy_trajectory_lpn(AllData * ad, int tr_id, int base_tr_id);

void copy_stream_lpn(AllData * ad, int tr_id, int base_tr_id);

void copy_trajectory_data(AllData * ad, int tr_id, int base_tr_id);

void var_trajectory_lpn(AllData * ad, CoreBehavior * cb, int tr_id, int base_tr_id);

void var_first(
	MKL_Complex16 * phi_var,
	double * phi_var_double,
	MKL_Complex16 * phi_var_all,
	MKL_Complex16 * scalar_mults_all,
	AllData * ad,
	CoreBehavior *cb
);

void var_first_with_history(
	MKL_Complex16 * phi_var,
	double * phi_var_double,
	MKL_Complex16 * phi_var_all,
	MKL_Complex16 * scalar_mults_all,
	AllData * ad,
	CoreBehavior *cb
);

void var_not_first(
	int tr_id,
	MKL_Complex16 * phi_var,
	double * phi_var_double,
	MKL_Complex16 * phi_var_all,
	MKL_Complex16 * scalar_mults_all,
	AllData * ad,
	CoreBehavior *cb
);

void var_not_first_with_history(
	int tr_id,
	MKL_Complex16 * phi_var,
	double * phi_var_double,
	MKL_Complex16 * phi_var_all,
	MKL_Complex16 * scalar_mults_all,
	AllData * ad,
	CoreBehavior *cb
);

void gs_orth_init(AllData * ad, CoreBehavior *cb);

void only_orth(AllData * ad, CoreBehavior *cb, MKL_Complex16 * phi_var_all);

void gs_orth_evo(AllData * ad, CoreBehavior *cb, MKL_Complex16 *phi_var_all);

void lambda_lpn(AllData * ad, CoreBehavior *cb, int tr_id, int base_tr_id);

void lambda_lpn_per_periods(AllData * ad, CoreBehavior *cb, int tr_id, int base_tr_id, int num_steps_T, int curr_step, int num_periods);

void lambda_lpn_all(AllData * ad, CoreBehavior *cb);

void trans_process_single(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb, int tr_id, int thread_id);

void trans_process_single_deep(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb, int tr_id, int thread_id);