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

MKL_Complex16 get_spec(AllData * ad, int tr_id);

void resresh_times(AllData * ad, int tr_id);

void copy_trajectory_lpn(AllData * ad, int tr_id);

void copy_stream_lpn(AllData * ad, int tr_id);

void copy_trajectory_data(AllData * ad, int tr_id);

void var_trajectory_lpn(AllData * ad, CoreBehavior * cb, int tr_id);

void lambda_lpn(AllData * ad, CoreBehavior *cb, int tr_id);

void trans_process_single(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb, int tr_id, int thread_id);

void trans_process_single_deep(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb, int tr_id, int thread_id);