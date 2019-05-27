#pragma once
#include "config.h"
#include "data.h"
#include "propagator.h"
#include "core.h"

class ExpNewDelBehavior
{
public:
	virtual void init_data(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const = 0;
	virtual void free_data(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const = 0;
};

class LpnNewDelBehaviour : public ExpNewDelBehavior
{
public:
	virtual void init_data(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const;
	virtual void free_data(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const;
};

class LpnMultNewDelBehaviour : public ExpNewDelBehavior
{
public:
	virtual void init_data(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const;
	virtual void free_data(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const;
};

class StdNewDelBehaviour : public ExpNewDelBehavior
{
public:
	virtual void init_data(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const;
	virtual void free_data(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const;
};

class CorrDimNewDelBehaviour : public ExpNewDelBehavior
{
public:
	virtual void init_data(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const;
	virtual void free_data(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const;
};

class SigmaNewDelBehaviour : public ExpNewDelBehavior
{
public:
	virtual void init_data(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const;
	virtual void free_data(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const;
};

class StdDeepNewDelBehaviour : public ExpNewDelBehavior
{
public:
	virtual void init_data(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const;
	virtual void free_data(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const;
};

class LpnDeepNewDelBehaviour : public ExpNewDelBehavior
{
public:
	virtual void init_data(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const;
	virtual void free_data(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const;
};

void init_streams(AllData * ad);
void leap_frog_single_stream(AllData * ad, int tr_id);
void leap_frog_all_streams(AllData * ad);
void copy_half_streams(AllData * ad);
void copy_streams(AllData * ad);
void init_streams_var(AllData * ad);
void init_basic_data(AllData * ad);
void init_dump_periods(AllData * ad);
void init_dump_periods_deep(AllData * ad);
void init_obs_std(AllData * ad);
void init_obs_lpn(AllData * ad);
void init_obs_cd(AllData * ad);
void init_start_state(AllData * ad, int tr_id);

void free_streams(AllData * ad);
void free_streams_var(AllData * ad);
void free_basic_data(AllData * ad);
void free_dump_priods(AllData * ad);
void free_obs_std(AllData * ad);
void free_obs_lpn(AllData * ad);
void free_obs_cd(AllData * ad);


