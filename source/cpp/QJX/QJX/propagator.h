#pragma once
#include "config.h"
#include "data.h"
#include "core.h"

class PropagateBehavior
{
public:

	virtual void init_prop_data(AllData * ad, CoreBehavior * cb) const = 0;
	virtual void free_prop_data(AllData * ad, CoreBehavior * cb) const = 0;
	virtual void init_prop_data_deep(AllData * ad, CoreBehavior * cb) const = 0;
	virtual void free_prop_data_deep(AllData * ad, CoreBehavior * cb) const = 0;

	virtual void one_period(AllData * ad, CoreBehavior * cb, int tr_id, int th_id, int period_id) const = 0;
	virtual void one_period_trp_deep(AllData * ad, CoreBehavior * cb, int tr_id, int th_id, int period_id) const = 0;
	virtual void one_period_obs_deep(AllData * ad, CoreBehavior * cb, int tr_id, int th_id, int period_id) const = 0;
	virtual void one_period_obs_deep_lpn(AllData * ad, CoreBehavior * cb, int tr_id, int th_id, int period_id) const = 0;
	virtual void one_period_obs_deep_cd(AllData * ad, CoreBehavior * cb, int tr_id, int th_id, int period_id) const = 0;
	virtual void one_period_obs_deep_sigma(AllData * ad, CoreBehavior * cb, int tr_id, int th_id, int period_id) const = 0;
};

class QJPropagateBehavior : public PropagateBehavior
{
public:

	virtual void init_prop_data(AllData * ad, CoreBehavior * cb) const;
	virtual void free_prop_data(AllData * ad, CoreBehavior * cb) const;
	virtual void init_prop_data_deep(AllData * ad, CoreBehavior * cb) const;
	virtual void free_prop_data_deep(AllData * ad, CoreBehavior * cb) const;

	virtual void one_period(AllData * ad, CoreBehavior * cb, int tr_id, int th_id, int period_id) const;
	virtual void one_period_trp_deep(AllData * ad, CoreBehavior * cb, int tr_id, int th_id, int period_id) const;
	virtual void one_period_obs_deep(AllData * ad, CoreBehavior * cb, int tr_id, int th_id, int period_id) const;
	virtual void one_period_obs_deep_lpn(AllData * ad, CoreBehavior * cb, int tr_id, int th_id, int period_id) const;
	virtual void one_period_obs_deep_cd(AllData * ad, CoreBehavior * cb, int tr_id, int th_id, int period_id) const;
	virtual void one_period_obs_deep_sigma(AllData * ad, CoreBehavior * cb, int tr_id, int th_id, int period_id) const;
};

class RKPropagateBehavior : public PropagateBehavior
{
public:

	virtual void init_prop_data(AllData * ad, CoreBehavior * cb) const;
	virtual void free_prop_data(AllData * ad, CoreBehavior * cb) const;
	virtual void init_prop_data_deep(AllData * ad, CoreBehavior * cb) const;
	virtual void free_prop_data_deep(AllData * ad, CoreBehavior * cb) const;

	virtual void one_period(AllData * ad, CoreBehavior * cb, int tr_id, int th_id, int period_id) const;
	virtual void one_period_trp_deep(AllData * ad, CoreBehavior * cb, int tr_id, int th_id, int period_id) const;
	virtual void one_period_obs_deep(AllData * ad, CoreBehavior * cb, int tr_id, int th_id, int period_id) const;
	virtual void one_period_obs_deep_lpn(AllData * ad, CoreBehavior * cb, int tr_id, int th_id, int period_id) const;
	virtual void one_period_obs_deep_cd(AllData * ad, CoreBehavior * cb, int tr_id, int th_id, int period_id) const;
	virtual void one_period_obs_deep_sigma(AllData * ad, CoreBehavior * cb, int tr_id, int th_id, int period_id) const;
};


