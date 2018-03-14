#pragma once
#include "config.h"
#include "data.h"

class PropagateBehavior
{
public:

	virtual void init_prop_data(AllData * ad) const = 0;
	virtual void free_prop_data(AllData * ad) const = 0;
	virtual void init_prop_data_deep(AllData * ad) const = 0;
	virtual void free_prop_data_deep(AllData * ad) const = 0;

	virtual void one_period(AllData * ad, int tr_id, int thread_id, int period_id) const = 0;
	virtual void one_period_tp_deep(AllData * ad, int tr_id, int thread_id, int period_id) const = 0;
	virtual void one_period_obs_deep(AllData * ad, int tr_id, int thread_id, int period_id) const = 0;
	virtual void one_period_obs_cd(AllData * ad, int tr_id, int thread_id, int period_id) const = 0;
	virtual void one_period_obs_sigma(AllData * ad, int tr_id, int thread_id, int period_id) const = 0;
};

class QJPropagateBehavior : public PropagateBehavior
{
public:

	virtual void init_prop_data(AllData * ad) const;
	virtual void free_prop_data(AllData * ad) const;
	virtual void init_prop_data_deep(AllData * ad) const;
	virtual void free_prop_data_deep(AllData * ad) const;

	virtual void one_period(AllData * ad, int tr_id, int thread_id, int period_id) const;
	virtual void one_period_tp_deep(AllData * ad, int tr_id, int thread_id, int period_id) const;
	virtual void one_period_obs_deep(AllData * ad, int tr_id, int thread_id, int period_id) const;
	virtual void one_period_obs_cd(AllData * ad, int tr_id, int thread_id, int period_id) const;
	virtual void one_period_obs_sigma(AllData * ad, int tr_id, int thread_id, int period_id) const;
};

class RKPropagateBehavior : public PropagateBehavior
{
public:

	virtual void init_prop_data(AllData * ad) const;
	virtual void free_prop_data(AllData * ad) const;
	virtual void init_prop_data_deep(AllData * ad) const;
	virtual void free_prop_data_deep(AllData * ad) const;

	virtual void one_period(AllData * ad, int tr_id, int thread_id, int period_id) const;
	virtual void one_period_tp_deep(AllData * ad, int tr_id, int thread_id, int period_id) const;
	virtual void one_period_obs_deep(AllData * ad, int tr_id, int thread_id, int period_id) const;
	virtual void one_period_obs_cd(AllData * ad, int tr_id, int thread_id, int period_id) const;
	virtual void one_period_obs_sigma(AllData * ad, int tr_id, int thread_id, int period_id) const;
};


