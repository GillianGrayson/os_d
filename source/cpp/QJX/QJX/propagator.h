#pragma once
#include "config.h"
#include "data.h"

class PropagateBehavior
{
public:

	virtual void init_prop_data(AllData * ad) const = 0;
	virtual void free_prop_data(AllData * ad) const = 0;
	virtual void init_prop_data_cd(AllData * ad) const = 0;
	virtual void free_prop_data_cd(AllData * ad) const = 0;

	virtual void one_period(AllData * ad, int tr_id, int thread_id, int period_id) const = 0;
	virtual void one_period_cd_tp(AllData * ad, int tr_id, int thread_id, int period_id) const = 0;
	virtual void one_period_cd_obs(AllData * ad, int tr_id, int thread_id, int period_id) const = 0;
	virtual void one_period_sigma_obs(AllData * ad, int tr_id, int thread_id, int period_id) const = 0;
};

class QJPropagateBehavior : public PropagateBehavior
{
public:

	virtual void init_prop_data(AllData * ad) const;
	virtual void free_prop_data(AllData * ad) const;
	virtual void init_prop_data_cd(AllData * ad) const;
	virtual void free_prop_data_cd(AllData * ad) const;

	virtual void one_period(AllData * ad, int tr_id, int thread_id, int period_id) const;
	virtual void one_period_cd_tp(AllData * ad, int tr_id, int thread_id, int period_id) const;
	virtual void one_period_cd_obs(AllData * ad, int tr_id, int thread_id, int period_id) const;
	virtual void one_period_sigma_obs(AllData * ad, int tr_id, int thread_id, int period_id) const;
};

class RKPropagateBehavior : public PropagateBehavior
{
public:

	virtual void init_prop_data(AllData * ad) const;
	virtual void free_prop_data(AllData * ad) const;
	virtual void init_prop_data_cd(AllData * ad) const;
	virtual void free_prop_data_cd(AllData * ad) const;

	virtual void one_period(AllData * ad, int tr_id, int thread_id, int period_id) const;
	virtual void one_period_cd_tp(AllData * ad, int tr_id, int thread_id, int period_id) const;
	virtual void one_period_cd_obs(AllData * ad, int tr_id, int thread_id, int period_id) const;
	virtual void one_period_sigma_obs(AllData * ad, int tr_id, int thread_id, int period_id) const;
};


