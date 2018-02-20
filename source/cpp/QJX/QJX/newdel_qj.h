#pragma once
#include "config.h"
#include "data.h"

class QJNewDelBehavior
{
public:
	virtual void init_data(AllData * ad) const = 0;
	virtual void free_data(AllData * ad) const = 0;
};

class LpnNewDelBehaviour : public QJNewDelBehavior
{
public:
	virtual void init_data(AllData * ad) const;
	virtual void free_data(AllData * ad) const;
};

class StdNewDelBehaviour : public QJNewDelBehavior
{
public:
	virtual void init_data(AllData * ad) const;
	virtual void free_data(AllData * ad) const;
};

class CorrDimNewDelBehaviour : public QJNewDelBehavior
{
public:
	virtual void init_data(AllData * ad) const;
	virtual void free_data(AllData * ad) const;
};

class SigmaNewDelBehaviour : public QJNewDelBehavior
{
public:
	virtual void init_data(AllData * ad) const;
	virtual void free_data(AllData * ad) const;
};

void init_splits_cd(AllData * ad);
void init_splits(AllData * ad);
void init_streams(AllData * ad);
void leap_frog_single_stream(AllData * ad, int tr_id);
void leap_frog_all_streams(AllData * ad);
void copy_streams(AllData * ad);
void init_streams_var(AllData * ad);
void init_basic_data(AllData * ad);
void init_dump_periods_cd_deep(AllData * ad);
void init_dump_periods(AllData * ad);
void init_obs_std(AllData * ad);
void init_obs_lpn(AllData * ad);
void init_obs_cd(AllData * ad);
void init_start_state(AllData * ad, int tr_id);

void free_splits_deep(AllData * ad);
void free_splits(AllData * ad);
void free_streams(AllData * ad);
void free_streams_var(AllData * ad);
void free_basic_data(AllData * ad);
void free_dump_priods(AllData * ad);
void free_obs_std(AllData * ad);
void free_obs_lpn(AllData * ad);
void free_obs_cd(AllData * ad);


