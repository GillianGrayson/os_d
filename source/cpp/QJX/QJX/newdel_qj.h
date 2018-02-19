#pragma once
#include "config.h"
#include "data.h"
#include "data_qj.h"

class QJNewDelBehavior
{
public:
	virtual void init_data(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd) const = 0;
	virtual void free_data(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd) const = 0;
};

class LpnNewDelBehaviour : public QJNewDelBehavior
{
public:
	virtual void init_data(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd) const;
	virtual void free_data(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd) const;
};

class StdNewDelBehaviour : public QJNewDelBehavior
{
public:
	virtual void init_data(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd) const;
	virtual void free_data(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd) const;
};

class CorrDimNewDelBehaviour : public QJNewDelBehavior
{
public:
	virtual void init_data(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd) const;
	virtual void free_data(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd) const;
};

class SigmaNewDelBehaviour : public QJNewDelBehavior
{
public:
	virtual void init_data(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd) const;
	virtual void free_data(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd) const;
};

void init_splits_cd(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd);
void init_splits(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd);
void init_streams(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd);
void leap_frog_single_stream(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int tr_id);
void leap_frog_all_streams(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd);
void copy_streams(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd);
void init_streams_var(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd);
void init_basic_data(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd);
void init_dump_periods_cd_deep(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd);
void init_dump_periods(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd);
void init_obs_std(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd);
void init_obs_lpn(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd);
void init_obs_cd(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd);
void init_start_state(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int tr_id);

void free_splits_deep(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd);
void free_splits(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd);
void free_streams(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd);
void free_streams_var(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd);
void free_basic_data(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd);
void free_dump_priods(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd);
void free_obs_std(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd);
void free_obs_lpn(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd);
void free_obs_cd(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd);


