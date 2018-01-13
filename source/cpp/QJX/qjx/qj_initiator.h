#pragma once
#include "config.h"
#include "data.h"
#include "qj_data.h"

class QJInitBehavior
{
public:
	virtual void init_data(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd) const = 0;
};

class LyapunovMCInitBehaviour : public QJInitBehavior
{
public:
	virtual void init_data(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd) const;
};

void init_splits(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd);

void init_streams(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd);

void init_streams_var(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd);

void init_basic_data(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd);

void init_dump_periods(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd);

void init_obs_std(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd);

void init_obs_lpn(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd);

void init_start_state(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int tr_id);


