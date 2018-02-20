#pragma once
#include "config.h"
#include "data.h"
#include "experiment.h"

class PropagateBehavior
{
public:
	virtual void one_period(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int tr_id, int thread_id) const = 0;

	virtual void one_period_cd_tp(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int tr_id, int thread_id) const = 0;
	
	virtual void one_period_cd_obs(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int tr_id, int thread_id, int period_id) const = 0;

	virtual void one_period_sigma_obs(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int tr_id, int thread_id, int period_id) const = 0;
};

class QJPropagateBehavior : public PropagateBehavior
{
public:
	virtual void one_period(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int tr_id, int thread_id) const;

	virtual void one_period_cd_tp(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int tr_id, int thread_id) const;
	
	virtual void one_period_cd_obs(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int tr_id, int thread_id, int period_id) const;

	virtual void one_period_sigma_obs(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int tr_id, int thread_id, int period_id) const;
};

