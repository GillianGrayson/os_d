#pragma once
#include "config.h"
#include "data.h"
#include "qj_data.h"

class QJFreeBehavior
{
public:
	virtual void free_data(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd) const = 0;
};

class LpnFreeBehaviour : public QJFreeBehavior
{
public:
	virtual void free_data(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd) const;
};

class StdFreeBehaviour : public QJFreeBehavior
{
public:
	virtual void free_data(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd) const;
};

class CorrDimFreeBehaviour : public QJFreeBehavior
{
public:
	virtual void free_data(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd) const;
};

void free_splits_deep(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd);

void free_splits(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd);

void free_streams(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd);

void free_streams_var(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd);

void free_basic_data(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd);

void free_dump_priods(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd);

void free_obs_std(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd);

void free_obs_lpn(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd);

void free_obs_cd(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd);