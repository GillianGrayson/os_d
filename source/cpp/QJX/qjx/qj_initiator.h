#pragma once
#include "config.h"
#include "data.h"
#include "qj_data.h"

class QJInitBehavior
{
public:
	virtual void init_data(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd) const = 0;
};

class LyapunovMCBehaviour : public QJInitBehavior
{
public:
	virtual void init_data(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd) const;
};
