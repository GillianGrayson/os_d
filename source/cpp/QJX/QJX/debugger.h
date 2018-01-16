#pragma once
#include "config.h"
#include "data.h"

class DebugBehavior
{
public:
	virtual void save(RunParam * rp, ConfigParam * cp, MainData * md) const = 0;
};

class DimerDebugBehaviour : public DebugBehavior
{
public:
	virtual void save(RunParam * rp, ConfigParam * cp, MainData * md) const;
};