#pragma once
#include "config.h"
#include "data.h"

class InitBehavior
{
public:
	virtual void init_sizes(RunParam * rp, ConfigParam * cp, MainData * md) const = 0;
	virtual void init_hamiltonians(RunParam * rp, ConfigParam * cp, MainData * md) const = 0;
	virtual void init_dissipators(RunParam * rp, ConfigParam * cp, MainData * md) const = 0;
	virtual void init_hamiltonians_qj(RunParam * rp, ConfigParam * cp, MainData * md) const = 0;
};

class DimerInitBehaviour : public InitBehavior
{
public:
	virtual void init_sizes(RunParam * rp, ConfigParam * cp, MainData * md) const;
	virtual void init_hamiltonians(RunParam * rp, ConfigParam * cp, MainData * md) const;
	virtual void init_dissipators(RunParam * rp, ConfigParam * cp, MainData * md) const;
	virtual void init_hamiltonians_qj(RunParam * rp, ConfigParam * cp, MainData * md) const;
};
