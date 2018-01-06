#pragma once
#include "config.h"
#include "data.h"

class FreeBehavior
{
public:
	virtual void free_hamiltonians(RunParam * rp, ConfigParam * cp, MainData * md) const = 0;
	virtual void free_dissipators(RunParam * rp, ConfigParam * cp, MainData * md) const = 0;
	virtual void free_hamiltonians_qj(RunParam * rp, ConfigParam * cp, MainData * md) const = 0;
};

class DimerFreeBehaviour : public FreeBehavior
{
public:
	virtual void free_hamiltonians(RunParam * rp, ConfigParam * cp, MainData * md) const;
	virtual void free_dissipators(RunParam * rp, ConfigParam * cp, MainData * md) const;
	virtual void free_hamiltonians_qj(RunParam * rp, ConfigParam * cp, MainData * md) const;
};