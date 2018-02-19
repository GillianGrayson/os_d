#pragma once
#include "config.h"
#include "data.h"

class NewDelBehavior
{
public:
	virtual void init_sizes(RunParam * rp, ConfigParam * cp, MainData * md) const = 0;
	virtual void init_hamiltonians(RunParam * rp, ConfigParam * cp, MainData * md) const = 0;
	virtual void init_dissipators(RunParam * rp, ConfigParam * cp, MainData * md) const = 0;
	virtual void init_hamiltonians_qj(RunParam * rp, ConfigParam * cp, MainData * md) const = 0;

	virtual void free_hamiltonians(RunParam * rp, ConfigParam * cp, MainData * md) const = 0;
	virtual void free_dissipators(RunParam * rp, ConfigParam * cp, MainData * md) const = 0;
	virtual void free_hamiltonians_qj(RunParam * rp, ConfigParam * cp, MainData * md) const = 0;
};

class DimerNewDelBehaviour : public NewDelBehavior
{
public:
	virtual void init_sizes(RunParam * rp, ConfigParam * cp, MainData * md) const;
	virtual void init_hamiltonians(RunParam * rp, ConfigParam * cp, MainData * md) const;
	virtual void init_dissipators(RunParam * rp, ConfigParam * cp, MainData * md) const;
	virtual void init_hamiltonians_qj(RunParam * rp, ConfigParam * cp, MainData * md) const;

	virtual void free_hamiltonians(RunParam * rp, ConfigParam * cp, MainData * md) const;
	virtual void free_dissipators(RunParam * rp, ConfigParam * cp, MainData * md) const;
	virtual void free_hamiltonians_qj(RunParam * rp, ConfigParam * cp, MainData * md) const;
};
