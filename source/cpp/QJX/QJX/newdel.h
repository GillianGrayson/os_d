#pragma once
#include "config.h"
#include "data.h"
#include "propagator.h"

class NewDelBehavior
{
public:
	virtual void init_sizes(AllData * ad) const = 0;
	virtual void init_hamiltonians(AllData * ad) const = 0;
	virtual void init_dissipators(AllData * ad) const = 0;
	virtual void init_hamiltonians_qj(AllData * ad) const = 0;

	virtual void free_hamiltonians(AllData * ad) const = 0;
	virtual void free_dissipators(AllData * ad) const = 0;
	virtual void free_hamiltonians_qj(AllData * ad) const = 0;
};

class DimerNewDelBehaviour : public NewDelBehavior
{
public:
	virtual void init_sizes(AllData * ad) const;
	virtual void init_hamiltonians(AllData * ad) const;
	virtual void init_dissipators(AllData * ad) const;
	virtual void init_hamiltonians_qj(AllData * ad) const;

	virtual void free_hamiltonians(AllData * ad) const;
	virtual void free_dissipators(AllData * ad) const;
	virtual void free_hamiltonians_qj(AllData * ad) const;
};

class JCSNewDelBehaviour : public NewDelBehavior
{
public:
	virtual void init_sizes(AllData * ad) const;
	virtual void init_hamiltonians(AllData * ad) const;
	virtual void init_dissipators(AllData * ad) const;
	virtual void init_hamiltonians_qj(AllData * ad) const;

	virtual void free_hamiltonians(AllData * ad) const;
	virtual void free_dissipators(AllData * ad) const;
	virtual void free_hamiltonians_qj(AllData * ad) const;
};


class PSNewDelBehaviour : public NewDelBehavior
{
public:
	virtual void init_sizes(AllData * ad) const;
	virtual void init_hamiltonians(AllData * ad) const;
	virtual void init_dissipators(AllData * ad) const;
	virtual void init_hamiltonians_qj(AllData * ad) const;

	virtual void free_hamiltonians(AllData * ad) const;
	virtual void free_dissipators(AllData * ad) const;
	virtual void free_hamiltonians_qj(AllData * ad) const;
};

class MBLNewDelBehaviour : public NewDelBehavior
{
public:
	virtual void init_sizes(AllData * ad) const;
	virtual void init_hamiltonians(AllData * ad) const;
	virtual void init_dissipators(AllData * ad) const;
	virtual void init_hamiltonians_qj(AllData * ad) const;

	virtual void free_hamiltonians(AllData * ad) const;
	virtual void free_dissipators(AllData * ad) const;
	virtual void free_hamiltonians_qj(AllData * ad) const;
};