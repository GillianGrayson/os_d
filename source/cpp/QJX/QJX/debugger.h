#pragma once
#include "config.h"
#include "data.h"

class DebugBehavior
{
public:
	virtual void save(AllData * ad) const = 0;
};

class DimerDebugBehaviour : public DebugBehavior
{
public:
	virtual void save(AllData * ad) const;
};

class JCSDebugBehaviour : public DebugBehavior
{
public:
	virtual void save(AllData * ad) const;
};

class PSDebugBehaviour : public DebugBehavior
{
public:
	virtual void save(AllData * ad) const;
};

class MBLDebugBehaviour : public DebugBehavior
{
public:
	virtual void save(AllData * ad) const;
};

class LndHamDebugBehaviour : public DebugBehavior
{
public:
	virtual void save(AllData* ad) const;
};

class IntegrableDebugBehaviour : public DebugBehavior
{
public:
	virtual void save(AllData* ad) const;
};

void save_hamiltonian_and_dissipation(AllData* ad, bool save_diss = true);
