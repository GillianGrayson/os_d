#pragma once
#include "config.h"
#include "data.h"

class CoreBehavior
{
public:
	virtual void init_splits(AllData * ad) const = 0;
	virtual void free_splits(AllData * ad) const = 0;

	virtual void init_splits_deep(AllData * ad) const = 0;
	virtual void free_splits_deep(AllData * ad) const = 0;
};

class DimerCoreBehaviour : public CoreBehavior
{
public:
	virtual void init_splits(AllData * ad) const;
	virtual void free_splits(AllData * ad) const;

	virtual void init_splits_deep(AllData * ad) const;
	virtual void free_splits_deep(AllData * ad) const;
};

class JCSCoreBehaviour : public CoreBehavior
{
public:
	virtual void init_splits(AllData * ad) const;
	virtual void free_splits(AllData * ad) const;

	virtual void init_splits_deep(AllData * ad) const;
	virtual void free_splits_deep(AllData * ad) const;
};

Split * init_split_structure_dimer(AllData * ad);
Split * init_split_structure_dimer_deep(AllData * ad);
Split * init_split_structure_jcs(AllData * ad);
Split * init_split_structure_jcs_deep(AllData * ad);