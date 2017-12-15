#pragma once
#include "config.h"
#include "data.h"

class OuputBehavior
{
public:
	virtual void init_suffix(ConfigParam * cp, int precision) const = 0;
};

class DimerIOuputBehavior : public OuputBehavior
{
public:
	virtual void init_suffix(ConfigParam * cp, int precision) const;
};