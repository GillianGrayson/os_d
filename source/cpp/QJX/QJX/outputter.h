#pragma once
#include "config.h"
#include "data.h"

class OuputBehavior
{
public:
	virtual void suffix_param(RunParam * rp, ConfigParam * cp, int precision) const = 0;
};

class DimerIOuputBehavior : public OuputBehavior
{
public:
	virtual void suffix_param(RunParam * rp, ConfigParam * cp, int precision) const;
};

string suffix_qj(RunParam * rp, ConfigParam * cp, int precision);

string extension();