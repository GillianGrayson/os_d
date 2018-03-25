#pragma once
#include "config.h"
#include "data.h"

class OutputBehavior
{
public:
	virtual void suffix_param(RunParam * rp, ConfigParam * cp, int precision) const = 0;
};

class DimerOutputBehavior : public OutputBehavior
{
public:
	virtual void suffix_param(RunParam * rp, ConfigParam * cp, int precision) const;
};

class JCSOutputBehavior : public OutputBehavior
{
public:
	virtual void suffix_param(RunParam * rp, ConfigParam * cp, int precision) const;
};

string suffix_qj(RunParam * rp, ConfigParam * cp, int precision);

string extension();