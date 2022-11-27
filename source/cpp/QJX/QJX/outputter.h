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

class DimerSyncOutputBehavior : public OutputBehavior
{
public:
	virtual void suffix_param(RunParam* rp, ConfigParam* cp, int precision) const;
};

class JCSOutputBehavior : public OutputBehavior
{
public:
	virtual void suffix_param(RunParam * rp, ConfigParam * cp, int precision) const;
};

class PSOutputBehavior : public OutputBehavior
{
public:
	virtual void suffix_param(RunParam * rp, ConfigParam * cp, int precision) const;
};

class MBLOutputBehavior : public OutputBehavior
{
public:
	virtual void suffix_param(RunParam * rp, ConfigParam * cp, int precision) const;
};

class LndHamOutputBehavior : public OutputBehavior
{
public:
	virtual void suffix_param(RunParam* rp, ConfigParam* cp, int precision) const;
};

class IntegrableOutputBehavior : public OutputBehavior
{
public:
	virtual void suffix_param(RunParam* rp, ConfigParam* cp, int precision) const;
};

string suffix_qj(RunParam * rp, ConfigParam * cp, int precision);

string suffix_setup(RunParam * rp);

string suffix_lpn(RunParam * rp, ConfigParam * cp);

string extension();