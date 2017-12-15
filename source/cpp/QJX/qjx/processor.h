#pragma once
#include "config.h"
#include "data.h"
#include "debugger.h"
#include "outputter.h"
#include "initiator.h"

class Processor
{
private:
	RunParam * rp;
	ConfigParam * cp;
	MainData *md;

	DebugBehavior * db;
	OuputBehavior * ob;
	InitBehavior * ib;

public:
	Processor(RunParam * rp, ConfigParam * cp, MainData *md)
	{
		this->rp = rp;
		this->cp = cp;
		this->md = md;
	}

	void set_debug_behaviour(DebugBehavior* db) 
	{
		this->db = db;
	}

	void set_output_behaviour(OuputBehavior* ob) 
	{
		this->ob = ob;
	}

	void set_init_behaviour(InitBehavior* ib) 
	{
		this->ib = ib;
	}

	void debug_save()
	{
		db->save(rp, cp, md);
	}

	void init_suffix()
	{
		ob->init_suffix(cp, 4);
	}
};