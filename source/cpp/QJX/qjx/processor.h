#pragma once
#include "config.h"
#include "data.h"
#include "qj_data.h"
#include "debugger.h"
#include "outputter.h"
#include "initiator.h"
#include "destructor.h"
#include "qj_initiator.h"
#include "qj_destructor.h"

class Processor
{
private:
	RunParam * rp;
	ConfigParam * cp;
	MainData * md;
	QJData * qjd;

	DebugBehavior * db;
	OuputBehavior * ob;
	InitBehavior * ib;
	FreeBehavior * fb;
	QJInitBehavior * qj_ib;
	QJFreeBehavior * qj_fb;
	
public:
	Processor(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd)
	{
		this->rp = rp;
		this->cp = cp;
		this->md = md;
		this->qjd = qjd;
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

	void set_free_behaviour(FreeBehavior* fb)
	{
		this->fb = fb;
	}

	void set_qj_init_behaviour(QJInitBehavior* qj_ib)
	{
		this->qj_ib = qj_ib;
	}

	void set_qj_free_behaviour(QJFreeBehavior* fb)
	{
		this->qj_fb = qj_fb;
	}
};