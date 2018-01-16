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
#include "qj_experiment.h"

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
	QJExperimentBehavior * qj_eb;

public:
	Processor(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd);

	void set_debug_behaviour(DebugBehavior* db);

	void set_output_behaviour(OuputBehavior* ob);

	void set_init_behaviour(InitBehavior* ib);

	void set_free_behaviour(FreeBehavior* fb);

	void set_qj_init_behaviour(QJInitBehavior* qj_ib);

	void set_qj_free_behaviour(QJFreeBehavior* fb);

	void set_qj_experiment_behaviour(QJExperimentBehavior* eb);

	void process();
};