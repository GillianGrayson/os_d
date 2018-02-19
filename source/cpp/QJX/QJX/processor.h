#pragma once
#include "config.h"
#include "data.h"
#include "data_qj.h"
#include "debugger.h"
#include "outputter.h"
#include "newdel.h"
#include "newdel_qj.h"
#include "experiment.h"

class Processor
{
private:
	RunParam * rp;
	ConfigParam * cp;
	MainData * md;
	QJData * qjd;

	DebugBehavior * db;
	OuputBehavior * ob;
	NewDelBehavior * ndb;
	QJNewDelBehavior * ndb_qj;
	ExperimentBehavior * eb;

public:
	Processor(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd);

	void set_debug_behaviour(DebugBehavior* db);

	void set_output_behaviour(OuputBehavior* ob);

	void set_init_behaviour(NewDelBehavior* ndb);

	void set_init_behaviour_qj(QJNewDelBehavior* ndb_qj);

	void set_experiment_behaviour(ExperimentBehavior* eb);

	void process();
};