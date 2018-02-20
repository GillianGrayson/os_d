#pragma once
#include "config.h"
#include "data.h"
#include "debugger.h"
#include "outputter.h"
#include "newdel.h"
#include "newdel_qj.h"
#include "experiment.h"

class Processor
{
private:
	AllData * ad;

	DebugBehavior * db;
	OuputBehavior * ob;
	NewDelBehavior * ndb;
	QJNewDelBehavior * ndb_qj;
	ExperimentBehavior * eb;

public:
	Processor(AllData * ad);

	void set_debug_behaviour(DebugBehavior* db);

	void set_output_behaviour(OuputBehavior* ob);

	void set_init_behaviour(NewDelBehavior* ndb);

	void set_init_behaviour_qj(QJNewDelBehavior* ndb_qj);

	void set_experiment_behaviour(ExperimentBehavior* eb);

	void process();
};