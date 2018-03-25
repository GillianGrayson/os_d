#pragma once
#include "config.h"
#include "data.h"
#include "debugger.h"
#include "outputter.h"
#include "newdel.h"
#include "newdel_exp.h"
#include "experiment.h"
#include "propagator.h"
#include "core.h"

class Processor
{
private:
	AllData * ad;

	DebugBehavior * db;
	OuputBehavior * ob;
	NewDelBehavior * ndb;
	ExpNewDelBehavior * ndb_exp;
	ExperimentBehavior * eb;
	PropagateBehavior * pb;
	CoreBehavior * cb;

public:
	Processor(AllData * ad);

	void set_debug_behaviour(DebugBehavior* db);

	void set_output_behaviour(OuputBehavior* ob);

	void set_newdel_behaviour(NewDelBehavior* ndb);

	void set_newdel_exp_behaviour(ExpNewDelBehavior* ndb_exp);

	void set_experiment_behaviour(ExperimentBehavior* eb);

	void set_propagate_behaviour(PropagateBehavior * pb);

	void set_core_behaviour(CoreBehavior * cb);

	void process();
};