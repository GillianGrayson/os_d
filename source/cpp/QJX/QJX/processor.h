#pragma once
#include "config.h"
#include "data.h"
#include "debugger.h"
#include "outputter.h"
#include "newdel.h"
#include "newdel_exp.h"
#include "experiment.h"
#include "propagator.h"

class Processor
{
private:
	AllData * ad;

	DebugBehavior * db;
	OuputBehavior * ob;
	NewDelBehavior * ndb;
	ExpNewDelBehavior * ndb_qj;
	ExperimentBehavior * eb;
	PropagateBehavior * pb;

public:
	Processor(AllData * ad);

	void set_debug_behaviour(DebugBehavior* db);

	void set_output_behaviour(OuputBehavior* ob);

	void set_init_behaviour(NewDelBehavior* ndb);

	void set_init_behaviour_qj(ExpNewDelBehavior* ndb_qj);

	void set_experiment_behaviour(ExperimentBehavior* eb);

	void set_propagate_behaviour(PropagateBehavior * pb);

	void process();
};