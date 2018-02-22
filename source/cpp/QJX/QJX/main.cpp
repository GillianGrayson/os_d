#include "processor.h"

int main()
{
	double start_time = omp_get_wtime();

	RunParam * rp = new RunParam();
	ConfigParam * cp = new ConfigParam();
	MainData * md = new MainData();
	ExpData * ed = new ExpData();

	AllData * ad = new AllData();

	init_params(rp, cp, "config.txt", "params.txt");

	DebugBehavior * db;
	OuputBehavior * ob;
	NewDelBehavior * ndb;
	ExpNewDelBehavior * ndb_exp;
	ExperimentBehavior * eb;
	PropagateBehavior * pb;

	if (rp->sys_id == DIMER_SYS_ID)
	{
		db = new DimerDebugBehaviour();
		ob = new DimerIOuputBehavior();
		ndb = new DimerNewDelBehaviour();
	}
	else
	{
		stringstream msg;
		msg << "Error: Wrong sys_id: " << rp->sys_id << endl;
		Error(msg.str());
	}

	if (rp->task_id == LPN_TASK_ID)
	{
		ndb_exp = new LpnNewDelBehaviour();
		eb = new LpnExperimentBehaviour();
	}
	else if (rp->task_id == STD_TASK_ID)
	{
		ndb_exp = new StdNewDelBehaviour();
		eb = new StdExperimentBehaviour();
	}
	else if (rp->task_id == CD_TASK_ID)
	{
		ndb_exp = new CorrDimNewDelBehaviour();
		eb = new CorrDimExperimentBehaviour();
	}
	else if (rp->task_id == SIGMA_TASK_ID)
	{
		ndb_exp = new SigmaNewDelBehaviour();
		eb = new SigmaExperimentBehaviour();
	}
	else
	{
		stringstream msg;
		msg << "Error: Wrong task_id: " << rp->task_id << endl;
		Error(msg.str());
	}

	if (rp->prop_id == QJ_PROP_TYPE)
	{
		pb = new QJPropagateBehavior();
	}
	else
	{
		stringstream msg;
		msg << "Error: Wrong prop_id: " << rp->prop_id << endl;
		Error(msg.str());
	}

	ad->rp = rp;
	ad->cp = cp;
	ad->md = md;
	ad->ed = ed;

	Processor * p = new Processor(ad);

	p->set_debug_behaviour(db);
	p->set_output_behaviour(ob);
	p->set_newdel_behaviour(ndb);

	p->set_newdel_exp_behaviour(ndb_exp);
	p->set_experiment_behaviour(eb);
	p->set_propagate_behaviour(pb);

	p->process();

	delete db;
	delete ob;
	delete ndb;

	delete ndb_exp;
	delete eb;

	delete rp;
	delete cp;
	delete md;
	delete ed;

	double time = omp_get_wtime() - start_time;
	cout << "Elapsed time: " << time << endl;

	return 0;
}