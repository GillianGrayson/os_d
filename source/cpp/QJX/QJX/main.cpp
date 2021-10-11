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

	DebugBehavior * db = nullptr;
	OutputBehavior * ob = nullptr;
	NewDelBehavior * ndb = nullptr;
	ExpNewDelBehavior * ndb_exp = nullptr;
	ExperimentBehavior * eb = nullptr;
	PropagateBehavior * pb = nullptr;
	CoreBehavior * cb = nullptr;

	if (rp->sys_id == DIMER_SYS_ID)
	{
		db = new DimerDebugBehaviour();
		ob = new DimerOutputBehavior();
		ndb = new DimerNewDelBehaviour();
		cb = new DimerCoreBehaviour();
	}
	else if(rp->sys_id == JCS_SYS_ID)
	{
		db = new JCSDebugBehaviour();
		ob = new JCSOutputBehavior();
		ndb = new JCSNewDelBehaviour();
		cb = new JCSCoreBehaviour();
	}
	else if (rp->sys_id == PS_SYS_ID)
	{
		db = new PSDebugBehaviour();
		ob = new PSOutputBehavior();
		ndb = new PSNewDelBehaviour();
		cb = new PSCoreBehaviour();
	}
	else if (rp->sys_id == MBL_SYS_ID)
	{
		db = new MBLDebugBehaviour();
		ob = new MBLOutputBehavior();
		ndb = new MBLNewDelBehaviour();
		cb = new MBLCoreBehaviour();
	}
	else if (rp->sys_id == LNDHAM_SYS_ID)
	{
		db = new LndHamDebugBehaviour();
		ob = new LndHamOutputBehavior();
		ndb = new LndHamNewDelBehaviour();
		cb = new LndHamCoreBehaviour();
	}
	else if (rp->sys_id == INTEGRABLE_SYS_ID)
	{
		db = new IntegrableDebugBehaviour();
		ob = new IntegrableOutputBehavior();
		ndb = new IntegrableNewDelBehaviour();
		cb = new IntegrableCoreBehaviour();
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
	else if (rp->task_id == LPN_MULT_TASK_ID)
	{
		ndb_exp = new LpnMultNewDelBehaviour();
		eb = new LpnMultExperimentBehaviour();
	}
	else if (rp->task_id == LPN_MULT_DEEP_TASK_ID)
	{
		ndb_exp = new LpnMultDeepNewDelBehaviour();
		eb = new LpnMultDeepExperimentBehaviour();
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
	else if (rp->task_id == STD_DEEP_TASK_ID)
	{
		ndb_exp = new StdDeepNewDelBehaviour();
		eb = new StdDeepExperimentBehaviour();
	}
	else if (rp->task_id == LPN_DEEP_TASK_ID)
	{
		ndb_exp = new LpnDeepNewDelBehaviour();
		eb = new LpnDeepExperimentBehaviour();
	}
	else if (rp->task_id == LPN_DEEP_TASK_ID_PER_1T)
	{
		ndb_exp = new LpnDeepNewDelBehaviour();
		eb = new LpnDeepPer1TExperimentBehaviour();
	}
	else if (rp->task_id == LPN_ALL_TASK_ID)
	{
		ndb_exp = new LpnNewDelBehaviour();
		eb = new LpnAllExperimentBehaviour();
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
	else if (rp->prop_id == RK_PROP_TYPE)
	{
		pb = new RKPropagateBehavior();
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
	p->set_core_behaviour(cb);

	p->process();

	delete db;
	delete ob;
	delete ndb;
	delete ndb_exp;
	delete eb;
	delete pb;
	delete cb;

	delete rp;
	delete cp;
	delete md;
	delete ed;

	double time = omp_get_wtime() - start_time;
	cout << "Elapsed time: " << time << endl;

	return 0;
}