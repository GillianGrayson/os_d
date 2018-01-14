#include "processor.h"

int main()
{
	RunParam * rp = new RunParam();
	ConfigParam * cp = new ConfigParam();
	MainData * md = new MainData();
	QJData * qjd = new QJData();

	init_params(rp, cp, "config.txt", "params.txt");

	DebugBehavior * db;
	OuputBehavior * ob;
	InitBehavior * ib;
	FreeBehavior * fb;
	QJInitBehavior * qj_ib;
	QJFreeBehavior * qj_fb;
	QJExperimentBehavior * qj_eb;

	if (rp->sys_id == DIMER_SYS_ID)
	{
		db = new DimerDebugBehaviour();
		ob = new DimerIOuputBehavior();
		ib = new DimerInitBehaviour();
		fb = new DimerFreeBehaviour();
	}
	else
	{
		stringstream msg;
		msg << "Error: Wrong sys_id: " << rp->sys_id << endl;
		Error(msg.str());
	}

	if (rp->task_id == LPN_TASK_ID)
	{
		qj_ib = new LpnInitBehaviour();
		qj_fb = new LpnFreeBehaviour();
		qj_eb = new LpnExperimentBehaviour();
	}
	else
	{
		stringstream msg;
		msg << "Error: Wrong task_id: " << rp->task_id << endl;
		Error(msg.str());
	}

	Processor * p = new Processor(rp, cp, md, qjd);

	p->set_debug_behaviour(db);
	p->set_output_behaviour(ob);
	p->set_init_behaviour(ib);
	p->set_free_behaviour(fb);

	p->set_qj_init_behaviour(qj_ib);
	p->set_qj_free_behaviour(qj_fb);
	p->set_qj_experiment_behaviour(qj_eb);

	p->process();

	delete db;
	delete ob;
	delete ib;
	delete fb;

	delete qj_ib;
	delete qj_fb;
	delete qj_eb;

	delete rp;
	delete cp;
	delete md;
	delete qjd;

	return 0;
}