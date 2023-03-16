#include "processor.h"

//int main()
//{
//	double start_time = omp_get_wtime();
//
//	RunParam * rp = new RunParam();
//	ConfigParam * cp = new ConfigParam();
//	init_params(rp, cp, "config.txt", "params.txt");
//
//	MainData* md = new MainData();
//	ExpData* ed = new ExpData();
//	AllData* ad = new AllData();
//
//	DebugBehavior * db = nullptr;
//	OutputBehavior * ob = nullptr;
//	NewDelBehavior * ndb = nullptr;
//	ExpNewDelBehavior * ndb_exp = nullptr;
//	ExperimentBehavior * eb = nullptr;
//	PropagateBehavior * pb = nullptr;
//	CoreBehavior * cb = nullptr;
//
//	if (rp->sys_id == DIMER_SYS_ID)
//	{
//		db = new DimerDebugBehaviour();
//		ob = new DimerOutputBehavior();
//		ndb = new DimerNewDelBehaviour();
//		cb = new DimerCoreBehaviour();
//	}
//	else if (rp->sys_id == DIMERSYNC_SYS_ID)
//	{
//		db = new DimerSyncDebugBehaviour();
//		ob = new DimerSyncOutputBehavior();
//		ndb = new DimerSyncNewDelBehaviour();
//		cb = new DimerSyncCoreBehaviour();
//	}
//	else if(rp->sys_id == JCS_SYS_ID)
//	{
//		db = new JCSDebugBehaviour();
//		ob = new JCSOutputBehavior();
//		ndb = new JCSNewDelBehaviour();
//		cb = new JCSCoreBehaviour();
//	}
//	else if (rp->sys_id == PS_SYS_ID)
//	{
//		db = new PSDebugBehaviour();
//		ob = new PSOutputBehavior();
//		ndb = new PSNewDelBehaviour();
//		cb = new PSCoreBehaviour();
//	}
//	else if (rp->sys_id == MBL_SYS_ID)
//	{
//		db = new MBLDebugBehaviour();
//		ob = new MBLOutputBehavior();
//		ndb = new MBLNewDelBehaviour();
//		cb = new MBLCoreBehaviour();
//	}
//	else if (rp->sys_id == LNDHAM_SYS_ID)
//	{
//		db = new LndHamDebugBehaviour();
//		ob = new LndHamOutputBehavior();
//		ndb = new LndHamNewDelBehaviour();
//		cb = new LndHamCoreBehaviour();
//	}
//	else if (rp->sys_id == INTEGRABLE_SYS_ID)
//	{
//		db = new IntegrableDebugBehaviour();
//		ob = new IntegrableOutputBehavior();
//		ndb = new IntegrableNewDelBehaviour();
//		cb = new IntegrableCoreBehaviour();
//	}
//	else if (rp->sys_id == FLOQ_2_SPINS_SYS_ID)
//	{
//		db = new Floq2SpinsDebugBehaviour();
//		ob = new Floq2SpinsOutputBehavior();
//		ndb = new Floq2SpinsNewDelBehaviour();
//		cb = new Floq2SpinsCoreBehaviour();
//	}
//	else
//	{
//		stringstream msg;
//		msg << "Error: Wrong sys_id: " << rp->sys_id << endl;
//		Error(msg.str());
//	}
//
//	if (rp->task_id == LPN_TASK_ID)
//	{
//		ndb_exp = new LpnNewDelBehaviour();
//		eb = new LpnExperimentBehaviour();
//	}
//	else if (rp->task_id == STD_TASK_ID)
//	{
//		ndb_exp = new StdNewDelBehaviour();
//		eb = new StdExperimentBehaviour();
//	}
//	else if (rp->task_id == LPN_MULT_TASK_ID)
//	{
//		ndb_exp = new LpnMultNewDelBehaviour();
//		eb = new LpnMultExperimentBehaviour();
//	}
//	else if (rp->task_id == LPN_MULT_DEEP_TASK_ID)
//	{
//		ndb_exp = new LpnMultDeepNewDelBehaviour();
//		eb = new LpnMultDeepExperimentBehaviour();
//	}
//	else if (rp->task_id == CD_TASK_ID)
//	{
//		ndb_exp = new CorrDimNewDelBehaviour();
//		eb = new CorrDimExperimentBehaviour();
//	}
//	else if (rp->task_id == SIGMA_TASK_ID)
//	{
//		ndb_exp = new SigmaNewDelBehaviour();
//		eb = new SigmaExperimentBehaviour();
//	}
//	else if (rp->task_id == STD_DEEP_TASK_ID)
//	{
//		ndb_exp = new StdDeepNewDelBehaviour();
//		eb = new StdDeepExperimentBehaviour();
//	}
//	else if (rp->task_id == LPN_DEEP_TASK_ID)
//	{
//		ndb_exp = new LpnDeepNewDelBehaviour();
//		eb = new LpnDeepExperimentBehaviour();
//	}
//	else if (rp->task_id == LPN_DEEP_TASK_ID_PER_1T)
//	{
//		ndb_exp = new LpnDeepNewDelBehaviour();
//		eb = new LpnDeepPer1TExperimentBehaviour();
//	}
//	else if (rp->task_id == LPN_ALL_TASK_ID)
//	{
//		ndb_exp = new LpnNewDelBehaviour();
//		eb = new LpnAllExperimentBehaviour();
//	}
//	else
//	{
//		stringstream msg;
//		msg << "Error: Wrong task_id: " << rp->task_id << endl;
//		Error(msg.str());
//	}
//
//	if (rp->prop_id == QJ_PROP_TYPE)
//	{
//		pb = new QJPropagateBehavior();
//	}
//	else if (rp->prop_id == RK_PROP_TYPE)
//	{
//		pb = new RKPropagateBehavior();
//	}
//	else
//	{
//		stringstream msg;
//		msg << "Error: Wrong prop_id: " << rp->prop_id << endl;
//		Error(msg.str());
//	}
//
//	ad->rp = rp;
//	ad->cp = cp;
//	ad->md = md;
//	ad->ed = ed;
//
//	Processor * p = new Processor(ad);
//
//	p->set_debug_behaviour(db);
//	p->set_output_behaviour(ob);
//	p->set_newdel_behaviour(ndb);
//	p->set_newdel_exp_behaviour(ndb_exp);
//	p->set_experiment_behaviour(eb);
//	p->set_propagate_behaviour(pb);
//	p->set_core_behaviour(cb);
//
//	p->process();
//
//	delete db;
//	delete ob;
//	delete ndb;
//	delete ndb_exp;
//	delete eb;
//	delete pb;
//	delete cb;
//
//	delete rp;
//	delete cp;
//	delete md;
//	delete ed;
//
//	double time = omp_get_wtime() - start_time;
//	cout << "Elapsed time: " << time << endl;
//
//	return 0;
//}


int main()
{
	double start_time = omp_get_wtime();

	RunParam* rp = new RunParam();
	ConfigParam* cp = new ConfigParam();
	init_params(rp, cp, "config.txt", "params.txt");

	double ampl_start = 1.0;
	double ampl_shift = 1.0;
	int ampl_num = 100;

	double omega_start = 0.1;
	double omega_shift = 0.1;
	int omega_num = 100;

	for (int ampl_id = 0; ampl_id < ampl_num; ampl_id++)
	{
		double ampl = ampl_start + ampl_shift * ampl_id;

		for (int omega_id = 0; omega_id < omega_num; omega_id++)
		{
			double omega = omega_start + omega_shift * omega_id;
			double T = 2.0 * PI / omega;

			std::cout << endl << "ampl: " << ampl << endl;
			std::cout << "omega: " << omega << endl;
			
			cp->rk_ns = max(1000, static_cast<int>(T / (0.005)));

			double step = T / cp->rk_ns;

			std::cout << "step: " << step << endl ;
			std::cout << "rk_ns: " << cp->rk_ns << endl << endl;

			auto it = cp->params.find("flq2sp_ampl");
			if (it != cp->params.end())
				it->second = ampl;

			it = cp->params.find("flq2sp_freq");
			if (it != cp->params.end())
				it->second = omega;

			MainData* md = new MainData();
			ExpData* ed = new ExpData();
			AllData* ad = new AllData();

			DebugBehavior* db = nullptr;
			OutputBehavior* ob = nullptr;
			NewDelBehavior* ndb = nullptr;
			ExpNewDelBehavior* ndb_exp = nullptr;
			ExperimentBehavior* eb = nullptr;
			PropagateBehavior* pb = nullptr;
			CoreBehavior* cb = nullptr;

			if (rp->sys_id == DIMER_SYS_ID)
			{
				db = new DimerDebugBehaviour();
				ob = new DimerOutputBehavior();
				ndb = new DimerNewDelBehaviour();
				cb = new DimerCoreBehaviour();
			}
			else if (rp->sys_id == DIMERSYNC_SYS_ID)
			{
				db = new DimerSyncDebugBehaviour();
				ob = new DimerSyncOutputBehavior();
				ndb = new DimerSyncNewDelBehaviour();
				cb = new DimerSyncCoreBehaviour();
			}
			else if (rp->sys_id == JCS_SYS_ID)
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
			else if (rp->sys_id == FLOQ_2_SPINS_SYS_ID)
			{
				db = new Floq2SpinsDebugBehaviour();
				ob = new Floq2SpinsOutputBehavior();
				ndb = new Floq2SpinsNewDelBehaviour();
				cb = new Floq2SpinsCoreBehaviour();
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

			Processor* p = new Processor(ad);

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

			delete md;
			delete ed;
		}
	}

	delete rp;
	delete cp;

	double time = omp_get_wtime() - start_time;
	cout << "Elapsed time: " << time << endl;

	return 0;
}