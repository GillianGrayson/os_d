#include "experiment.h"

void LpnExperimentBehaviour::trans_process(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const
{
	ConfigParam * cp = ad->cp;

	int num_trajectories = cp->num_trajectories;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);
	int dump_evo_avg = int(cp->params.find("dump_evo_avg")->second);

#pragma omp parallel for
	for (int tr_id = 0; tr_id < 1; tr_id++)
	{
		int thread_id = omp_get_thread_num();
		trans_process_single(ad, pb, cb, tr_id, thread_id);
	}

	cb->calc_chars_std_start(ad, 0);
	cb->calc_chars_lpn_start(ad, 0, 0);

#pragma omp parallel for
	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		if (tr_id > 0)
		{
			copy_trajectory_lpn(ad, tr_id, 0);
			var_trajectory_lpn(ad, cb, tr_id, 0);
		}

		resresh_times(ad, tr_id);

		cb->calc_chars_std_start(ad, tr_id);
		cb->calc_chars_lpn_start(ad, tr_id, 0);

		cb->evo_chars_std(ad, tr_id, 0);
		cb->evo_chars_lpn(ad, tr_id, 0);

		if (dump_evo_sep == 1)
		{
			dump_adr_single(ad, tr_id, false);
		}
	}

	if (dump_evo_avg == 1)
	{
		dump_adr_avg(ad, false);
	}
}

void LpnExperimentBehaviour::obser_process(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const
{	
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	ExpData * ed = ad->ed;
	MainData * md = ad->md;

	ed->is_obs = 1;

	int num_trajectories = cp->num_trajectories;

	int dump_num_total = ed->dump_num_total;
	int * dump_periods = ed->dump_periods;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);
	int dump_evo_avg = int(cp->params.find("dump_evo_avg")->second);

	int begin_period_id = 0;
	int end_period_id = 0;
	for (int dump_id = 1; dump_id < dump_num_total; dump_id++)
	{
		if (rp->is_pp == 1)
		{
			cout << "dump_id: " << dump_id << endl;
		}

		begin_period_id = dump_periods[dump_id - 1];
		end_period_id = dump_periods[dump_id];

		for (int period_id = begin_period_id; period_id < end_period_id; period_id++)
		{
#pragma omp parallel for
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				int thread_id = omp_get_thread_num();
				pb->one_period(ad, cb, tr_id, thread_id, period_id);
				cb->calc_chars_std(ad, tr_id);
			}

			cb->calc_chars_lpn(ad, 0, 0);

			ed->curr_time = (period_id + 1) * cb->calc_T(ad);

#pragma omp parallel for
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				if (tr_id > 0)
				{
					lambda_lpn(ad, cb, tr_id, 0);
				}
			}
		}

#pragma omp parallel for
		for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
		{
			cb->evo_chars_std(ad, tr_id, dump_id);
			cb->evo_chars_lpn(ad, tr_id, dump_id);

			if (dump_evo_sep == 1)
			{
				dump_adr_single(ad, tr_id, true);
			}
		}

		if (dump_evo_avg == 1)
		{
			dump_adr_avg(ad, true);
		}
	}

	cb->dump_std(ad);
	cb->dump_lpn(ad);

	if (dump_evo_sep == 1)
	{
		cb->dump_std_evo(ad);
		cb->dump_lpn_evo(ad);
	}

	if (dump_evo_avg == 0)
	{
		dump_adr_avg(ad, true);
	}
}


void LpnMultExperimentBehaviour::trans_process(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const
{
	ConfigParam * cp = ad->cp;

	int num_trajectories = cp->num_trajectories;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);
	int dump_evo_avg = int(cp->params.find("dump_evo_avg")->second);

#pragma omp parallel for
	for (int tr_id = 0; tr_id < num_trajectories / 2; tr_id++)
	{
		int thread_id = omp_get_thread_num();
		trans_process_single(ad, pb, cb, tr_id, thread_id);
		cb->calc_chars_std_start(ad, tr_id);
		cb->calc_chars_lpn_start(ad, tr_id, tr_id);
	}

#pragma omp parallel for
	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		if (tr_id >= num_trajectories / 2)
		{
			copy_trajectory_lpn(ad, tr_id, tr_id - num_trajectories / 2);
			var_trajectory_lpn(ad, cb, tr_id, tr_id - num_trajectories / 2);
		}

		resresh_times(ad, tr_id);

		cb->calc_chars_std_start(ad, tr_id);
		if (tr_id >= num_trajectories / 2)
		{
			cb->calc_chars_lpn_start(ad, tr_id, tr_id - num_trajectories / 2);
		}
		else
		{
			cb->calc_chars_lpn_start(ad, tr_id, tr_id);
		}

		cb->evo_chars_std(ad, tr_id, 0);
		cb->evo_chars_lpn(ad, tr_id, 0);

		if (dump_evo_sep == 1)
		{
			dump_adr_single(ad, tr_id, false);
		}
	}

	if (dump_evo_avg == 1)
	{
		dump_adr_avg(ad, false);
	}
}


void LpnMultExperimentBehaviour::obser_process(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	ExpData * ed = ad->ed;
	MainData * md = ad->md;

	ed->is_obs = 1;

	int num_trajectories = cp->num_trajectories;

	int dump_num_total = ed->dump_num_total;
	int * dump_periods = ed->dump_periods;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);
	int dump_evo_avg = int(cp->params.find("dump_evo_avg")->second);

	int lambda_per_periods = int(cp->params.find("lambda_per_periods")->second);

	int begin_period_id = 0;
	int end_period_id = 0;
	for (int dump_id = 1; dump_id < dump_num_total; dump_id++)
	{
		if (rp->is_pp == 1)
		{
			cout << "dump_id: " << dump_id << endl;
		}

		begin_period_id = dump_periods[dump_id - 1];
		end_period_id = dump_periods[dump_id];

		for (int period_id = begin_period_id; period_id < end_period_id; period_id++)
		{
#pragma omp parallel for
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				int thread_id = omp_get_thread_num();
				pb->one_period(ad, cb, tr_id, thread_id, period_id);
				cb->calc_chars_std(ad, tr_id);
			}

			for (int tr_id = 0; tr_id < num_trajectories / 2; tr_id++)
			{
				cb->calc_chars_lpn(ad, tr_id, tr_id);
			}

			ed->curr_time = (period_id + 1) * cb->calc_T(ad);

#pragma omp parallel for
			for (int tr_id = num_trajectories / 2; tr_id < num_trajectories; tr_id++)
			{
				if (lambda_per_periods == 0)
				{
					lambda_lpn(ad, cb, tr_id, tr_id - num_trajectories / 2);
				}
				else
				{
					lambda_lpn_now(ad, cb, tr_id, tr_id - num_trajectories / 2);
				}
			}
		}

#pragma omp parallel for
		for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
		{
			cb->evo_chars_std(ad, tr_id, dump_id);
			cb->evo_chars_lpn(ad, tr_id, dump_id);

			if (dump_evo_sep == 1)
			{
				dump_adr_single(ad, tr_id, true);
			}
		}

		if (dump_evo_avg == 1)
		{
			dump_adr_avg_mult(ad, true, 0, num_trajectories / 2);
		}
	}

	cb->dump_std(ad);
	cb->dump_lpn(ad);

	if (dump_evo_sep == 1)
	{
		cb->dump_std_evo(ad);
		cb->dump_lpn_evo(ad);
	}

	if (dump_evo_avg == 0)
	{
		dump_adr_avg_mult(ad, true, 0, num_trajectories / 2);
	}
}


void LpnMultDeepExperimentBehaviour::trans_process(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const
{
	ConfigParam * cp = ad->cp;

	int num_trajectories = cp->num_trajectories;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);
	int dump_evo_avg = int(cp->params.find("dump_evo_avg")->second);

#pragma omp parallel for
	for (int tr_id = 0; tr_id < num_trajectories / 2; tr_id++)
	{
		int thread_id = omp_get_thread_num();
		trans_process_single_deep(ad, pb, cb, tr_id, thread_id);
		cb->calc_chars_std_start(ad, tr_id);
		cb->calc_chars_lpn_start(ad, tr_id, tr_id);
	}

#pragma omp parallel for
	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		if (tr_id >= num_trajectories / 2)
		{
			copy_trajectory_lpn(ad, tr_id, tr_id - num_trajectories / 2);
			var_trajectory_lpn(ad, cb, tr_id, tr_id - num_trajectories / 2);
		}

		resresh_times(ad, tr_id);

		cb->calc_chars_std_start(ad, tr_id);
		if (tr_id >= num_trajectories / 2)
		{
			cb->calc_chars_lpn_start(ad, tr_id, tr_id - num_trajectories / 2);
		}
		else
		{
			cb->calc_chars_lpn_start(ad, tr_id, tr_id);
		}

		cb->evo_chars_std(ad, tr_id, 0);
		cb->evo_chars_lpn(ad, tr_id, 0);

		if (dump_evo_sep == 1)
		{
			dump_adr_single(ad, tr_id, false);
		}
	}

	if (dump_evo_avg == 1)
	{
		dump_adr_avg(ad, false);
	}
}


void LpnMultDeepExperimentBehaviour::obser_process(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	ExpData * ed = ad->ed;
	MainData * md = ad->md;

	ed->is_obs = 1;

	int num_branches = md->num_ham_qj;
	int num_sub_steps = num_branches * int(cp->params.find("deep_num_steps")->second);

	int num_trajectories = cp->num_trajectories;

	int dump_num_total = cp->num_obs_periods + 1;

	int * dump_periods = ed->dump_periods;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);
	int dump_evo_avg = int(cp->params.find("dump_evo_avg")->second);

	int begin_period_id = 0;
	int end_period_id = 0;
	for (int dump_id = 1; dump_id < dump_num_total; dump_id++)
	{
		if (rp->is_pp == 1)
		{
			cout << "dump_id: " << dump_id << endl;
		}

		begin_period_id = dump_id - 1;
		end_period_id = dump_id;

		for (int period_id = begin_period_id; period_id < end_period_id; period_id++)
		{
			pb->one_period_obs_deep_mult_lpn(ad, cb, period_id);
		}

		if (dump_evo_avg == 1)
		{
			dump_adr_avg_mult(ad, true, 0, num_trajectories / 2);
		}
	}

	cb->dump_std(ad);
	cb->dump_lpn(ad);

	if (dump_evo_sep == 1)
	{
		cb->dump_std_evo(ad);
		cb->dump_lpn_evo(ad);
	}

	if (dump_evo_avg == 0)
	{
		dump_adr_avg_mult(ad, true, 0, num_trajectories / 2);
	}
}





void StdExperimentBehaviour::trans_process(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const
{
	ConfigParam * cp = ad->cp;

	int num_trajectories = cp->num_trajectories;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);
	int dump_evo_avg = int(cp->params.find("dump_evo_avg")->second);

#pragma omp parallel for
	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		int thread_id = omp_get_thread_num();
		trans_process_single(ad, pb, cb, tr_id, thread_id);
	}

#pragma omp parallel for
	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		resresh_times(ad, tr_id);

		cb->calc_chars_std_start(ad, tr_id);

		cb->evo_chars_std(ad, tr_id, 0);

		if (dump_evo_sep == 1)
		{
			dump_adr_single(ad, tr_id, false);
		}
	}

	if (dump_evo_avg == 1)
	{
		dump_adr_avg(ad, false);
	}
}

void StdExperimentBehaviour::obser_process(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	ExpData * ed = ad->ed;

	ed->is_obs = 1;

	int num_trajectories = cp->num_trajectories;

	int dump_num_total = ed->dump_num_total;
	int * dump_periods = ed->dump_periods;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);
	int dump_evo_avg = int(cp->params.find("dump_evo_avg")->second);

	int jumps_counts_lim = int(cp->params.find("jumps_counts")->second);

	dump_phi_evo(ad, false);

	int begin_period_id = 0;
	int end_period_id = 0;
	bool all_traj_finished = false;
	for (int dump_id = 1; dump_id < dump_num_total; dump_id++)
	{
		if (rp->is_pp == 1)
		{
			cout << "dump_id: " << dump_id << endl;
		}

		begin_period_id = dump_periods[dump_id - 1];
		end_period_id = dump_periods[dump_id];

		for (int period_id = begin_period_id; period_id < end_period_id; period_id++)
		{
			all_traj_finished = true;
#pragma omp parallel for
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				int thread_id = omp_get_thread_num();
				pb->one_period(ad, cb, tr_id, thread_id, period_id);
				
				if (ed->jumps_counts[tr_id] < jumps_counts_lim)
				{
					all_traj_finished = false;
				}
			}

			if (all_traj_finished == true && jumps_counts_lim > 0)
			{
				break;
			}
		}

#pragma omp parallel for
		for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
		{
			cb->calc_chars_std(ad, tr_id);

			cb->evo_chars_std(ad, tr_id, dump_id);

			if (dump_evo_sep == 1)
			{
				dump_adr_single(ad, tr_id, true);
			}
		}

		dump_phi_evo(ad, true);

		if (dump_evo_avg == 1)
		{
			dump_adr_avg(ad, true);
		}

		if (all_traj_finished == true && jumps_counts_lim > 0)
		{
			break;
		}
	}

	cb->dump_std(ad);

	dump_phi(ad);
	dump_phi_evo(ad, true);

	if (dump_evo_sep == 1)
	{
		cb->dump_std_evo(ad);
	}

	if (dump_evo_avg == 0)
	{
		dump_adr_avg(ad, false);
	}
}

void CorrDimExperimentBehaviour::trans_process(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const
{
	ConfigParam * cp = ad->cp;

	int num_trajectories = cp->num_trajectories;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);
	int dump_evo_avg = int(cp->params.find("dump_evo_avg")->second);

#pragma omp parallel for
	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		int thread_id = omp_get_thread_num();
		trans_process_single_deep(ad, pb, cb, tr_id, thread_id);
	}

#pragma omp parallel for
	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		resresh_times(ad, tr_id);

		cb->calc_chars_std_start(ad, tr_id);

		cb->evo_chars_std(ad, tr_id, 0);

		if (dump_evo_sep == 1)
		{
			dump_adr_single(ad, tr_id, false);
		}
	}

	if (dump_evo_avg == 1)
	{
		dump_adr_avg(ad, false);
	}
}

void CorrDimExperimentBehaviour::obser_process(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	ExpData * ed = ad->ed;

	ed->is_obs = 1;

	int num_trajectories = cp->num_trajectories;

	int dump_num_total = ed->dump_num_total;

	int * dump_periods = ed->dump_periods;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);
	int dump_evo_avg = int(cp->params.find("dump_evo_avg")->second);

	int begin_period_id = 0;
	int end_period_id = 0;
	for (int dump_id = 1; dump_id < dump_num_total; dump_id++)
	{
		if (rp->is_pp == 1)
		{
			cout << "dump_id: " << dump_id << endl;
		}

		begin_period_id = dump_periods[dump_id - 1];
		end_period_id = dump_periods[dump_id];

		for (int period_id = begin_period_id; period_id < end_period_id; period_id++)
		{
#pragma omp parallel for
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				int thread_id = omp_get_thread_num();
				pb->one_period_obs_deep_cd(ad, cb, tr_id, thread_id, period_id);
			}
		}

#pragma omp parallel for
		for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
		{
			cb->evo_chars_std(ad, tr_id, dump_id);

			if (dump_evo_sep == 1)
			{
				dump_adr_single(ad, tr_id, true);
			}
		}

		if (dump_evo_avg == 1)
		{
			dump_adr_avg(ad, true);
		}
	}

#pragma omp parallel for
	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		cb->calc_ci(ad, tr_id);
	}

	cb->dump_std(ad);

	dump_cd(ad);

	if (dump_evo_sep == 1)
	{
		cb->dump_std_evo(ad);
	}

	if (dump_evo_avg == 0)
	{
		dump_adr_avg(ad, false);
	}
}

void SigmaExperimentBehaviour::trans_process(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const
{
	ConfigParam * cp = ad->cp;

	int num_trajectories = cp->num_trajectories;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);
	int dump_evo_avg = int(cp->params.find("dump_evo_avg")->second);

#pragma omp parallel for
	for (int tr_id = 0; tr_id < 1; tr_id++)
	{
		int thread_id = omp_get_thread_num();
		trans_process_single_deep(ad, pb, cb, tr_id, thread_id);
	}

#pragma omp parallel for
	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		if (tr_id > 0)
		{
			copy_trajectory_data(ad, tr_id, 0);
		}

		resresh_times(ad, tr_id);

		cb->calc_chars_std_start(ad, tr_id);

		cb->evo_chars_std(ad, tr_id, 0);

		if (dump_evo_sep == 1)
		{
			dump_adr_single(ad, tr_id, false);
		}
	}

	if (dump_evo_avg == 1)
	{
		dump_adr_avg(ad, false);
	}
}

void SigmaExperimentBehaviour::obser_process(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	ExpData * ed = ad->ed;

	ed->is_obs = 1;

	int num_trajectories = cp->num_trajectories;

	int dump_num_total = ed->dump_num_total;

	int * dump_periods = ed->dump_periods;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);
	int dump_evo_avg = int(cp->params.find("dump_evo_avg")->second);

	int begin_period_id = 0;
	int end_period_id = 0;
	for (int dump_id = 1; dump_id < dump_num_total; dump_id++)
	{
		if (rp->is_pp == 1)
		{
			cout << "dump_id: " << dump_id << endl;
		}

		begin_period_id = dump_periods[dump_id - 1];
		end_period_id = dump_periods[dump_id];

		for (int period_id = begin_period_id; period_id < end_period_id; period_id++)
		{
#pragma omp parallel for
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				int thread_id = omp_get_thread_num();
				pb->one_period_obs_deep_sigma(ad, cb, tr_id, thread_id, period_id);
			}
		}

#pragma omp parallel for
		for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
		{
			cb->evo_chars_std(ad, tr_id, dump_id);

			if (dump_evo_sep == 1)
			{
				dump_adr_single(ad, tr_id, true);
			}
		}

		if (dump_evo_avg == 1)
		{
			dump_adr_avg(ad, true);
		}
	}

	cb->dump_std(ad);

	if (dump_evo_sep == 1)
	{
		cb->dump_std_evo(ad);
	}

	if (dump_evo_avg == 0)
	{
		dump_adr_avg(ad, false);
	}
}

void StdDeepExperimentBehaviour::trans_process(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const
{
	ConfigParam * cp = ad->cp;

	int num_trajectories = cp->num_trajectories;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);
	int dump_evo_avg = int(cp->params.find("dump_evo_avg")->second);

#pragma omp parallel for
	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		int thread_id = omp_get_thread_num();
		trans_process_single_deep(ad, pb, cb, tr_id, thread_id);
	}

#pragma omp parallel for
	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		resresh_times(ad, tr_id);

		cb->calc_chars_std_start(ad, tr_id);

		cb->evo_chars_std(ad, tr_id, 0);

		if (dump_evo_sep == 1)
		{
			dump_adr_single(ad, tr_id, false);
		}
	}

	if (dump_evo_avg == 1)
	{
		dump_adr_avg(ad, false);
	}
}

void StdDeepExperimentBehaviour::obser_process(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	ExpData * ed = ad->ed;

	ed->is_obs = 1;

	int num_trajectories = cp->num_trajectories;

	int dump_num_total = cp->num_obs_periods + 1;

	int * dump_periods = ed->dump_periods;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);
	int dump_evo_avg = int(cp->params.find("dump_evo_avg")->second);

	int begin_period_id = 0;
	int end_period_id = 0;
	for (int dump_id = 1; dump_id < dump_num_total; dump_id++)
	{
		if (rp->is_pp == 1)
		{
			cout << "dump_id: " << dump_id << endl;
		}

		begin_period_id = dump_id - 1;
		end_period_id = dump_id;

		for (int period_id = begin_period_id; period_id < end_period_id; period_id++)
		{
#pragma omp parallel for
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				int thread_id = omp_get_thread_num();
				pb->one_period_obs_deep(ad, cb, tr_id, thread_id, period_id);
			}
		}

		if (dump_evo_avg == 1)
		{
			dump_adr_avg(ad, true);
		}
	}

	cb->dump_std(ad);

	if (dump_evo_sep == 1)
	{
		cb->dump_std_evo(ad);
	}

	if (dump_evo_avg == 0)
	{
		dump_adr_avg(ad, false);
	}
}

void LpnDeepExperimentBehaviour::trans_process(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const
{
	ConfigParam * cp = ad->cp;

	int num_trajectories = cp->num_trajectories;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);
	int dump_evo_avg = int(cp->params.find("dump_evo_avg")->second);

#pragma omp parallel for
	for (int tr_id = 0; tr_id < 1; tr_id++)
	{
		int thread_id = omp_get_thread_num();
		trans_process_single_deep(ad, pb, cb, tr_id, thread_id);
	}

	cb->calc_chars_std_start(ad, 0);
	cb->calc_chars_lpn_start(ad, 0, 0);

#pragma omp parallel for
	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		if (tr_id > 0)
		{
			copy_trajectory_lpn(ad, tr_id, 0);
			var_trajectory_lpn(ad, cb, tr_id, 0);
		}

		resresh_times(ad, tr_id);

		cb->calc_chars_std_start(ad, tr_id);
		cb->calc_chars_lpn_start(ad, tr_id, 0);

		cb->evo_chars_std(ad, tr_id, 0);
		cb->evo_chars_lpn(ad, tr_id, 0);

		if (dump_evo_sep == 1)
		{
			dump_adr_single(ad, tr_id, false);
		}
	}

	if (dump_evo_avg == 1)
	{
		dump_adr_avg(ad, false);
	}
}

void LpnDeepExperimentBehaviour::obser_process(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	ExpData * ed = ad->ed;
	MainData * md = ad->md;

	ed->is_obs = 1;

	int num_branches = md->num_ham_qj;
	int num_sub_steps = num_branches * int(cp->params.find("deep_num_steps")->second);

	int num_trajectories = cp->num_trajectories;

	int dump_num_total = cp->num_obs_periods + 1;

	int * dump_periods = ed->dump_periods;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);
	int dump_evo_avg = int(cp->params.find("dump_evo_avg")->second);

	int begin_period_id = 0;
	int end_period_id = 0;
	for (int dump_id = 1; dump_id < dump_num_total; dump_id++)
	{
		if (rp->is_pp == 1)
		{
			cout << "dump_id: " << dump_id << endl;
		}

		begin_period_id = dump_id - 1;
		end_period_id = dump_id;

		for (int period_id = begin_period_id; period_id < end_period_id; period_id++)
		{
			pb->one_period_obs_deep_lpn(ad, cb, period_id);
		}

		if (dump_evo_avg == 1)
		{
			dump_adr_avg(ad, true);
		}
	}

	cb->dump_std(ad);
	cb->dump_lpn(ad);

	if (dump_evo_sep == 1)
	{
		cb->dump_std_evo(ad);
		cb->dump_lpn_evo(ad);
	}

	if (dump_evo_avg == 0)
	{
		dump_adr_avg(ad, true);
	}
}





void LpnDeepPer1TExperimentBehaviour::trans_process(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const
{
	ConfigParam * cp = ad->cp;

	int num_trajectories = cp->num_trajectories;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);
	int dump_evo_avg = int(cp->params.find("dump_evo_avg")->second);

#pragma omp parallel for
	for (int tr_id = 0; tr_id < 1; tr_id++)
	{
		int thread_id = omp_get_thread_num();
		trans_process_single_deep(ad, pb, cb, tr_id, thread_id);
	}

	cb->calc_chars_std_start(ad, 0);
	cb->calc_chars_lpn_start(ad, 0, 0);

#pragma omp parallel for
	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		if (tr_id > 0)
		{
			copy_trajectory_lpn(ad, tr_id, 0);
			var_trajectory_lpn(ad, cb, tr_id, 0);
		}

		resresh_times(ad, tr_id);

		cb->calc_chars_std_start(ad, tr_id);
		cb->calc_chars_lpn_start(ad, tr_id, 0);

		cb->evo_chars_std(ad, tr_id, 0);
		cb->evo_chars_lpn(ad, tr_id, 0);

		if (dump_evo_sep == 1)
		{
			dump_adr_single(ad, tr_id, false);
		}
	}

	if (dump_evo_avg == 1)
	{
		dump_adr_avg(ad, false);
	}
}

void LpnDeepPer1TExperimentBehaviour::obser_process(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	ExpData * ed = ad->ed;
	MainData * md = ad->md;

	ed->is_obs = 1;

	int num_branches = md->num_ham_qj;
	int num_sub_steps = num_branches * int(cp->params.find("deep_num_steps")->second);

	int num_trajectories = cp->num_trajectories;

	int dump_num_total = cp->num_obs_periods + 1;

	int * dump_periods = ed->dump_periods;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);
	int dump_evo_avg = int(cp->params.find("dump_evo_avg")->second);

	int num_lambdas_periods = int(cp->params.find("num_lambdas_periods")->second);

	int begin_period_id = 0;
	int end_period_id = 0;
	for (int dump_id = 1; dump_id < dump_num_total; dump_id++)
	{
		if (rp->is_pp == 1)
		{
			cout << "dump_id: " << dump_id << endl;
		}

		begin_period_id = dump_id - 1;
		end_period_id = dump_id;

		for (int period_id = begin_period_id; period_id < end_period_id; period_id++)
		{
			pb->one_period_obs_deep_lpn_per_period(ad, cb, period_id, num_lambdas_periods);
		}

		if (dump_evo_avg == 1)
		{
			dump_adr_avg(ad, true);
		}
	}

	cb->dump_std(ad);
	cb->dump_lpn(ad);

	if (dump_evo_sep == 1)
	{
		cb->dump_std_evo(ad);
		cb->dump_lpn_evo(ad);
	}

	if (dump_evo_avg == 0)
	{
		dump_adr_avg(ad, true);
	}
}





void LpnAllExperimentBehaviour::trans_process(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const
{
	ConfigParam * cp = ad->cp;

	int num_trajectories = cp->num_trajectories;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);
	int dump_evo_avg = int(cp->params.find("dump_evo_avg")->second);

#pragma omp parallel for
	for (int tr_id = 0; tr_id < 1; tr_id++)
	{
		int thread_id = omp_get_thread_num();
		trans_process_single(ad, pb, cb, tr_id, thread_id);
	}

	cb->calc_chars_std_start(ad, 0);
	cb->calc_chars_lpn_start(ad, 0, 0);

#pragma omp parallel for
	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		if (tr_id > 0)
		{
			copy_trajectory_lpn(ad, tr_id, 0);
		}
	}

	gs_orth_init(ad, cb);

#pragma omp parallel for
	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		resresh_times(ad, tr_id);

		cb->calc_chars_std_start(ad, tr_id);
		cb->calc_chars_lpn_start(ad, tr_id, 0);

		cb->evo_chars_std(ad, tr_id, 0);
		cb->evo_chars_lpn(ad, tr_id, 0);

		if (dump_evo_sep == 1)
		{
			dump_adr_single(ad, tr_id, false);
		}
	}

	if (dump_evo_avg == 1)
	{
		dump_adr_avg(ad, false);
	}
}

void LpnAllExperimentBehaviour::obser_process(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	ExpData * ed = ad->ed;
	MainData * md = ad->md;

	ed->is_obs = 1;

	int num_trajectories = cp->num_trajectories;

	int dump_num_total = ed->dump_num_total;
	int * dump_periods = ed->dump_periods;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);
	int dump_evo_avg = int(cp->params.find("dump_evo_avg")->second);

	int begin_period_id = 0;
	int end_period_id = 0;
	for (int dump_id = 1; dump_id < dump_num_total; dump_id++)
	{
		if (rp->is_pp == 1)
		{
			cout << "dump_id: " << dump_id << endl;
		}

		begin_period_id = dump_periods[dump_id - 1];
		end_period_id = dump_periods[dump_id];

		for (int period_id = begin_period_id; period_id < end_period_id; period_id++)
		{
#pragma omp parallel for
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				int thread_id = omp_get_thread_num();
				pb->one_period(ad, cb, tr_id, thread_id, period_id);
				cb->calc_chars_std(ad, tr_id);
			}

			cb->calc_chars_lpn(ad, 0, 0);

			ed->curr_time = (period_id + 1) * cb->calc_T(ad);

			lambda_lpn_all(ad, cb);
		}

#pragma omp parallel for
		for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
		{
			cb->evo_chars_std(ad, tr_id, dump_id);
			cb->evo_chars_lpn(ad, tr_id, dump_id);

			if (dump_evo_sep == 1)
			{
				dump_adr_single(ad, tr_id, true);
			}
		}

		if (dump_evo_avg == 1)
		{
			dump_adr_avg(ad, true);
		}
	}

	cb->dump_std(ad);
	cb->dump_lpn(ad);

	if (dump_evo_sep == 1)
	{
		cb->dump_std_evo(ad);
		cb->dump_lpn_evo(ad);
	}

	if (dump_evo_avg == 0)
	{
		dump_adr_avg(ad, true);
	}
}

MKL_Complex16 mult_scalar_double(MKL_Complex16 a, double b)
{
	MKL_Complex16 res = { 0.0, 0.0 };
	res.real = a.real * b;
	res.imag = a.imag * b;
	return res;
}

MKL_Complex16 mult_scalar_complex(MKL_Complex16 * a, MKL_Complex16 * b, int N)
{
	MKL_Complex16 res = { 0.0, 0.0 };
	for (int i = 0; i < N; i++)
	{
		res.real += a[i].real * b[i].real + a[i].imag * b[i].imag;
		res.imag += b[i].real * a[i].imag - a[i].real * b[i].imag;
	}
	return res;
}

int is_norm_crossed(MKL_Complex16 * phi, double * eta, int sys_size)
{
	double norm = 0.0;
	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		norm += phi[st_id].real * phi[st_id].real + phi[st_id].imag * phi[st_id].imag;
	}

	if (norm > 1.0)
	{
		cout << "norm more than 1" << endl;
		throw std::invalid_argument("received negative value");
	}

	if (norm > *(eta))
	{
		return 0;
	}
	else
	{
		return 1;
	}
}

double norm_square(MKL_Complex16 * phi, int sys_size)
{
	double norm = 0.0;
	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		norm += phi[st_id].real * phi[st_id].real + phi[st_id].imag * phi[st_id].imag;
	}
	return norm;
}

void recovery(AllData * ad, Split * head, int tr_id)
{
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int jump = int(ad->cp->params.find("jump")->second);

	MKL_Complex16 ZERO = { 0.0, 0.0 };
	MKL_Complex16 ONE = { 1.0, 0.0 };

	int sys_size = md->sys_size;
	VSLStreamStatePtr * stream = &(ed->streams[tr_id]);
	MKL_Complex16 * phi = &(ed->phi_all_aux[tr_id * sys_size]);
	double * eta = &(ed->etas_all[tr_id]);
	double * g = head->g;
	MKL_Complex16 * A = head->matrix;
	int k = head->steps;

	MKL_Complex16 * res = new MKL_Complex16[sys_size];
	double norm = sqrt(norm_square(phi, sys_size));
	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		phi[st_id].real /= (norm);
		phi[st_id].imag /= (norm);
	}
	norm = 0.0;

	double * gnorms = new double[k];
	double tmp = 0.0;
	double ran = 0.0;
	vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, *stream, 1, &ran, 0.0, 1.0);

	for (int i = 0; i < k; i++)
	{
		cblas_zgemv(
			CblasRowMajor,
			CblasNoTrans,
			sys_size,
			sys_size,
			&ONE,
			&A[i *  sys_size * sys_size + 0],
			sys_size,
			phi,
			1,
			&ZERO,
			res,
			1
		);

		gnorms[i] = norm_square(res, sys_size);
		gnorms[i] *= g[i];
		tmp += gnorms[i];
	}

	ran *= tmp;

	int index = 0;
	while (ran - gnorms[index] > 0.0)
	{
		ran -= gnorms[index];
		index++;
		if (index == k - 1)
		{
			break;
		}
	}

	while (gnorms[index] == 0)
	{
		if (index == 0)
		{
			index++;
		}
		else
		{
			index--;
		}
	}

	memset(res, 0, sys_size * sizeof(MKL_Complex16));

	cblas_zgemv(
		CblasRowMajor,
		CblasNoTrans,
		sys_size,
		sys_size,
		&ONE,
		&A[index *  sys_size * sys_size + 0],
		sys_size,
		phi,
		1,
		&ZERO,
		res,
		1
	);

	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		phi[st_id].real = res[st_id].real / sqrt(gnorms[index] / g[index]);
		phi[st_id].imag = res[st_id].imag / sqrt(gnorms[index] / g[index]);
	}

	if (jump > 0 && ed->is_obs == 1)
	{
		ed->diss_types[tr_id].push_back(index);
	}

	delete[] res;
	delete[] gnorms;
}

double get_norm_cd(double * vec, int size)
{
	double sum = 0.0;
	for (int i = 0; i < size; i++)
	{
		sum += vec[i] * vec[i];
	}

	return sqrt(sum);
}

double get_mean_simple(double * adr, int sys_size)
{
	double mean = 0.0;
	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		mean += double(st_id) * adr[st_id];
	}

	return mean;
}

double get_dispersion_simple(double mean_curr, double mean_start)
{
	double dispersion = (mean_curr - mean_start) * (mean_curr - mean_start);
	return dispersion;
}

double get_m2(double * adr, int sys_size, double mean)
{
	double m2 = 0.0;
	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		m2 += (double(st_id) - mean) * (double(st_id) - mean) * adr[st_id];
	}

	return m2;
}

double get_energy(AllData * ad, int tr_id)
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;

	double dimer_prm_E = double(cp->params.find("dimer_prm_E")->second);

	MKL_Complex16 * phi = &(ed->phi_all[tr_id * sys_size]);
	double * hamiltonian = md->hamiltonian;
	double * hamiltonian_drv = md->hamiltonian_drv;

	double norm_2 = norm_square(phi, sys_size);

	MKL_Complex16 * sm = new MKL_Complex16[sys_size];

	for (int st_id_1 = 0; st_id_1 < sys_size; st_id_1++)
	{
		MKL_Complex16 tmp;
		tmp.real = 0;
		tmp.imag = 0;
		for (int st_id_2 = 0; st_id_2 < sys_size; st_id_2++)
		{
			int index = st_id_1 * sys_size + st_id_2;

			double ham_val = (hamiltonian[index] + dimer_prm_E * hamiltonian_drv[index]);

			tmp.real += (ham_val * phi[st_id_2].real / sqrt(norm_2) - (0.0) * phi[st_id_2].imag / sqrt(norm_2));
			tmp.imag += (ham_val * phi[st_id_2].imag / sqrt(norm_2) + (0.0) * phi[st_id_2].real / sqrt(norm_2));
		}

		sm[st_id_1].real = tmp.real;
		sm[st_id_1].imag = tmp.imag;
	}

	double energy = 0.0;
	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		energy += (phi[st_id].real / sqrt(norm_2) * sm[st_id].real + phi[st_id].imag / sqrt(norm_2) * sm[st_id].imag);
	}

	delete[] sm;

	return energy;
}

MKL_Complex16 get_spec_jcs(AllData * ad, int tr_id)
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;

	double alpha = double(cp->params.find("jcs_prm_alpha")->second);

	MKL_Complex16 * phi = &(ed->phi_all[tr_id * sys_size]);
	MKL_Complex16 * phi_normed = new MKL_Complex16[sys_size];
	MKL_Complex16 * phi_normed_conj = new MKL_Complex16[sys_size];
	double norm = sqrt(norm_square(phi, sys_size));
	MKL_Complex16 * mult_tmp = new MKL_Complex16[sys_size];
	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		mult_tmp[st_id].real = 0.0;
		mult_tmp[st_id].imag = 0.0;

		phi_normed[st_id].real = phi[st_id].real / norm;
		phi_normed[st_id].imag = phi[st_id].imag / norm;

		phi_normed_conj[st_id].real = phi_normed[st_id].real;
		phi_normed_conj[st_id].imag = -phi_normed[st_id].imag;
	}

	MKL_Complex16 ZERO = { 0.0, 0.0 };
	MKL_Complex16 ONE = { 1.0, 0.0 };
	cblas_zgemv(CblasRowMajor, CblasNoTrans, sys_size, sys_size, &ONE, md->special, sys_size, phi_normed, 1, &ZERO, mult_tmp, 1);

	MKL_Complex16 result = { 0.0, 0.0 };
	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		result.real += (phi_normed_conj[st_id].real * mult_tmp[st_id].real - phi_normed_conj[st_id].imag * mult_tmp[st_id].imag);
		result.imag += (phi_normed_conj[st_id].imag * mult_tmp[st_id].real + phi_normed_conj[st_id].real * mult_tmp[st_id].imag);
	}

	result.real /= alpha;
	result.imag /= alpha;

	delete[] mult_tmp;
	delete[] phi_normed;
	delete[] phi_normed_conj;

	return result;
}

MKL_Complex16 get_spec_ps(AllData * ad, int tr_id)
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;

	double alpha = double(cp->params.find("ps_prm_alpha")->second);

	MKL_Complex16 * phi = &(ed->phi_all[tr_id * sys_size]);
	MKL_Complex16 * phi_normed = new MKL_Complex16[sys_size];
	MKL_Complex16 * phi_normed_conj = new MKL_Complex16[sys_size];
	double norm = sqrt(norm_square(phi, sys_size));
	MKL_Complex16 * mult_tmp = new MKL_Complex16[sys_size];
	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		mult_tmp[st_id].real = 0.0;
		mult_tmp[st_id].imag = 0.0;

		phi_normed[st_id].real = phi[st_id].real / norm;
		phi_normed[st_id].imag = phi[st_id].imag / norm;

		phi_normed_conj[st_id].real = phi_normed[st_id].real;
		phi_normed_conj[st_id].imag = -phi_normed[st_id].imag;
	}

	MKL_Complex16 ZERO = { 0.0, 0.0 };
	MKL_Complex16 ONE = { 1.0, 0.0 };
	cblas_zgemv(CblasRowMajor, CblasNoTrans, sys_size, sys_size, &ONE, md->special, sys_size, phi_normed, 1, &ZERO, mult_tmp, 1);

	MKL_Complex16 result = { 0.0, 0.0 };
	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		result.real += (phi_normed_conj[st_id].real * mult_tmp[st_id].real - phi_normed_conj[st_id].imag * mult_tmp[st_id].imag);
		result.imag += (phi_normed_conj[st_id].imag * mult_tmp[st_id].real + phi_normed_conj[st_id].real * mult_tmp[st_id].imag);
	}

	result.real /= alpha;
	result.imag /= alpha;

	delete[] mult_tmp;
	delete[] phi_normed;
	delete[] phi_normed_conj;

	return result;
}

MKL_Complex16 get_spec_2_ps(AllData * ad, int tr_id)
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;

	double alpha = double(cp->params.find("ps_prm_alpha")->second);

	MKL_Complex16 * phi = &(ed->phi_all[tr_id * sys_size]);
	MKL_Complex16 * phi_normed = new MKL_Complex16[sys_size];
	MKL_Complex16 * phi_normed_conj = new MKL_Complex16[sys_size];
	double norm = sqrt(norm_square(phi, sys_size));
	MKL_Complex16 * mult_tmp = new MKL_Complex16[sys_size];
	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		mult_tmp[st_id].real = 0.0;
		mult_tmp[st_id].imag = 0.0;

		phi_normed[st_id].real = phi[st_id].real / norm;
		phi_normed[st_id].imag = phi[st_id].imag / norm;

		phi_normed_conj[st_id].real = phi_normed[st_id].real;
		phi_normed_conj[st_id].imag = -phi_normed[st_id].imag;
	}

	MKL_Complex16 ZERO = { 0.0, 0.0 };
	MKL_Complex16 ONE = { 1.0, 0.0 };
	cblas_zgemv(CblasRowMajor, CblasNoTrans, sys_size, sys_size, &ONE, md->special_2, sys_size, phi_normed, 1, &ZERO, mult_tmp, 1);

	MKL_Complex16 result = { 0.0, 0.0 };
	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		result.real += (phi_normed_conj[st_id].real * mult_tmp[st_id].real - phi_normed_conj[st_id].imag * mult_tmp[st_id].imag);
		result.imag += (phi_normed_conj[st_id].imag * mult_tmp[st_id].real + phi_normed_conj[st_id].real * mult_tmp[st_id].imag);
	}

	delete[] mult_tmp;
	delete[] phi_normed;
	delete[] phi_normed_conj;

	return result;
}


MKL_Complex16 get_spec_3_ps(AllData * ad, int tr_id)
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;

	MKL_Complex16 * phi = &(ed->phi_all[tr_id * sys_size]);
	MKL_Complex16 * phi_normed = new MKL_Complex16[sys_size];
	MKL_Complex16 * phi_normed_conj = new MKL_Complex16[sys_size];
	double norm = sqrt(norm_square(phi, sys_size));
	MKL_Complex16 * mult_tmp = new MKL_Complex16[sys_size];
	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		mult_tmp[st_id].real = 0.0;
		mult_tmp[st_id].imag = 0.0;

		phi_normed[st_id].real = phi[st_id].real / norm;
		phi_normed[st_id].imag = phi[st_id].imag / norm;

		phi_normed_conj[st_id].real = phi_normed[st_id].real;
		phi_normed_conj[st_id].imag = -phi_normed[st_id].imag;
	}

	MKL_Complex16 ZERO = { 0.0, 0.0 };
	MKL_Complex16 ONE = { 1.0, 0.0 };
	cblas_zgemv(CblasRowMajor, CblasNoTrans, sys_size, sys_size, &ONE, md->special_3, sys_size, phi_normed, 1, &ZERO, mult_tmp, 1);

	MKL_Complex16 result = { 0.0, 0.0 };
	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		result.real += (phi_normed_conj[st_id].real * mult_tmp[st_id].real - phi_normed_conj[st_id].imag * mult_tmp[st_id].imag);
		result.imag += (phi_normed_conj[st_id].imag * mult_tmp[st_id].real + phi_normed_conj[st_id].real * mult_tmp[st_id].imag);
	}

	delete[] mult_tmp;
	delete[] phi_normed;
	delete[] phi_normed_conj;

	return result;
}

MKL_Complex16 get_spec_mbl(AllData * ad, int tr_id)
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;

	MKL_Complex16 * phi = &(ed->phi_all[tr_id * sys_size]);
	MKL_Complex16 * phi_normed = new MKL_Complex16[sys_size];
	MKL_Complex16 * phi_normed_conj = new MKL_Complex16[sys_size];
	double norm = sqrt(norm_square(phi, sys_size));
	MKL_Complex16 * mult_tmp = new MKL_Complex16[sys_size];
	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		mult_tmp[st_id].real = 0.0;
		mult_tmp[st_id].imag = 0.0;

		phi_normed[st_id].real = phi[st_id].real / norm;
		phi_normed[st_id].imag = phi[st_id].imag / norm;

		phi_normed_conj[st_id].real = phi_normed[st_id].real;
		phi_normed_conj[st_id].imag = -phi_normed[st_id].imag;
	}

	MKL_Complex16 ZERO = { 0.0, 0.0 };
	MKL_Complex16 ONE = { 1.0, 0.0 };
	cblas_zgemv(CblasRowMajor, CblasNoTrans, sys_size, sys_size, &ONE, md->special, sys_size, phi_normed, 1, &ZERO, mult_tmp, 1);

	MKL_Complex16 result = { 0.0, 0.0 };
	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		result.real += (phi_normed_conj[st_id].real * mult_tmp[st_id].real - phi_normed_conj[st_id].imag * mult_tmp[st_id].imag);
		result.imag += (phi_normed_conj[st_id].imag * mult_tmp[st_id].real + phi_normed_conj[st_id].real * mult_tmp[st_id].imag);
	}

	delete[] mult_tmp;
	delete[] phi_normed;
	delete[] phi_normed_conj;

	return result;
}

double get_imbalance_mbl(AllData * ad, double * adr)
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;

	double * n_part = new double[md->mbl_Nc];

	for (int cell_id = 0; cell_id < md->mbl_Nc; cell_id++)
	{
		n_part[cell_id] = 0.0;
	}

	for (int state_id = 0; state_id < sys_size; state_id++)
	{
		vector<int> vb = convert_int_to_vector_of_bits(md->mbl_id_to_x[state_id], md->mbl_Nc);

		for (int cell_id = 0; cell_id < md->mbl_Nc; cell_id++)
		{
			n_part[cell_id] += adr[state_id] * double(vb[cell_id]);
		}
	}

	double sum_odd = 0.0;
	double sum_even = 0.0;
	double sum_all = 0.0;

	for (int cell_id = 0; cell_id < md->mbl_Nc; cell_id++)
	{
		if (cell_id % 2 == 0)
		{
			sum_odd += n_part[cell_id];
		}
		else
		{
			sum_even += n_part[cell_id];
		}
	}
	sum_all = sum_even + sum_odd;

	double imbalance = (sum_odd - sum_even) / sum_all;

	delete[] n_part;

	return imbalance;
}


vector<complex<double>> get_random_obs(AllData * ad, int tr_id)
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int num_random_obs = int(cp->params.find("num_random_obs")->second);

	vector<complex<double>> random_obs;

	if (num_random_obs > 0)
	{ 
		int sys_size = md->sys_size;

		MKL_Complex16 * phi = &(ed->phi_all[tr_id * sys_size]);
		MKL_Complex16 * phi_normed = new MKL_Complex16[sys_size];
		MKL_Complex16 * phi_normed_conj = new MKL_Complex16[sys_size];
		double norm = sqrt(norm_square(phi, sys_size));
		MKL_Complex16 * mult_tmp = new MKL_Complex16[sys_size];
		for (int st_id = 0; st_id < sys_size; st_id++)
		{
			phi_normed[st_id].real = phi[st_id].real / norm;
			phi_normed[st_id].imag = phi[st_id].imag / norm;

			phi_normed_conj[st_id].real = phi_normed[st_id].real;
			phi_normed_conj[st_id].imag = -phi_normed[st_id].imag;
		}

		MKL_Complex16 ZERO = { 0.0, 0.0 };
		MKL_Complex16 ONE = { 1.0, 0.0 };
		for (int obs_id = 0; obs_id < num_random_obs; obs_id++)
		{
			cblas_zgemv(CblasRowMajor, CblasNoTrans, sys_size, sys_size, &ONE, md->random_obs_mtxs[obs_id], sys_size, phi_normed, 1, &ZERO, mult_tmp, 1);
			complex<double> result(0.0, 0.0);

			for (int st_id = 0; st_id < sys_size; st_id++)
			{
				complex<double> tmp(
					(phi_normed_conj[st_id].real * mult_tmp[st_id].real - phi_normed_conj[st_id].imag * mult_tmp[st_id].imag),
					(phi_normed_conj[st_id].imag * mult_tmp[st_id].real + phi_normed_conj[st_id].real * mult_tmp[st_id].imag)

				);

				result += tmp;
			}

			random_obs.push_back(result);
		}

		delete[] mult_tmp;
		delete[] phi_normed;
		delete[] phi_normed_conj;
	}

	return random_obs;
}

MKL_Complex16 get_num_photons_jcs(AllData * ad, int tr_id)
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;

	double alpha = double(cp->params.find("jcs_prm_alpha")->second);

	MKL_Complex16 * phi = &(ed->phi_all[tr_id * sys_size]);
	MKL_Complex16 * phi_normed = new MKL_Complex16[sys_size];
	MKL_Complex16 * phi_normed_conj = new MKL_Complex16[sys_size];
	double norm = sqrt(norm_square(phi, sys_size));
	MKL_Complex16 * mult_tmp = new MKL_Complex16[sys_size];

	double * a_std = new double[md->sys_size * md->sys_size];
	double * a_dag = new double[md->sys_size * md->sys_size];
	double * n_mtx_double = new double[md->sys_size * md->sys_size];
	MKL_Complex16 * n_mtx = new MKL_Complex16[md->sys_size * md->sys_size];

	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		mult_tmp[st_id].real = 0.0;
		mult_tmp[st_id].imag = 0.0;

		phi_normed[st_id].real = phi[st_id].real / norm;
		phi_normed[st_id].imag = phi[st_id].imag / norm;

		phi_normed_conj[st_id].real = phi_normed[st_id].real;
		phi_normed_conj[st_id].imag = -phi_normed[st_id].imag;
	}

	for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
	{
		for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
		{
			int index = st_id_1 * md->sys_size + st_id_2;

			a_std[index] = 0.0;
			a_dag[index] = 0.0;
			n_mtx_double[index] = 0.0;
			n_mtx[index].real = 0.0;
			n_mtx[index].imag = 0.0;
		}
	}

	for (int st_id = 0; st_id < sys_size - 1; st_id++)
	{
		int index_std = st_id * md->sys_size + (st_id + 1);
		int index_dag = (st_id + 1) * md->sys_size + st_id;
		a_std[index_std] = sqrt(double(st_id + 1));
		a_dag[index_dag] = sqrt(double(st_id + 1));
	}

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, sys_size, sys_size, sys_size, 1.0, a_dag, sys_size, a_std, sys_size, 0.0, n_mtx_double, sys_size);
	for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
	{
		for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
		{
			int index = st_id_1 * md->sys_size + st_id_2;

			n_mtx[index].real = n_mtx_double[index];
			n_mtx[index].imag = 0.0;
		}
	}

	MKL_Complex16 ZERO = { 0.0, 0.0 };
	MKL_Complex16 ONE = { 1.0, 0.0 };
	cblas_zgemv(CblasRowMajor, CblasNoTrans, sys_size, sys_size, &ONE, n_mtx, sys_size, phi_normed, 1, &ZERO, mult_tmp, 1);

	MKL_Complex16 result = { 0.0, 0.0 };
	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		result.real += (phi_normed_conj[st_id].real * mult_tmp[st_id].real - phi_normed_conj[st_id].imag * mult_tmp[st_id].imag);
		result.imag += (phi_normed_conj[st_id].imag * mult_tmp[st_id].real + phi_normed_conj[st_id].real * mult_tmp[st_id].imag);
	}

	delete[] mult_tmp;
	delete[] phi_normed;
	delete[] phi_normed_conj;

	delete[] a_std;
	delete[] a_dag;
	delete[] n_mtx_double;
	delete[] n_mtx;

	return result;
}

MKL_Complex16 get_num_photons_ps(AllData * ad, int tr_id)
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;

	MKL_Complex16 * phi = &(ed->phi_all[tr_id * sys_size]);
	MKL_Complex16 * phi_normed = new MKL_Complex16[sys_size];
	MKL_Complex16 * phi_normed_conj = new MKL_Complex16[sys_size];
	double norm = sqrt(norm_square(phi, sys_size));
	MKL_Complex16 * mult_tmp = new MKL_Complex16[sys_size];

	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		mult_tmp[st_id].real = 0.0;
		mult_tmp[st_id].imag = 0.0;

		phi_normed[st_id].real = phi[st_id].real / norm;
		phi_normed[st_id].imag = phi[st_id].imag / norm;

		phi_normed_conj[st_id].real = phi_normed[st_id].real;
		phi_normed_conj[st_id].imag = -phi_normed[st_id].imag;
	}

	MKL_Complex16 ZERO = { 0.0, 0.0 };
	MKL_Complex16 ONE = { 1.0, 0.0 };
	cblas_zgemv(CblasRowMajor, CblasNoTrans, sys_size, sys_size, &ONE, md->special_4, sys_size, phi_normed, 1, &ZERO, mult_tmp, 1);

	MKL_Complex16 result = { 0.0, 0.0 };
	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		result.real += (phi_normed_conj[st_id].real * mult_tmp[st_id].real - phi_normed_conj[st_id].imag * mult_tmp[st_id].imag);
		result.imag += (phi_normed_conj[st_id].imag * mult_tmp[st_id].real + phi_normed_conj[st_id].real * mult_tmp[st_id].imag);
	}

	delete[] mult_tmp;
	delete[] phi_normed;
	delete[] phi_normed_conj;

	return result;
}

void resresh_times(AllData * ad, int tr_id)
{
	ExpData * ed = ad->ed;

	ed->times_all[tr_id] = 0.0;
}

void copy_trajectory_lpn(AllData * ad, int tr_id, int base_tr_id)
{
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;

	VSLStreamStatePtr * streams = ed->streams;
	MKL_Complex16 * phi = &(ed->phi_all[tr_id * sys_size]);
	MKL_Complex16 * phi_original = &(ed->phi_all[base_tr_id * sys_size]);
	double * etas = ed->etas_all;

	vslCopyStream(&streams[tr_id], streams[base_tr_id]);

	etas[tr_id] = etas[base_tr_id];

	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		phi[st_id].real = phi_original[st_id].real;
		phi[st_id].imag = phi_original[st_id].imag;
	}
}

void copy_stream_lpn(AllData * ad, int tr_id, int base_tr_id)
{
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	VSLStreamStatePtr * streams = ed->streams;
	double * etas = ed->etas_all;

	vslCopyStream(&streams[tr_id], streams[base_tr_id]);

	etas[tr_id] = etas[base_tr_id];
}

void copy_trajectory_data(AllData * ad, int tr_id, int base_tr_id)
{
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;

	MKL_Complex16 * phi = &(ed->phi_all[tr_id * sys_size]);
	MKL_Complex16 * phi_original = &(ed->phi_all[base_tr_id]);

	double * etas = ed->etas_all;

	etas[tr_id] = etas[base_tr_id];

	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		phi[st_id].real = phi_original[st_id].real;
		phi[st_id].imag = phi_original[st_id].imag;
	}
}

void var_trajectory_lpn(AllData * ad, CoreBehavior *cb, int tr_id, int base_tr_id)
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;

	double lpn_eps_low = double(cp->params.find("lpn_eps_low")->second);
	double lpn_eps_high = double(cp->params.find("lpn_eps_high")->second);
	double lpn_eps_error = double(cp->params.find("lpn_eps_error")->second);
	int lpn_eps_deep = int(cp->params.find("lpn_eps_deep")->second);

	double curr_eps_low = lpn_eps_low;
	double curr_eps_high = lpn_eps_high;
	double curr_eps_diff = curr_eps_high - curr_eps_low;
	int curr_eps_deep = 0;
	double lpn_eps_start = pow(10.0, (curr_eps_low + curr_eps_high) * 0.5);
	double curr_lpn_eps = lpn_eps_start;

	double lpn_delta_s = double(cp->params.find("lpn_delta_s")->second);

	VSLStreamStatePtr * streams_var = ed->streams_var;
	MKL_Complex16 * phi_original = &(ed->phi_all[base_tr_id * sys_size]);
	MKL_Complex16 * phi = &(ed->phi_all[tr_id * sys_size]);

	MKL_Complex16 * phi_copy = new MKL_Complex16[sys_size];
	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		phi_copy[st_id].real = phi[st_id].real;
		phi_copy[st_id].imag = phi[st_id].imag;
	}

	MKL_Complex16 * phi_var = new MKL_Complex16[sys_size];
	double * phi_var_double = new double[2 * sys_size];


	vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, streams_var[tr_id], 2 * sys_size, phi_var_double, -1.0, 1.0);
	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		phi_var[st_id].real = phi_var_double[0 * sys_size + st_id];
		phi_var[st_id].imag = phi_var_double[1 * sys_size + st_id];
	}

	double norm_var_2 = norm_square(phi_var, sys_size);
	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		phi_var[st_id].real = phi_var[st_id].real / sqrt(norm_var_2);
		phi_var[st_id].imag = phi_var[st_id].imag / sqrt(norm_var_2);
	}

	double norm_2 = norm_square(phi_original, sys_size);

	while (curr_eps_deep < lpn_eps_deep)
	{
		for (int st_id = 0; st_id < sys_size; st_id++)
		{
			phi[st_id].real = phi_copy[st_id].real;
			phi[st_id].imag = phi_copy[st_id].imag;
		}

		for (int st_id = 0; st_id < sys_size; st_id++)
		{
			phi[st_id].real = phi_original[st_id].real + curr_lpn_eps * phi_var[st_id].real;
			phi[st_id].imag = phi_original[st_id].imag + curr_lpn_eps * phi_var[st_id].imag;
		}

		double norm_mod_2 = norm_square(phi, sys_size);
		for (int st_id = 0; st_id < sys_size; st_id++)
		{
			phi[st_id].real /= sqrt(norm_mod_2);
			phi[st_id].imag /= sqrt(norm_mod_2);
		}

		for (int st_id = 0; st_id < sys_size; st_id++)
		{
			phi[st_id].real *= sqrt(norm_2);
			phi[st_id].imag *= sqrt(norm_2);
		}

		cb->calc_chars_lpn(ad, tr_id, base_tr_id);

		double delta_s = cb->calc_delta_s(ad, tr_id, base_tr_id);

		if (fabs(delta_s - lpn_delta_s) < lpn_eps_error)
		{
			break;
		}

		if (delta_s > lpn_delta_s)
		{
			curr_eps_high -= 0.5 * curr_eps_diff;
			curr_eps_diff = curr_eps_high - curr_eps_low;
			curr_eps_deep++;
			curr_lpn_eps = pow(10.0, (curr_eps_low + curr_eps_high) * 0.5);
		}

		if (delta_s < lpn_delta_s)
		{
			curr_eps_low += 0.5 * curr_eps_diff;
			curr_eps_diff = curr_eps_high - curr_eps_low;
			curr_eps_deep++;
			curr_lpn_eps = pow(10.0, (curr_eps_low + curr_eps_high) * 0.5);
		}
	}

	if (rp->is_pp == 1)
	{
		cout << "lpn: " << tr_id + 1 << endl;
		cout << "curr_lpn_eps: " << log10(curr_lpn_eps) << endl;
		cout << "lpn_deep: " << curr_eps_deep << endl;
		cout << "delta_s: " << cb->calc_delta_s(ad, tr_id, base_tr_id) << endl;
		cout << endl;
	}

	delete[] phi_var_double;
	delete[] phi_var;
	delete[] phi_copy;
}

void var_first(
	MKL_Complex16 * phi_var,
	double * phi_var_double,
	MKL_Complex16 * phi_var_all,
	MKL_Complex16 * scalar_mults_all,
	AllData * ad,
	CoreBehavior *cb
)
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;
	int num_trajectories = cp->num_trajectories;

	VSLStreamStatePtr * streams_var = ed->streams_var;

	vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, streams_var[1], 2 * sys_size, phi_var_double, -1.0, 1.0);
	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		phi_var[st_id].real = phi_var_double[0 * sys_size + st_id];
		phi_var[st_id].imag = phi_var_double[1 * sys_size + st_id];
	}

	double norm_var_2 = norm_square(phi_var, sys_size);
	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		phi_var[st_id].real = phi_var[st_id].real / sqrt(norm_var_2);
		phi_var[st_id].imag = phi_var[st_id].imag / sqrt(norm_var_2);

		int index = 0 * sys_size + st_id;
		phi_var_all[index].real = phi_var[st_id].real;
		phi_var_all[index].imag = phi_var[st_id].imag;
	}
}

void var_first_with_history(
	MKL_Complex16 * phi_var,
	double * phi_var_double,
	MKL_Complex16 * phi_var_all,
	MKL_Complex16 * scalar_mults_all,
	AllData * ad,
	CoreBehavior *cb
)
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;

	MKL_Complex16 * phi_original = &(ed->phi_all[0]);
	MKL_Complex16 * phi = &(ed->phi_all[1 * sys_size]);

	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		phi_var[st_id].real = phi[st_id].real - phi_original[st_id].real;
		phi_var[st_id].imag = phi[st_id].imag - phi_original[st_id].imag;
	}

	double norm_var_2 = norm_square(phi_var, sys_size);
	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		phi_var[st_id].real = phi_var[st_id].real / sqrt(norm_var_2);
		phi_var[st_id].imag = phi_var[st_id].imag / sqrt(norm_var_2);

		int index = 0 * sys_size + st_id;
		phi_var_all[index].real = phi_var[st_id].real;
		phi_var_all[index].imag = phi_var[st_id].imag;
	}
}

void var_not_first(
	int tr_id,
	MKL_Complex16 * phi_var,
	double * phi_var_double,
	MKL_Complex16 * phi_var_all,
	MKL_Complex16 * scalar_mults_all,
	AllData * ad,
	CoreBehavior *cb
)
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;
	int num_trajectories = cp->num_trajectories;

	int lpn_id = tr_id - 1;
	
	VSLStreamStatePtr * streams_var = ed->streams_var;

	vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, streams_var[tr_id], 2 * sys_size, phi_var_double, -1.0, 1.0);
	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		phi_var[st_id].real = phi_var_double[0 * sys_size + st_id];
		phi_var[st_id].imag = phi_var_double[1 * sys_size + st_id];
	}

	double norm_var_2 = norm_square(phi_var, sys_size);
	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		phi_var[st_id].real = phi_var[st_id].real / sqrt(norm_var_2);
		phi_var[st_id].imag = phi_var[st_id].imag / sqrt(norm_var_2);

		int index = lpn_id * sys_size + st_id;
		phi_var_all[index].real = phi_var[st_id].real;
		phi_var_all[index].imag = phi_var[st_id].imag;
	}

	// orth
	for (int lpn_id_tmp = 0; lpn_id_tmp < lpn_id; lpn_id_tmp++)
	{
		scalar_mults_all[lpn_id_tmp].real = 0.0;
		scalar_mults_all[lpn_id_tmp].imag = 0.0;

		for (int lpn_st_id = 0; lpn_st_id < sys_size; lpn_st_id++)
		{
			int index = lpn_id * sys_size + lpn_st_id;
			int index_tmp = lpn_id_tmp * sys_size + lpn_st_id;
			scalar_mults_all[lpn_id_tmp].real += (phi_var_all[index].real * phi_var_all[index_tmp].real + phi_var_all[index].imag * phi_var_all[index_tmp].imag);
			scalar_mults_all[lpn_id_tmp].imag += (phi_var_all[index].imag * phi_var_all[index_tmp].real - phi_var_all[index].real * phi_var_all[index_tmp].imag);
		}

		for (int lpn_st_id = 0; lpn_st_id < sys_size; lpn_st_id++)
		{
			int index = lpn_id * sys_size + lpn_st_id;
			int index_tmp = lpn_id_tmp * sys_size + lpn_st_id;

			phi_var_all[index].real -= (scalar_mults_all[lpn_id_tmp].real * phi_var_all[index_tmp].real - scalar_mults_all[lpn_id_tmp].imag * phi_var_all[index_tmp].imag);
			phi_var_all[index].imag -= (scalar_mults_all[lpn_id_tmp].imag * phi_var_all[index_tmp].real + scalar_mults_all[lpn_id_tmp].real * phi_var_all[index_tmp].imag);
		}
	}

	// back orth to var
	double norm_test = norm_square(&phi_var_all[lpn_id * sys_size], sys_size);
	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		int index = lpn_id * sys_size + st_id;

		phi_var_all[index].real = phi_var_all[index].real / sqrt(norm_test);
		phi_var_all[index].imag = phi_var_all[index].imag / sqrt(norm_test);

		phi_var[st_id].real = phi_var_all[index].real;
		phi_var[st_id].imag = phi_var_all[index].imag;
	}
}

void var_not_first_with_history(
	int tr_id,
	MKL_Complex16 * phi_var,
	double * phi_var_double,
	MKL_Complex16 * phi_var_all,
	MKL_Complex16 * scalar_mults_all,
	AllData * ad,
	CoreBehavior *cb
)
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;
	int num_trajectories = cp->num_trajectories;

	int lpn_id = tr_id - 1;

	MKL_Complex16 * phi_original = &(ed->phi_all[0]);
	MKL_Complex16 * phi = &(ed->phi_all[tr_id * sys_size]);

	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		phi_var[st_id].real = phi[st_id].real - phi_original[st_id].real;
		phi_var[st_id].imag = phi[st_id].imag - phi_original[st_id].imag;
	}

	double norm_var_2 = norm_square(phi_var, sys_size);
	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		phi_var[st_id].real = phi_var[st_id].real / sqrt(norm_var_2);
		phi_var[st_id].imag = phi_var[st_id].imag / sqrt(norm_var_2);

		int index = lpn_id * sys_size + st_id;
		phi_var_all[index].real = phi_var[st_id].real;
		phi_var_all[index].imag = phi_var[st_id].imag;
	}

	// orth
	for (int lpn_id_tmp = 0; lpn_id_tmp < lpn_id; lpn_id_tmp++)
	{
		scalar_mults_all[lpn_id_tmp].real = 0.0;
		scalar_mults_all[lpn_id_tmp].imag = 0.0;

		for (int lpn_st_id = 0; lpn_st_id < sys_size; lpn_st_id++)
		{
			int index = lpn_id * sys_size + lpn_st_id;
			int index_tmp = lpn_id_tmp * sys_size + lpn_st_id;
			scalar_mults_all[lpn_id_tmp].real += (phi_var_all[index].real * phi_var_all[index_tmp].real + phi_var_all[index].imag * phi_var_all[index_tmp].imag);
			scalar_mults_all[lpn_id_tmp].imag += (phi_var_all[index].imag * phi_var_all[index_tmp].real - phi_var_all[index].real * phi_var_all[index_tmp].imag);
		}

		for (int lpn_st_id = 0; lpn_st_id < sys_size; lpn_st_id++)
		{
			int index = lpn_id * sys_size + lpn_st_id;
			int index_tmp = lpn_id_tmp * sys_size + lpn_st_id;

			phi_var_all[index].real -= (scalar_mults_all[lpn_id_tmp].real * phi_var_all[index_tmp].real - scalar_mults_all[lpn_id_tmp].imag * phi_var_all[index_tmp].imag);
			phi_var_all[index].imag -= (scalar_mults_all[lpn_id_tmp].imag * phi_var_all[index_tmp].real + scalar_mults_all[lpn_id_tmp].real * phi_var_all[index_tmp].imag);
		}
	}

	// back orth to var
	double norm_test = norm_square(&phi_var_all[lpn_id * sys_size], sys_size);
	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		int index = lpn_id * sys_size + st_id;

		phi_var_all[index].real = phi_var_all[index].real / sqrt(norm_test);
		phi_var_all[index].imag = phi_var_all[index].imag / sqrt(norm_test);

		phi_var[st_id].real = phi_var_all[index].real;
		phi_var[st_id].imag = phi_var_all[index].imag;
	}
}

void gs_orth_init(AllData * ad, CoreBehavior *cb)
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;
	int num_trajectories = cp->num_trajectories;

	int num_lpns = num_trajectories - 1;

	double lpn_eps_low = double(cp->params.find("lpn_eps_low")->second);
	double lpn_eps_high = double(cp->params.find("lpn_eps_high")->second);
	int lpn_eps_deep = int(cp->params.find("lpn_eps_deep")->second);

	double curr_eps_low = lpn_eps_low;
	double curr_eps_high = lpn_eps_high;
	double curr_eps_diff = lpn_eps_high - lpn_eps_low;
	int curr_eps_deep = 0;
	double lpn_eps_start = pow(10.0, (curr_eps_low + curr_eps_high) * 0.5);
	double curr_lpn_eps = lpn_eps_start;

	double lpn_delta_s_high = double(cp->params.find("lpn_delta_s_high")->second);
	double lpn_delta_s_low = double(cp->params.find("lpn_delta_s_low")->second);

	MKL_Complex16 * phi_copy = new MKL_Complex16[sys_size];

	MKL_Complex16 * phi_var = new MKL_Complex16[sys_size];
	double * phi_var_double = new double[2 * sys_size];
	MKL_Complex16 * phi_var_all = new MKL_Complex16[num_lpns * sys_size];
	MKL_Complex16 * scalar_mults_all = new MKL_Complex16[num_lpns];

	VSLStreamStatePtr * streams_var = ed->streams_var;
	MKL_Complex16 * phi_original = &(ed->phi_all[0]);
	MKL_Complex16 * phi = NULL;

	// ===== First lpn =====

	phi = &(ed->phi_all[1 * sys_size]);

	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		phi_copy[st_id].real = phi[st_id].real;
		phi_copy[st_id].imag = phi[st_id].imag;
	}

	double delta_s = lpn_delta_s_high + 1.0;

	var_first(phi_var, phi_var_double, phi_var_all, scalar_mults_all, ad, cb);

	while ((delta_s > lpn_delta_s_high || delta_s < lpn_delta_s_low) && (curr_eps_deep < lpn_eps_deep))
	{
		for (int st_id = 0; st_id < sys_size; st_id++)
		{
			phi[st_id].real = phi_copy[st_id].real;
			phi[st_id].imag = phi_copy[st_id].imag;
		}

		double norm_2 = norm_square(phi_original, sys_size);
		for (int st_id = 0; st_id < sys_size; st_id++)
		{
			phi[st_id].real = phi_original[st_id].real + curr_lpn_eps * phi_var[st_id].real;
			phi[st_id].imag = phi_original[st_id].imag + curr_lpn_eps * phi_var[st_id].imag;
		}

		double norm_mod_2 = norm_square(phi, sys_size);
		for (int st_id = 0; st_id < sys_size; st_id++)
		{
			phi[st_id].real /= sqrt(norm_mod_2);
			phi[st_id].imag /= sqrt(norm_mod_2);
		}

		for (int st_id = 0; st_id < sys_size; st_id++)
		{
			phi[st_id].real *= sqrt(norm_2);
			phi[st_id].imag *= sqrt(norm_2);
		}

		cb->calc_chars_lpn(ad, 1, 0);

		delta_s = cb->calc_delta_s(ad, 1, 0);

		if (rp->is_pp == 1)
		{
			cout << "lpn: " << 1 << endl;
			cout << "curr_lpn_eps: " << curr_lpn_eps << endl;
			cout << "lpn_deep: " << curr_eps_deep << endl;
			cout << "delta_s: " << delta_s << endl;
			cout << endl;
		}

		if (delta_s > lpn_delta_s_high)
		{
			curr_eps_high -= 0.5 * curr_eps_diff;
			curr_eps_diff = curr_eps_high - curr_eps_low;
			curr_eps_deep++;
			curr_lpn_eps = pow(10.0, (curr_eps_low + curr_eps_high) * 0.5);
			var_first(phi_var, phi_var_double, phi_var_all, scalar_mults_all, ad, cb);
		}

		if (delta_s < lpn_delta_s_low)
		{
			curr_eps_low += 0.5 * curr_eps_diff;
			curr_eps_diff = curr_eps_high - curr_eps_low;
			curr_eps_deep++;
			curr_lpn_eps = pow(10.0, (curr_eps_low + curr_eps_high) * 0.5);
			var_first(phi_var, phi_var_double, phi_var_all, scalar_mults_all, ad, cb);
		}
	}

	// ==== All lpns ====
	for (int tr_id = 2; tr_id < num_trajectories; tr_id++)
	{
		curr_eps_low = lpn_eps_low;
		curr_eps_high = lpn_eps_high;
		curr_eps_diff = curr_eps_high - curr_eps_low;
		curr_eps_deep = 0;
		lpn_eps_start = pow(10.0, (curr_eps_low + curr_eps_high) * 0.5);
		curr_lpn_eps = lpn_eps_start;

		int lpn_id = tr_id - 1;

		phi = &(ed->phi_all[tr_id * sys_size]);

		for (int st_id = 0; st_id < sys_size; st_id++)
		{
			phi_copy[st_id].real = phi[st_id].real;
			phi_copy[st_id].imag = phi[st_id].imag;
		}

		double delta_s = lpn_delta_s_high + 1.0;

		var_not_first(tr_id, phi_var, phi_var_double, phi_var_all, scalar_mults_all, ad, cb);

		while ((delta_s > lpn_delta_s_high || delta_s < lpn_delta_s_low) && (curr_eps_deep < lpn_eps_deep))
		{
			for (int st_id = 0; st_id < sys_size; st_id++)
			{
				phi[st_id].real = phi_copy[st_id].real;
				phi[st_id].imag = phi_copy[st_id].imag;
			}

			double norm_2 = norm_square(phi_original, sys_size);
			for (int st_id = 0; st_id < sys_size; st_id++)
			{
				phi[st_id].real = phi_original[st_id].real + curr_lpn_eps * phi_var[st_id].real;
				phi[st_id].imag = phi_original[st_id].imag + curr_lpn_eps * phi_var[st_id].imag;
			}

			double norm_mod_2 = norm_square(phi, sys_size);
			for (int st_id = 0; st_id < sys_size; st_id++)
			{
				phi[st_id].real /= sqrt(norm_mod_2);
				phi[st_id].imag /= sqrt(norm_mod_2);
			}

			for (int st_id = 0; st_id < sys_size; st_id++)
			{
				phi[st_id].real *= sqrt(norm_2);
				phi[st_id].imag *= sqrt(norm_2);
			}

			cb->calc_chars_lpn(ad, tr_id, 0);

			delta_s = cb->calc_delta_s(ad, tr_id, 0);

			if (rp->is_pp == 1)
			{
				cout << "lpn: " << tr_id << endl;
				cout << "curr_lpn_eps: " << curr_lpn_eps << endl;
				cout << "lpn_deep: " << curr_eps_deep << endl;
				cout << "delta_s: " << delta_s << endl;
				cout << endl;
			}

			if (delta_s > lpn_delta_s_high)
			{
				curr_eps_high -= 0.5 * curr_eps_diff;
				curr_eps_diff = curr_eps_high - curr_eps_low;
				curr_eps_deep++;
				curr_lpn_eps = pow(10.0, (curr_eps_low + curr_eps_high) * 0.5);
				var_not_first(tr_id, phi_var, phi_var_double, phi_var_all, scalar_mults_all, ad, cb);
			}

			if (delta_s < lpn_delta_s_low)
			{
				curr_eps_low += 0.5 * curr_eps_diff;
				curr_eps_diff = curr_eps_high - curr_eps_low;
				curr_eps_deep++;
				curr_lpn_eps = pow(10.0, (curr_eps_low + curr_eps_high) * 0.5);
				var_not_first(tr_id, phi_var, phi_var_double, phi_var_all, scalar_mults_all, ad, cb);
			}
		}
	}

	delete[] phi_var_double;
	delete[] phi_var;
	delete[] phi_var_all;
	delete[] phi_copy;

	delete[] scalar_mults_all;
}

void only_orth(AllData * ad, CoreBehavior *cb, MKL_Complex16 * phi_var_all)
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;
	int num_trajectories = cp->num_trajectories;

	int num_lpns = num_trajectories - 1;

	MKL_Complex16 * phi_var = new MKL_Complex16[sys_size];
	double * phi_var_double = new double[2 * sys_size];
	MKL_Complex16 * scalar_mults_all = new MKL_Complex16[num_lpns];

	// ===== First lpn =====
	var_first_with_history(phi_var, phi_var_double, phi_var_all, scalar_mults_all, ad, cb);


	// ==== All lpns ====
	for (int tr_id = 2; tr_id < num_trajectories; tr_id++)
	{
		var_not_first_with_history(tr_id, phi_var, phi_var_double, phi_var_all, scalar_mults_all, ad, cb);
	}

	delete[] phi_var_double;
	delete[] phi_var;
	delete[] scalar_mults_all;
}

void gs_orth_evo(AllData * ad, CoreBehavior *cb, MKL_Complex16 *phi_var_all)
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;
	int num_trajectories = cp->num_trajectories;

	int num_lpns = num_trajectories - 1;

	double lpn_eps_low = double(cp->params.find("lpn_eps_low")->second);
	double lpn_eps_high = double(cp->params.find("lpn_eps_high")->second);
	int lpn_eps_deep = int(cp->params.find("lpn_eps_deep")->second);

	double curr_eps_low = lpn_eps_low;
	double curr_eps_high = lpn_eps_high;
	double curr_eps_diff = lpn_eps_high - lpn_eps_low;
	int curr_eps_deep = 0;
	double lpn_eps_start = pow(10.0, (curr_eps_low + curr_eps_high) * 0.5);
	double curr_lpn_eps = lpn_eps_start;

	double lpn_delta_s_high = double(cp->params.find("lpn_delta_s_high")->second);
	double lpn_delta_s_low = double(cp->params.find("lpn_delta_s_low")->second);

	MKL_Complex16 * phi_copy = new MKL_Complex16[sys_size];
	MKL_Complex16 * phi_var = new MKL_Complex16[sys_size];
	VSLStreamStatePtr * streams_var = ed->streams_var;
	MKL_Complex16 * phi_original = &(ed->phi_all[0]);
	MKL_Complex16 * phi = NULL;

	// ===== First lpn =====
	phi = &(ed->phi_all[1 * sys_size]);

	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		phi_copy[st_id].real = phi[st_id].real;
		phi_copy[st_id].imag = phi[st_id].imag;
	}

	double delta_s = lpn_delta_s_high + 1.0;

	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		int index = 0 * sys_size + st_id;
		phi_var[st_id].real = phi_var_all[index].real;
		phi_var[st_id].imag = phi_var_all[index].imag;
	}

	while ((delta_s > lpn_delta_s_high || delta_s < lpn_delta_s_low) && (curr_eps_deep < lpn_eps_deep))
	{
		for (int st_id = 0; st_id < sys_size; st_id++)
		{
			phi[st_id].real = phi_copy[st_id].real;
			phi[st_id].imag = phi_copy[st_id].imag;
		}

		double norm_2 = norm_square(phi_original, sys_size);
		for (int st_id = 0; st_id < sys_size; st_id++)
		{
			phi[st_id].real = phi_original[st_id].real + curr_lpn_eps * phi_var[st_id].real;
			phi[st_id].imag = phi_original[st_id].imag + curr_lpn_eps * phi_var[st_id].imag;
		}

		double norm_mod_2 = norm_square(phi, sys_size);
		for (int st_id = 0; st_id < sys_size; st_id++)
		{
			phi[st_id].real /= sqrt(norm_mod_2);
			phi[st_id].imag /= sqrt(norm_mod_2);
		}

		for (int st_id = 0; st_id < sys_size; st_id++)
		{
			phi[st_id].real *= sqrt(norm_2);
			phi[st_id].imag *= sqrt(norm_2);
		}

		cb->calc_chars_lpn(ad, 1, 0);

		delta_s = cb->calc_delta_s(ad, 1, 0);

		if (rp->is_pp == 1)
		{
			cout << "lpn: " << 1 << endl;
			cout << "curr_lpn_eps: " << curr_lpn_eps << endl;
			cout << "lpn_deep: " << curr_eps_deep << endl;
			cout << "delta_s: " << delta_s << endl;
			cout << endl;
		}

		if (delta_s > lpn_delta_s_high)
		{
			curr_eps_high -= 0.5 * curr_eps_diff;
			curr_eps_diff = curr_eps_high - curr_eps_low;
			curr_eps_deep++;
			curr_lpn_eps = pow(10.0, (curr_eps_low + curr_eps_high) * 0.5);
		}

		if (delta_s < lpn_delta_s_low)
		{
			curr_eps_low += 0.5 * curr_eps_diff;
			curr_eps_diff = curr_eps_high - curr_eps_low;
			curr_eps_deep++;
			curr_lpn_eps = pow(10.0, (curr_eps_low + curr_eps_high) * 0.5);
		}
	}

	// ==== All lpns ====
	for (int tr_id = 2; tr_id < num_trajectories; tr_id++)
	{
		curr_eps_low = lpn_eps_low;
		curr_eps_high = lpn_eps_high;
		curr_eps_diff = curr_eps_high - curr_eps_low;
		curr_eps_deep = 0;
		lpn_eps_start = pow(10.0, (curr_eps_low + curr_eps_high) * 0.5);
		curr_lpn_eps = lpn_eps_start;

		int lpn_id = tr_id - 1;

		phi = &(ed->phi_all[tr_id * sys_size]);

		for (int st_id = 0; st_id < sys_size; st_id++)
		{
			phi_copy[st_id].real = phi[st_id].real;
			phi_copy[st_id].imag = phi[st_id].imag;
		}

		double delta_s = lpn_delta_s_high + 1.0;

		for (int st_id = 0; st_id < sys_size; st_id++)
		{
			int index = lpn_id * sys_size + st_id;
			phi_var[st_id].real = phi_var_all[index].real;
			phi_var[st_id].imag = phi_var_all[index].imag;
		}

		while ((delta_s > lpn_delta_s_high || delta_s < lpn_delta_s_low) && (curr_eps_deep < lpn_eps_deep))
		{
			for (int st_id = 0; st_id < sys_size; st_id++)
			{
				phi[st_id].real = phi_copy[st_id].real;
				phi[st_id].imag = phi_copy[st_id].imag;
			}

			double norm_2 = norm_square(phi_original, sys_size);
			for (int st_id = 0; st_id < sys_size; st_id++)
			{
				phi[st_id].real = phi_original[st_id].real + curr_lpn_eps * phi_var[st_id].real;
				phi[st_id].imag = phi_original[st_id].imag + curr_lpn_eps * phi_var[st_id].imag;
			}

			double norm_mod_2 = norm_square(phi, sys_size);
			for (int st_id = 0; st_id < sys_size; st_id++)
			{
				phi[st_id].real /= sqrt(norm_mod_2);
				phi[st_id].imag /= sqrt(norm_mod_2);
			}

			for (int st_id = 0; st_id < sys_size; st_id++)
			{
				phi[st_id].real *= sqrt(norm_2);
				phi[st_id].imag *= sqrt(norm_2);
			}

			cb->calc_chars_lpn(ad, tr_id, 0);

			delta_s = cb->calc_delta_s(ad, tr_id, 0);

			if (rp->is_pp == 1)
			{
				cout << "lpn: " << tr_id << endl;
				cout << "curr_lpn_eps: " << curr_lpn_eps << endl;
				cout << "lpn_deep: " << curr_eps_deep << endl;
				cout << "delta_s: " << delta_s << endl;
				cout << endl;
			}

			if (delta_s > lpn_delta_s_high)
			{
				curr_eps_high -= 0.5 * curr_eps_diff;
				curr_eps_diff = curr_eps_high - curr_eps_low;
				curr_eps_deep++;
				curr_lpn_eps = pow(10.0, (curr_eps_low + curr_eps_high) * 0.5);
			}

			if (delta_s < lpn_delta_s_low)
			{
				curr_eps_low += 0.5 * curr_eps_diff;
				curr_eps_diff = curr_eps_high - curr_eps_low;
				curr_eps_deep++;
				curr_lpn_eps = pow(10.0, (curr_eps_low + curr_eps_high) * 0.5);
			}
		}
	}

	delete[] phi_var;
	delete[] phi_copy;
}

void lambda_lpn(AllData * ad, CoreBehavior *cb, int tr_id, int base_tr_id)
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;

	double lpn_delta_f_high = double(cp->params.find("lpn_delta_f_high")->second);
	double lpn_delta_f_low = double(cp->params.find("lpn_delta_f_low")->second);

	double delta_f = cb->calc_delta_f(ad, tr_id, base_tr_id);

	if ((delta_f > max(lpn_delta_f_high, ed->delta_s[tr_id] * 10.0)) || (delta_f < lpn_delta_f_low))
	{
		ed->lambda[tr_id] += log(delta_f / ed->delta_s[tr_id] + 1.0e-16);
		ed->num_renorms[tr_id] += 1;

		int save_lambdas = int(cp->params.find("save_lambdas")->second);
		if (save_lambdas > 0)
		{
			ed->lambdas[tr_id].push_back(delta_f / ed->delta_s[tr_id]);
			ed->deltas_s[tr_id].push_back(ed->delta_s[tr_id]);
			ed->deltas_f[tr_id].push_back(delta_f);
		}

		copy_stream_lpn(ad, tr_id, base_tr_id);
		var_trajectory_lpn(ad, cb, tr_id, base_tr_id);

		cb->calc_chars_lpn(ad, tr_id, base_tr_id);

		ed->delta_s[tr_id] = cb->calc_delta_s(ad, tr_id, base_tr_id);

		ed->lambda_now[tr_id] = ed->lambda[tr_id] / ed->curr_time;
	}
	else
	{
		cb->calc_chars_lpn(ad, tr_id, base_tr_id);
		ed->lambda_now[tr_id] = (ed->lambda[tr_id] + log(delta_f / ed->delta_s[tr_id] + 1.0e-16)) / ed->curr_time;
	}
}

void lambda_lpn_now(AllData* ad, CoreBehavior* cb, int tr_id, int base_tr_id)
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;
	
	double delta_f = cb->calc_delta_f(ad, tr_id, base_tr_id);

	ed->lambda[tr_id] += log(delta_f / ed->delta_s[tr_id] + 1.0e-16);
	ed->num_renorms[tr_id] += 1;

	int save_lambdas = int(cp->params.find("save_lambdas")->second);
	if (save_lambdas > 0)
	{
		ed->lambdas[tr_id].push_back(delta_f / ed->delta_s[tr_id]);
		ed->deltas_s[tr_id].push_back(ed->delta_s[tr_id]);
		ed->deltas_f[tr_id].push_back(delta_f);
	}

	copy_stream_lpn(ad, tr_id, base_tr_id);
	var_trajectory_lpn(ad, cb, tr_id, base_tr_id);

	cb->calc_chars_lpn(ad, tr_id, base_tr_id);

	ed->delta_s[tr_id] = cb->calc_delta_s(ad, tr_id, base_tr_id);

	ed->lambda_now[tr_id] = ed->lambda[tr_id] / ed->curr_time;
}

void lambda_lpn_per_periods(AllData * ad, CoreBehavior *cb, int tr_id, int base_tr_id, int num_steps_T, int curr_step, int num_periods)
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;

	double delta_f = cb->calc_delta_f(ad, tr_id, base_tr_id);

	if ((curr_step >= (num_steps_T * num_periods)) && (curr_step % (num_steps_T * num_periods) == 0))
	{
		ed->lambda[tr_id] += log(delta_f / ed->delta_s[tr_id] + 1.0e-16);
		ed->num_renorms[tr_id] += 1;

		int save_lambdas = int(cp->params.find("save_lambdas")->second);
		if (save_lambdas > 0)
		{
			ed->lambdas[tr_id].push_back(delta_f / ed->delta_s[tr_id]);
			ed->deltas_s[tr_id].push_back(ed->delta_s[tr_id]);
			ed->deltas_f[tr_id].push_back(delta_f);
		}

		copy_stream_lpn(ad, tr_id, base_tr_id);
		var_trajectory_lpn(ad, cb, tr_id, base_tr_id);

		cb->calc_chars_lpn(ad, tr_id, base_tr_id);

		ed->delta_s[tr_id] = cb->calc_delta_s(ad, tr_id, base_tr_id);

		ed->lambda_now[tr_id] = ed->lambda[tr_id] / ed->curr_time;
	}
	else
	{
		cb->calc_chars_lpn(ad, tr_id, base_tr_id);
		ed->lambda_now[tr_id] = (ed->lambda[tr_id] + log(delta_f / ed->delta_s[tr_id] + 1.0e-16)) / ed->curr_time;
	}
}

void lambda_lpn_all(AllData * ad, CoreBehavior *cb)
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;
	int num_trajectories = cp->num_trajectories;

	double lpn_delta_f_high = double(cp->params.find("lpn_delta_f_high")->second);
	double lpn_delta_f_low = double(cp->params.find("lpn_delta_f_low")->second);

	int num_lpns = num_trajectories - 1;
	MKL_Complex16 * phi_var_all = new MKL_Complex16[num_lpns * sys_size];

	// Firstly we must orth
	only_orth(ad, cb, phi_var_all);

	// Secondly chech for reset
	bool is_reset = false;
	for (int tr_id = 1; tr_id < num_trajectories; tr_id++)
	{
		double delta_f = cb->calc_delta_f(ad, tr_id, 0);
		if ((delta_f > max(lpn_delta_f_high, ed->delta_s[tr_id] * 10.0)) || (delta_f < lpn_delta_f_low))
		{
			is_reset = true;
			break;
		}
	}
	
	// Reset (is nesessary) and calc chars
	for (int tr_id = 1; tr_id < num_trajectories; tr_id++)
	{
		double delta_f = cb->calc_delta_f(ad, tr_id, 0);

		ed->lambda_now[tr_id] = (ed->lambda[tr_id] + log(delta_f / ed->delta_s[tr_id] + 1.0e-16)) / ed->curr_time;

		if (is_reset)
		{
			ed->lambda[tr_id] += log(delta_f / ed->delta_s[tr_id] + 1.0e-16);
		}

		copy_stream_lpn(ad, tr_id, 0);
	}

	if (is_reset)
	{
		gs_orth_evo(ad, cb, phi_var_all);
	}

	for (int tr_id = 1; tr_id < num_trajectories; tr_id++)
	{
		if (is_reset)
		{
			ed->delta_s[tr_id] = cb->calc_delta_s(ad, tr_id, 0);
		}
		cb->calc_chars_lpn(ad, tr_id, 0);
	}

	delete[] phi_var_all;
}

void trans_process_single(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb, int tr_id, int thread_id)
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;
	int num_tp_periods = cp->num_tp_periods;

	MKL_Complex16 * phi = &(ed->phi_all[tr_id * sys_size]);
	double * eta = &(ed->etas_all[tr_id]);

	*eta = 0.0;
	while (*eta == 0.0)
	{
		vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, ed->streams[tr_id], 1, eta, 0.0, 1.0);
	}

	for (int period_id = 0; period_id < num_tp_periods; period_id++)
	{
		pb->one_period(ad, cb, tr_id, thread_id, period_id);
	}
}

void trans_process_single_deep(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb, int tr_id, int thread_id)
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;
	int num_tp_periods = cp->num_tp_periods;

	MKL_Complex16 * phi = &(ed->phi_all[tr_id * sys_size]);
	double * eta = &(ed->etas_all[tr_id]);

	*eta = 0.0;
	while (*eta == 0.0)
	{
		vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, ed->streams[tr_id], 1, eta, 0.0, 1.0);
	}

	for (int period_id = 0; period_id < num_tp_periods; period_id++)
	{
		pb->one_period_trp_deep(ad, cb, tr_id, thread_id, period_id);
	}
}

