#include "propagator.h"
#include "experiment.h"
#include "qj_proc.h"
#include "rk_proc.h"

void QJPropagateBehavior::init_prop_data(AllData * ad) const
{
	init_splits(ad);
}

void QJPropagateBehavior::free_prop_data(AllData * ad) const
{
	free_splits(ad);
}

void QJPropagateBehavior::init_prop_data_deep(AllData * ad) const
{
	init_splits_deep(ad);
}

void QJPropagateBehavior::free_prop_data_deep(AllData * ad) const
{
	free_splits_deep(ad);
}

void QJPropagateBehavior::one_period(AllData * ad, int tr_id, int thread_id, int period_id) const
{
	MainData * md = ad->md;

	int sys_size = md->sys_size;
	Split * head = &(md->splits)[thread_id];

	for (unsigned int b_id = 0; b_id < head->counter; b_id++)
	{
		one_period_branch(ad, head, tr_id, &(head->next)[b_id]);
	}
}

void QJPropagateBehavior::one_period_tp_deep(AllData * ad, int tr_id, int thread_id, int period_id) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;

	int num_branches = md->num_ham_qj;
	int num_sub_steps = int(cp->params.find("deep_num_steps")->second);

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		for (int sub_step_id = 0; sub_step_id < num_sub_steps; sub_step_id++)
		{
			one_sub_period_deep(ad, tr_id, part_id, thread_id);
		}
	}
}

void QJPropagateBehavior::one_period_obs_deep(AllData * ad, int tr_id, int thread_id, int period_id) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);

	int num_branches = md->num_ham_qj;
	int num_sub_steps_per_part = int(cp->params.find("deep_num_steps")->second);
	int num_sub_steps = num_branches * int(cp->params.find("deep_num_steps")->second);

	int dump_point_id = 0;
	int global_point_id = 0;

	int step_id = 0;

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		for (int sub_step_id = 0; sub_step_id < num_sub_steps_per_part; sub_step_id++)
		{
			step_id = part_id * num_sub_steps_per_part + sub_step_id;

			global_point_id = period_id * num_sub_steps + dump_point_id;
			dump_point_id++;

			one_sub_period_deep(ad, tr_id, part_id, thread_id);
			calc_chars_std(ad, tr_id);

			int dump_id = global_point_id + 1;

			evo_chars_std(ad, tr_id, dump_id);

			if (dump_evo_sep == 1)
			{
				dump_adr_single(ad, tr_id, true);
			}
		}
	}
}

void QJPropagateBehavior::one_period_obs_deep_lpn(AllData * ad, int tr_id, int thread_id, int period_id) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);

	int num_branches = md->num_ham_qj;
	int num_sub_steps_per_part = int(cp->params.find("deep_num_steps")->second);
	int num_sub_steps = num_branches * int(cp->params.find("deep_num_steps")->second);

	int dump_point_id = 0;
	int global_point_id = 0;

	int step_id = 0;

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		for (int sub_step_id = 0; sub_step_id < num_sub_steps_per_part; sub_step_id++)
		{
			step_id = part_id * num_sub_steps_per_part + sub_step_id;

			global_point_id = period_id * num_sub_steps + dump_point_id;
			dump_point_id++;

			one_sub_period_deep(ad, tr_id, part_id, thread_id);
			calc_chars_std(ad, tr_id);
			calc_chars_lpn(ad, tr_id);
			ed->energy[tr_id] = get_energy(ad, tr_id);

			int dump_id = global_point_id + 1;

			evo_chars_std(ad, tr_id, dump_id);
			evo_chars_lpn(ad, tr_id, dump_id);

			if (dump_evo_sep == 1)
			{
				dump_adr_single(ad, tr_id, true);
			}
		}
	}
}

void QJPropagateBehavior::one_period_obs_cd(AllData * ad, int tr_id, int thread_id, int period_id) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int deep_dump = int(cp->params.find("deep_dump")->second);
	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);

	int num_branches = md->num_ham_qj;
	int num_sub_steps_per_part = int(cp->params.find("deep_num_steps")->second);
	int num_sub_steps = num_branches * int(cp->params.find("deep_num_steps")->second);

	int dump_point_id = 0;
	int curr_point_id = 0;
	int global_point_id = 0;

	int step_id = 0;

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		for (int sub_step_id = 0; sub_step_id < num_sub_steps_per_part; sub_step_id++)
		{
			step_id = part_id * num_sub_steps_per_part + sub_step_id;

			global_point_id = period_id * num_sub_steps + dump_point_id;

			dump_point_id++;

			for (int cd_st_id = 0; cd_st_id < ed->cd_dim; cd_st_id++)
			{
				curr_point_id = global_point_id - cd_st_id;

				if (curr_point_id >= 0 && curr_point_id < ed->cd_num_points)
				{
					ed->cd_rec_data[tr_id][curr_point_id][cd_st_id] = ed->mean[tr_id];
				}
			}

			one_sub_period_deep(ad, tr_id, part_id, thread_id);
			calc_chars_std(ad, tr_id);

			if (deep_dump == 1)
			{
				int dump_id = global_point_id + 1;

				evo_chars_std(ad, tr_id, dump_id);

				if (dump_evo_sep == 1)
				{
					dump_adr_single(ad, tr_id, true);
				}
			}
		}
	}
}

void QJPropagateBehavior::one_period_obs_sigma(AllData * ad, int tr_id, int thread_id, int period_id) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;

	int deep_dump = int(cp->params.find("deep_dump")->second);
	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);

	int num_branches = md->num_ham_qj;
	int num_sub_steps_per_part = int(cp->params.find("deep_num_steps")->second);
	int num_sub_steps = num_branches * int(cp->params.find("deep_num_steps")->second);

	int dump_point_id = 0;
	int curr_point_id = 0;
	int global_point_id = 0;

	int step_id = 0;

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		for (int sub_step_id = 0; sub_step_id < num_sub_steps_per_part; sub_step_id++)
		{
			global_point_id = period_id * num_sub_steps + dump_point_id;

			dump_point_id++;

			one_sub_period_deep(ad, tr_id, part_id, thread_id);
			calc_chars_std(ad, tr_id);

			if (deep_dump == 1)
			{
				int dump_id = global_point_id + 1;

				evo_chars_std(ad, tr_id, dump_id);

				if (dump_evo_sep == 1)
				{
					dump_adr_single(ad, tr_id, true);
				}
			}
		}
	}
}

void RKPropagateBehavior::init_prop_data(AllData * ad) const
{
	init_rk(ad);
}

void RKPropagateBehavior::free_prop_data(AllData * ad) const
{
	free_rk(ad);
}

void RKPropagateBehavior::init_prop_data_deep(AllData * ad) const
{
	init_rk_cd(ad);
}

void RKPropagateBehavior::free_prop_data_deep(AllData * ad) const
{
	free_rk_cd(ad);
}

void RKPropagateBehavior::one_period(AllData * ad, int tr_id, int thread_id, int period_id) const
{
	rk_period(ad, tr_id, thread_id, period_id);
}

void RKPropagateBehavior::one_period_tp_deep(AllData * ad, int tr_id, int thread_id, int period_id) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int num_branches = md->num_ham_qj;
	int num_sub_steps_per_part = int(cp->params.find("deep_num_steps")->second);
	int num_sub_steps = num_branches * int(cp->params.find("deep_num_steps")->second);

	int step_id = 0;

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		for (int sub_step_id = 0; sub_step_id < num_sub_steps_per_part; sub_step_id++)
		{
			step_id = part_id * num_sub_steps_per_part + sub_step_id;

			double time = double(period_id) * md->T + double(step_id) * md->T / double(num_sub_steps);
			rk_period_deep(ad, tr_id, thread_id, time);
		}
	}
}

void RKPropagateBehavior::one_period_obs_deep(AllData * ad, int tr_id, int thread_id, int period_id) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);

	int num_branches = md->num_ham_qj;
	int num_sub_steps_per_part = int(cp->params.find("deep_num_steps")->second);
	int num_sub_steps = num_branches * int(cp->params.find("deep_num_steps")->second);

	int dump_point_id = 0;
	int global_point_id = 0;

	int step_id = 0;

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		for (int sub_step_id = 0; sub_step_id < num_sub_steps_per_part; sub_step_id++)
		{
			step_id = part_id * num_sub_steps_per_part + sub_step_id;

			global_point_id = period_id * num_sub_steps + dump_point_id;

			dump_point_id++;

			double time = double(period_id) * md->T + double(step_id) * md->T / double(num_sub_steps);
			rk_period_deep(ad, tr_id, thread_id, time);

			calc_chars_std(ad, tr_id);

			int dump_id = global_point_id + 1;

			evo_chars_std(ad, tr_id, dump_id);

			if (dump_evo_sep == 1)
			{
				dump_adr_single(ad, tr_id, true);
			}
		}
	}
}

void RKPropagateBehavior::one_period_obs_deep_lpn(AllData * ad, int tr_id, int thread_id, int period_id) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);

	int num_branches = md->num_ham_qj;
	int num_sub_steps_per_part = int(cp->params.find("deep_num_steps")->second);
	int num_sub_steps = num_branches * int(cp->params.find("deep_num_steps")->second);

	int dump_point_id = 0;
	int global_point_id = 0;

	int step_id = 0;

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		for (int sub_step_id = 0; sub_step_id < num_sub_steps_per_part; sub_step_id++)
		{
			step_id = part_id * num_sub_steps_per_part + sub_step_id;

			global_point_id = period_id * num_sub_steps + dump_point_id;
			dump_point_id++;

			double time = double(period_id) * md->T + double(step_id) * md->T / double(num_sub_steps);
			rk_period_deep(ad, tr_id, thread_id, time);

			calc_chars_std(ad, tr_id);
			calc_chars_lpn(ad, tr_id);
			ed->energy[tr_id] = get_energy(ad, tr_id);

			int dump_id = global_point_id + 1;

			evo_chars_std(ad, tr_id, dump_id);
			evo_chars_lpn(ad, tr_id, dump_id);

			if (dump_evo_sep == 1)
			{
				dump_adr_single(ad, tr_id, true);
			}
		}
	}
}

void RKPropagateBehavior::one_period_obs_cd(AllData * ad, int tr_id, int thread_id, int period_id) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int deep_dump = int(cp->params.find("deep_dump")->second);
	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);

	int num_branches = md->num_ham_qj;
	int num_sub_steps_per_part = int(cp->params.find("deep_num_steps")->second);
	int num_sub_steps = num_branches * int(cp->params.find("deep_num_steps")->second);

	int dump_point_id = 0;
	int curr_point_id = 0;
	int global_point_id = 0;

	int step_id = 0;

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		for (int sub_step_id = 0; sub_step_id < num_sub_steps_per_part; sub_step_id++)
		{
			step_id = part_id * num_sub_steps_per_part + sub_step_id;

			global_point_id = period_id * num_sub_steps + dump_point_id;

			dump_point_id++;

			for (int cd_st_id = 0; cd_st_id < ed->cd_dim; cd_st_id++)
			{
				curr_point_id = global_point_id - cd_st_id;

				if (curr_point_id >= 0 && curr_point_id < ed->cd_num_points)
				{
					ed->cd_rec_data[tr_id][curr_point_id][cd_st_id] = ed->mean[tr_id];
				}
			}

			double time = double(period_id) * md->T + double(step_id) * md->T / double(num_sub_steps);
			rk_period_deep(ad, tr_id, thread_id, time);

			calc_chars_std(ad, tr_id);

			if (deep_dump == 1)
			{
				int dump_id = global_point_id + 1;

				evo_chars_std(ad, tr_id, dump_id);

				if (dump_evo_sep == 1)
				{
					dump_adr_single(ad, tr_id, true);
				}
			}
		}
	}
}

void RKPropagateBehavior::one_period_obs_sigma(AllData * ad, int tr_id, int thread_id, int period_id) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;

	int deep_dump = int(cp->params.find("deep_dump")->second);
	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);

	int num_branches = md->num_ham_qj;
	int num_sub_steps_per_part = int(cp->params.find("deep_num_steps")->second);
	int num_sub_steps = num_branches * int(cp->params.find("deep_num_steps")->second);

	int dump_point_id = 0;
	int curr_point_id = 0;
	int global_point_id = 0;

	int step_id = 0;

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		for (int sub_step_id = 0; sub_step_id < num_sub_steps_per_part; sub_step_id++)
		{
			step_id = part_id * num_sub_steps_per_part + sub_step_id;
			global_point_id = period_id * num_sub_steps + dump_point_id;

			dump_point_id++;

			double time = double(period_id) * md->T + double(step_id) * md->T / double(num_sub_steps);
			rk_period_deep(ad, tr_id, thread_id, time);

			calc_chars_std(ad, tr_id);

			if (deep_dump == 1)
			{
				int dump_id = global_point_id + 1;

				evo_chars_std(ad, tr_id, dump_id);

				if (dump_evo_sep == 1)
				{
					dump_adr_single(ad, tr_id, true);
				}
			}
		}
	}
}
