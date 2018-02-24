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

void QJPropagateBehavior::init_prop_data_cd(AllData * ad) const
{
	init_splits_deep(ad);
}

void QJPropagateBehavior::free_prop_data_cd(AllData * ad) const
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

void QJPropagateBehavior::one_period_cd_tp(AllData * ad, int tr_id, int thread_id, int period_id) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;

	int num_branches = md->num_ham_qj;
	int num_sub_steps = int(cp->params.find("cd_num_sub_steps")->second);

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		for (int sub_step_id = 0; sub_step_id < num_sub_steps; sub_step_id++)
		{
			one_sub_period_cd(ad, tr_id, part_id, thread_id);
		}
	}
}

void QJPropagateBehavior::one_period_cd_obs(AllData * ad, int tr_id, int thread_id, int period_id) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int cd_dump_deep = int(cp->params.find("cd_dump_deep")->second);
	int is_evo_dump_sep = int(cp->params.find("is_evo_dump_sep")->second);

	int num_branches = md->num_ham_qj;
	int num_sub_steps_per_part = int(cp->params.find("cd_num_sub_steps")->second);
	int num_sub_steps = num_branches * int(cp->params.find("cd_num_sub_steps")->second);

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

			one_sub_period_cd(ad, tr_id, part_id, thread_id);
			calc_chars_std(ad, tr_id);

			if (cd_dump_deep == 1)
			{
				int dump_id = global_point_id + 1;

				evo_chars_std(ad, tr_id, dump_id);

				if (is_evo_dump_sep == 1)
				{
					dump_adr_single(ad, tr_id, true);
				}
			}
		}
	}
}

void QJPropagateBehavior::one_period_sigma_obs(AllData * ad, int tr_id, int thread_id, int period_id) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;

	int cd_dump_deep = int(cp->params.find("cd_dump_deep")->second);
	int is_evo_dump_sep = int(cp->params.find("is_evo_dump_sep")->second);

	int num_branches = md->num_ham_qj;
	int num_sub_steps_per_part = int(cp->params.find("cd_num_sub_steps")->second);
	int num_sub_steps = num_branches * int(cp->params.find("cd_num_sub_steps")->second);

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

			one_sub_period_cd(ad, tr_id, part_id, thread_id);
			calc_chars_std(ad, tr_id);

			if (cd_dump_deep == 1)
			{
				int dump_id = global_point_id + 1;

				evo_chars_std(ad, tr_id, dump_id);

				if (is_evo_dump_sep == 1)
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

void RKPropagateBehavior::init_prop_data_cd(AllData * ad) const
{
	init_rk_cd(ad);
}

void RKPropagateBehavior::free_prop_data_cd(AllData * ad) const
{
	free_rk_cd(ad);
}

void RKPropagateBehavior::one_period(AllData * ad, int tr_id, int thread_id, int period_id) const
{
	rk_period(ad, tr_id, thread_id, period_id);
}

void RKPropagateBehavior::one_period_cd_tp(AllData * ad, int tr_id, int thread_id, int period_id) const
{
	rk_period(ad, tr_id, thread_id, period_id);
}

void RKPropagateBehavior::one_period_cd_obs(AllData * ad, int tr_id, int thread_id, int period_id) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int cd_dump_deep = int(cp->params.find("cd_dump_deep")->second);
	int is_evo_dump_sep = int(cp->params.find("is_evo_dump_sep")->second);

	int num_branches = md->num_ham_qj;
	int num_sub_steps_per_part = int(cp->params.find("cd_num_sub_steps")->second);
	int num_sub_steps = num_branches * int(cp->params.find("cd_num_sub_steps")->second);

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
			rk_period_cd(ad, tr_id, thread_id, time);

			calc_chars_std(ad, tr_id);

			if (cd_dump_deep == 1)
			{
				int dump_id = global_point_id + 1;

				evo_chars_std(ad, tr_id, dump_id);

				if (is_evo_dump_sep == 1)
				{
					dump_adr_single(ad, tr_id, true);
				}
			}
		}
	}
}

void RKPropagateBehavior::one_period_sigma_obs(AllData * ad, int tr_id, int thread_id, int period_id) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;

	int cd_dump_deep = int(cp->params.find("cd_dump_deep")->second);
	int is_evo_dump_sep = int(cp->params.find("is_evo_dump_sep")->second);

	int num_branches = md->num_ham_qj;
	int num_sub_steps_per_part = int(cp->params.find("cd_num_sub_steps")->second);
	int num_sub_steps = num_branches * int(cp->params.find("cd_num_sub_steps")->second);

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
			rk_period_cd(ad, tr_id, thread_id, time);

			calc_chars_std(ad, tr_id);

			if (cd_dump_deep == 1)
			{
				int dump_id = global_point_id + 1;

				evo_chars_std(ad, tr_id, dump_id);

				if (is_evo_dump_sep == 1)
				{
					dump_adr_single(ad, tr_id, true);
				}
			}
		}
	}
}
