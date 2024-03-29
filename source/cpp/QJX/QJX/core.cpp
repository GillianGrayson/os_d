#include "core.h"
#include "qj_proc.h"
#include "rk_proc.h"
#include "experiment.h"

void DimerCoreBehaviour::init_splits(AllData * ad) const
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;

	int num_threads = rp->num_threads;

	md->structure = init_split_structure_dimer(ad);
	md->splits = new Split[num_threads];
	for (int th_id = 0; th_id < num_threads; th_id++)
	{
		copy_struct_not_member(md->structure, &(md->splits)[th_id]);
	}

	cout << "Split initialization complete" << endl;
}

void DimerCoreBehaviour::free_splits(AllData * ad) const
{
	RunParam * rp = ad->rp;
	MainData * md = ad->md;

	int num_threads = rp->num_threads;

	for (int i = 0; i < num_threads; i++)
	{
		delete_split_struct_not_member(&(md->splits[i]));
	}

	delete(md->splits);
	delete_split_struct(md->structure);
}

void DimerCoreBehaviour::init_splits_deep(AllData * ad) const
{
	RunParam * rp = ad->rp;
	MainData * md = ad->md;

	int num_branches = md->num_ham_qj;
	int num_threads = rp->num_threads;

	int num_total = num_threads * num_branches;

	md->structure = init_split_structure_dimer_deep(ad);
	md->splits = new Split[num_total];

	for (int b_id = 0; b_id < num_branches; b_id++)
	{
		for (int th_id = 0; th_id < num_threads; th_id++)
		{
			int index = b_id * num_threads + th_id;
			copy_struct_not_member(&(md->structure)[b_id], &(md->splits)[index]);
		}
	}
}

void DimerCoreBehaviour::free_splits_deep(AllData * ad) const
{
	RunParam * rp = ad->rp;
	MainData * md = ad->md;

	int num_branches = md->num_ham_qj;
	int num_threads = rp->num_threads;

	for (int b_id = 0; b_id < num_branches; b_id++)
	{
		for (int th_id = 0; th_id < num_threads; th_id++)
		{
			int index = b_id * num_threads + th_id;
			delete_split_struct_not_member(&(md->splits[index]));
		}
	}

	delete(md->splits);
	delete_split_struct(md->structure);
}

void DimerCoreBehaviour::ex_period(AllData * ad, int tr_id, int th_id, int period_id) const
{
	MainData * md = ad->md;

	int sys_size = md->sys_size;
	Split * head = &(md->splits)[th_id];

	for (int b_id = 0; b_id < head->counter; b_id++)
	{
		one_period_branch(ad, head, tr_id, &(head->next)[b_id]);
	}
}

void DimerCoreBehaviour::ex_period_trp_deep(AllData * ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;

	int num_branches = md->num_ham_qj;
	int num_sub_steps = int(cp->params.find("deep_num_steps")->second);

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		for (int sub_step_id = 0; sub_step_id < num_sub_steps; sub_step_id++)
		{
			one_sub_period_deep(ad, tr_id, part_id, th_id);
		}
	}
}

void DimerCoreBehaviour::ex_period_obs_deep(AllData * ad, int tr_id, int th_id, int period_id) const
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

			one_sub_period_deep(ad, tr_id, part_id, th_id);
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

void DimerCoreBehaviour::ex_period_obs_deep_lpn(AllData * ad, int period_id) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);

	int num_trajectories = cp->num_trajectories;

	int num_branches = md->num_ham_qj;
	int num_sub_steps_per_part = int(cp->params.find("deep_num_steps")->second);
	int num_sub_steps = num_branches * int(cp->params.find("deep_num_steps")->second);

	int dump_point_id = 0;
	int global_point_id = 0;

	int step_id = 0;

	CoreBehavior * tmp = new DimerCoreBehaviour;

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		for (int sub_step_id = 0; sub_step_id < num_sub_steps_per_part; sub_step_id++)
		{
			step_id = part_id * num_sub_steps_per_part + sub_step_id;

			global_point_id = period_id * num_sub_steps + dump_point_id;
			dump_point_id++;
			int dump_id = global_point_id + 1;

			ed->curr_time = double(dump_id) / double(num_sub_steps) * md->T;

			one_sub_period_deep(ad, 0, part_id, 0);
			calc_chars_std(ad, 0);
			calc_chars_lpn(ad, 0, 0);
			evo_chars_std(ad, 0, dump_id);
			evo_chars_lpn(ad, 0, dump_id);

#pragma omp parallel for
			for (int tr_id = 1; tr_id < num_trajectories; tr_id++)
			{
				int thread_id = omp_get_thread_num();
				one_sub_period_deep(ad, tr_id, part_id, thread_id);
				calc_chars_std(ad, tr_id);
				lambda_lpn(ad, tmp, tr_id, 0);
				evo_chars_std(ad, tr_id, dump_id);
				evo_chars_lpn(ad, tr_id, dump_id);
			}

#pragma omp parallel for
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				if (dump_evo_sep == 1)
				{
					dump_adr_single(ad, tr_id, true);
				}
			}
		}
	}

	delete tmp;
}

void DimerCoreBehaviour::ex_period_obs_deep_mult_lpn(AllData * ad, int period_id) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);

	int num_trajectories = cp->num_trajectories;

	int num_branches = md->num_ham_qj;
	int num_sub_steps_per_part = int(cp->params.find("deep_num_steps")->second);
	int num_sub_steps = num_branches * int(cp->params.find("deep_num_steps")->second);

	int lambda_per_periods = int(cp->params.find("lambda_per_periods")->second);

	int dump_point_id = 0;
	int global_point_id = 0;

	int step_id = 0;

	CoreBehavior * tmp = new DimerCoreBehaviour;

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		for (int sub_step_id = 0; sub_step_id < num_sub_steps_per_part; sub_step_id++)
		{
			step_id = part_id * num_sub_steps_per_part + sub_step_id;

			global_point_id = period_id * num_sub_steps + dump_point_id;
			dump_point_id++;
			int dump_id = global_point_id + 1;

			ed->curr_time = double(dump_id) / double(num_sub_steps) * md->T;

#pragma omp parallel for
			for (int tr_id = 0; tr_id < num_trajectories / 2; tr_id++)
			{
				int thread_id = omp_get_thread_num();
				one_sub_period_deep(ad, tr_id, part_id, thread_id);
				calc_chars_std(ad, tr_id);
				calc_chars_lpn(ad, tr_id, tr_id);
				evo_chars_std(ad, tr_id, dump_id);
				evo_chars_lpn(ad, tr_id, dump_id);
			}

#pragma omp parallel for
			for (int tr_id = num_trajectories / 2; tr_id < num_trajectories; tr_id++)
			{
				int thread_id = omp_get_thread_num();
				one_sub_period_deep(ad, tr_id, part_id, thread_id);
				calc_chars_std(ad, tr_id);
				if (lambda_per_periods == 0)
				{
					lambda_lpn(ad, tmp, tr_id, tr_id - num_trajectories / 2);
				}
				if (lambda_per_periods == 1)
				{
					if ((part_id == num_branches - 1) && (sub_step_id == num_sub_steps_per_part - 1))
					{
						lambda_lpn_now(ad, tmp, tr_id, tr_id - num_trajectories / 2);
					}
				}
				evo_chars_std(ad, tr_id, dump_id);
				evo_chars_lpn(ad, tr_id, dump_id);
			}

#pragma omp parallel for
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				if (dump_evo_sep == 1)
				{
					dump_adr_single(ad, tr_id, true);
				}
			}
		}
	}

	delete tmp;
}

void DimerCoreBehaviour::ex_period_obs_deep_lpn_per_period(struct AllData* ad, int period_id, int num_periods) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);

	int num_trajectories = cp->num_trajectories;

	int num_branches = md->num_ham_qj;
	int num_sub_steps_per_part = int(cp->params.find("deep_num_steps")->second);
	int num_sub_steps = num_branches * int(cp->params.find("deep_num_steps")->second);

	int dump_point_id = 0;
	int global_point_id = 0;

	int step_id = 0;

	CoreBehavior * tmp = new DimerCoreBehaviour;

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		for (int sub_step_id = 0; sub_step_id < num_sub_steps_per_part; sub_step_id++)
		{
			step_id = part_id * num_sub_steps_per_part + sub_step_id;

			global_point_id = period_id * num_sub_steps + dump_point_id;
			dump_point_id++;
			int dump_id = global_point_id + 1;

			ed->curr_time = double(dump_id) / double(num_sub_steps) * md->T;

			one_sub_period_deep(ad, 0, part_id, 0);
			calc_chars_std(ad, 0);
			calc_chars_lpn(ad, 0, 0);
			evo_chars_std(ad, 0, dump_id);
			evo_chars_lpn(ad, 0, dump_id);

#pragma omp parallel for
			for (int tr_id = 1; tr_id < num_trajectories; tr_id++)
			{
				int thread_id = omp_get_thread_num();
				one_sub_period_deep(ad, tr_id, part_id, thread_id);
				calc_chars_std(ad, tr_id);
				lambda_lpn_per_periods(ad, tmp, tr_id, 0, num_sub_steps, global_point_id, num_periods);
				evo_chars_std(ad, tr_id, dump_id);
				evo_chars_lpn(ad, tr_id, dump_id);
			}

#pragma omp parallel for
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				if (dump_evo_sep == 1)
				{
					dump_adr_single(ad, tr_id, true);
				}
			}
		}
	}

	delete tmp;
}

void DimerCoreBehaviour::ex_period_obs_deep_cd(AllData * ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

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

			one_sub_period_deep(ad, tr_id, part_id, th_id);
			calc_chars_std(ad, tr_id);
		}
	}
}

void DimerCoreBehaviour::ex_period_obs_deep_sigma(AllData * ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;

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

			one_sub_period_deep(ad, tr_id, part_id, th_id);
			calc_chars_std(ad, tr_id);
		}
	}
}

void DimerCoreBehaviour::rk_period(AllData * ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	for (int step_id = 0; step_id < cp->rk_ns; step_id++)
	{
		ed->times_all[tr_id] = double(period_id) * md->T + double(step_id) * ed->rk_step;
		rk_step_dimer(ad, tr_id, th_id, ed->rk_step);
	}
}

void DimerCoreBehaviour::rk_period_trp_deep(AllData * ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int num_branches = md->num_ham_qj;
	int num_sub_steps_per_part = int(cp->params.find("deep_num_steps")->second);

	int step_id = 0;

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		for (int sub_step_id = 0; sub_step_id < num_sub_steps_per_part; sub_step_id++)
		{
			step_id = part_id * num_sub_steps_per_part + sub_step_id;

			for (int in_step_id = 0; in_step_id < cp->rk_ns; in_step_id++)
			{
				ed->times_all[tr_id] = double(period_id) * md->T + double(step_id * cp->rk_ns + in_step_id) * ed->rk_step;
				rk_step_dimer(ad, tr_id, th_id, ed->rk_step);
			}
		}
	}
}

void DimerCoreBehaviour::rk_period_obs_deep(AllData * ad, int tr_id, int th_id, int period_id) const
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

			for (int in_step_id = 0; in_step_id < cp->rk_ns; in_step_id++)
			{
				ed->times_all[tr_id] = double(period_id) * md->T + double(step_id * cp->rk_ns + in_step_id) * ed->rk_step;
				rk_step_dimer(ad, tr_id, th_id, ed->rk_step);
			}

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

void DimerCoreBehaviour::rk_period_obs_deep_lpn(AllData * ad, int period_id) const
{
	stringstream msg;
	msg << "Error: need to add functionality" << endl;
	Error(msg.str());

	/*ConfigParam * cp = ad->cp;
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
			for (int in_step_id = 0; in_step_id < cp->rk_ns; in_step_id++)
			{
				ed->times_all[tr_id] = time + double(in_step_id) * ed->rk_step;
				rk_step_dimer(ad, tr_id, th_id, ed->rk_step);
			}

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
	}*/
}

void DimerCoreBehaviour::rk_period_obs_deep_cd(AllData * ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

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

			for (int in_step_id = 0; in_step_id < cp->rk_ns; in_step_id++)
			{
				ed->times_all[tr_id] = double(period_id) * md->T + double(step_id * cp->rk_ns + in_step_id) * ed->rk_step;
				rk_step_dimer(ad, tr_id, th_id, ed->rk_step);
			}

			calc_chars_std(ad, tr_id);
		}
	}
}

void DimerCoreBehaviour::rk_period_obs_deep_sigma(AllData * ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

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

			dump_point_id++;

			for (int in_step_id = 0; in_step_id < cp->rk_ns; in_step_id++)
			{
				ed->times_all[tr_id] = double(period_id) * md->T + double(step_id * cp->rk_ns + in_step_id) * ed->rk_step;
				rk_step_dimer(ad, tr_id, th_id, ed->rk_step);
			}

			calc_chars_std(ad, tr_id);
		}
	}
}

void DimerCoreBehaviour::calc_chars_std_start(AllData * ad, int tr_id) const
{
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;
	MKL_Complex16 * phi = &(ed->phi_all[tr_id * sys_size]);
	double * adr = &(ed->abs_diag_rho_all[tr_id * sys_size]);

	double norm = norm_square(phi, sys_size);

	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		adr[st_id] = mult_scalar_double(mult_scalar_complex(&phi[st_id], &phi[st_id], 1), 1.0 / norm).real;
	}

	double mean = get_mean_simple(adr, sys_size);
	double dispersion = get_dispersion_simple(mean, mean);
	double m2 = get_m2(adr, sys_size, mean);
	double energy = get_energy(ad, tr_id);

	ed->mean_start[tr_id] = mean;

	ed->norm[tr_id] = norm;
	ed->mean[tr_id] = mean;
	ed->dispersion[tr_id] = dispersion;
	ed->m2[tr_id] = m2;
	ed->energy[tr_id] = energy;

	calc_random_obs(ad, tr_id);
}

void DimerCoreBehaviour::calc_chars_std(AllData * ad, int tr_id) const
{
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;
	MKL_Complex16 * phi = &(ed->phi_all[tr_id * sys_size]);
	double * adr = &(ed->abs_diag_rho_all[tr_id * sys_size]);

	double norm = norm_square(phi, sys_size);

	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		adr[st_id] = mult_scalar_double(mult_scalar_complex(&phi[st_id], &phi[st_id], 1), 1.0 / norm).real;
	}

	double mean = get_mean_simple(adr, sys_size);
	double dispersion = get_dispersion_simple(mean, ed->mean_start[tr_id]);
	double m2 = get_m2(adr, sys_size, mean);
	double energy = get_energy(ad, tr_id);

	ed->norm[tr_id] = norm;
	ed->mean[tr_id] = mean;
	ed->dispersion[tr_id] = dispersion;
	ed->m2[tr_id] = m2;
	ed->energy[tr_id] = energy;

	calc_random_obs(ad, tr_id);
}

void DimerCoreBehaviour::calc_chars_lpn_start(AllData * ad, int tr_id, int base_tr_id) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;

	MKL_Complex16 * phi = &(ed->phi_all[tr_id * sys_size]);
	double * adr = &(ed->abs_diag_rho_all[tr_id * sys_size]);

	double norm_2 = norm_square(phi, sys_size);

	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		adr[st_id] = mult_scalar_double(mult_scalar_complex(&phi[st_id], &phi[st_id], 1), 1.0 / norm_2).real;
	}

	double lambda = 0.0;
	double lambda_now = 0.0;
	double mean_lpn = ed->mean[tr_id];
	double energy_lpn = ed->energy[tr_id];

	double delta_s = this->calc_delta_f(ad, tr_id, base_tr_id); // Important! Here we use calc_delta_f not calc_delta_s

	ed->lambda[tr_id] = lambda;
	ed->lambda_now[tr_id] = lambda_now;
	ed->mean_lpn[tr_id] = mean_lpn;
	ed->energy_lpn[tr_id] = energy_lpn;
	ed->delta_s[tr_id] = delta_s;

	calc_random_obs_lpn_start(ad, tr_id);
}

void DimerCoreBehaviour::calc_chars_lpn(AllData * ad, int tr_id, int base_tr_id) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;

	MKL_Complex16 * phi = &(ed->phi_all[tr_id * sys_size]);
	double * adr = &(ed->abs_diag_rho_all[tr_id * sys_size]);

	double norm_2 = norm_square(phi, sys_size);

	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		adr[st_id] = mult_scalar_double(mult_scalar_complex(&phi[st_id], &phi[st_id], 1), 1.0 / norm_2).real;
	}

	double mean_lpn = get_mean_simple(adr, sys_size);
	double energy_lpn = get_energy(ad, tr_id);

	ed->mean_lpn[tr_id] = mean_lpn;
	ed->energy_lpn[tr_id] = energy_lpn;

	calc_random_obs_lpn(ad, tr_id);
}

double DimerCoreBehaviour::calc_T(AllData * ad) const
{
	MainData * md = ad->md;
	double T = md->T;
	return T;
}

void DimerCoreBehaviour::evo_chars_std(AllData * ad, int tr_id, int dump_id) const
{
	ConfigParam * cp = ad->cp;
	ExpData * ed = ad->ed;

	int num_trajectories = cp->num_trajectories;
	int dump_num_total = ed->dump_num_total;

	int index = tr_id * dump_num_total + dump_id;

	ed->norm_evo[index] = ed->norm[tr_id];
	ed->mean_evo[index] = ed->mean[tr_id];
	ed->dispersion_evo[index] = ed->dispersion[tr_id];
	ed->m2_evo[index] = ed->m2[tr_id];
	ed->energy_evo[index] = ed->energy[tr_id];

	int num_random_obs = int(cp->params.find("num_random_obs")->second);
	for (int obs_id = 0; obs_id < num_random_obs; obs_id++)
	{
		ed->random_obs_evo[tr_id][obs_id].push_back(ed->random_obs[tr_id][obs_id]);
	}
}

void DimerCoreBehaviour::evo_chars_lpn(AllData * ad, int tr_id, int dump_id) const
{
	ConfigParam * cp = ad->cp;
	ExpData * ed = ad->ed;

	int num_trajectories = cp->num_trajectories;
	int dump_num_total = ed->dump_num_total;

	int index = tr_id * dump_num_total + dump_id;

	ed->lambda_evo[index] = ed->lambda_now[tr_id];
	ed->mean_lpn_evo[index] = ed->mean_lpn[tr_id];
	ed->energy_lpn_evo[index] = ed->energy_lpn[tr_id];

	int num_random_obs = int(cp->params.find("num_random_obs")->second);
	for (int obs_id = 0; obs_id < num_random_obs; obs_id++)
	{
		ed->random_obs_lpn_evo[tr_id][obs_id].push_back(ed->random_obs_lpn[tr_id][obs_id]);
	}
}

double DimerCoreBehaviour::calc_delta_s(AllData * ad, int tr_id, int base_tr_id) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;

	int lpn_type = int(cp->params.find("lpn_type")->second);
	int num_random_obs = int(ad->cp->params.find("num_random_obs")->second);

	double delta_s = 0.0;
	if (lpn_type == -1)
	{
		delta_s = fabs(ed->mean_lpn[tr_id] - ed->mean_lpn[base_tr_id]) / double(sys_size);
	}
	else if (lpn_type == -2)
	{
		delta_s = fabs(ed->energy_lpn[tr_id] - ed->energy_lpn[base_tr_id]) / ed->max_energy;
	}
	else if (lpn_type >= 0 && lpn_type < num_random_obs)
	{
		complex<double> base = ed->random_obs_lpn[base_tr_id][lpn_type];
		complex<double> var = ed->random_obs_lpn[tr_id][lpn_type];
		double tmp = pow((base.real() - var.real()), 2) + pow((base.imag() - var.imag()), 2);
		delta_s = sqrt(tmp);
	}
	else
	{
		stringstream msg;
		msg << "Error: Wrong lpn_type: " << lpn_type << endl;
		Error(msg.str());
	}

	return delta_s;
}

double DimerCoreBehaviour::calc_delta_f(AllData * ad, int tr_id, int base_tr_id) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;

	int lpn_type = int(cp->params.find("lpn_type")->second);
	int num_random_obs = int(ad->cp->params.find("num_random_obs")->second);

	double delta_f = 0.0;
	if (lpn_type == -1)
	{
		delta_f = fabs(ed->mean[tr_id] - ed->mean[base_tr_id]) / double(sys_size);
	}
	else if (lpn_type == -2)
	{
		delta_f = fabs(ed->energy[tr_id] - ed->energy[base_tr_id]) / ed->max_energy;
	}
	else if (lpn_type >= 0 && lpn_type < num_random_obs)
	{
		complex<double> base = ed->random_obs[base_tr_id][lpn_type];
		complex<double> var = ed->random_obs[tr_id][lpn_type];
		double tmp = pow((base.real() - var.real()), 2) + pow((base.imag() - var.imag()), 2);
		delta_f = sqrt(tmp);
	}
	else
	{
		stringstream msg;
		msg << "Error: Wrong lpn_type: " << lpn_type << endl;
		Error(msg.str());
	}

	return delta_f;
}

void DimerCoreBehaviour::calc_ci(AllData * ad, int tr_id) const
{
	calc_ci_double(ad, tr_id);
}

void DimerCoreBehaviour::dump_std(AllData * ad) const
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	ExpData * ed = ad->ed;

	int dump_obs = int(cp->params.find("dump_obs")->second);

	if (dump_obs == 1)
	{
		int num_trajectories = cp->num_trajectories;

		double * norm = ed->norm;
		double * mean_start = ed->mean_start;
		double * mean = ed->mean;
		double * dispersion = ed->dispersion;
		double * m2 = ed->m2;
		double * energy = ed->energy;

		string fn;

		fn = rp->path + "norm" + cp->fn_suffix;
		save_double_data(fn, norm, num_trajectories, 16, false);

		fn = rp->path + "mean_start" + cp->fn_suffix;
		save_double_data(fn, mean_start, num_trajectories, 16, false);

		fn = rp->path + "mean" + cp->fn_suffix;
		save_double_data(fn, mean, num_trajectories, 16, false);

		fn = rp->path + "dispersion" + cp->fn_suffix;
		save_double_data(fn, dispersion, num_trajectories, 16, false);

		fn = rp->path + "m2" + cp->fn_suffix;
		save_double_data(fn, m2, num_trajectories, 16, false);

		fn = rp->path + "energy" + cp->fn_suffix;
		save_double_data(fn, energy, num_trajectories, 16, false);
	}
}

void DimerCoreBehaviour::dump_lpn(AllData * ad) const
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	ExpData * ed = ad->ed;

	int num_trajectories = cp->num_trajectories;

	int dump_obs = int(cp->params.find("dump_obs")->second);

	if (dump_obs == 1)
	{
		double * lambda = ed->lambda_now;
		double * mean_lpn = ed->mean_lpn;
		double * energy_lpn = ed->energy_lpn;

		string fn;

		fn = rp->path + "lambda" + cp->fn_suffix;
		save_double_data(fn, lambda, num_trajectories, 16, false);

		int* num_renorms = ed->num_renorms;
		fn = rp->path + "num_renorms" + cp->fn_suffix;
		save_int_data(fn, num_renorms, num_trajectories, false);

		fn = rp->path + "mean_lpn" + cp->fn_suffix;
		save_double_data(fn, mean_lpn, num_trajectories, 16, false);

		fn = rp->path + "energy_lpn" + cp->fn_suffix;
		save_double_data(fn, energy_lpn, num_trajectories, 16, false);

		int save_lambdas = int(ad->cp->params.find("save_lambdas")->second);
		if (save_lambdas > 0)
		{
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				fn = rp->path + "delta_f_by_delta_s_" + to_string(tr_id) + cp->fn_suffix;
				save_double_vector(fn, ed->lambdas[tr_id], 16, false);

				fn = rp->path + "deltas_s_" + to_string(tr_id) + cp->fn_suffix;
				save_double_vector(fn, ed->deltas_s[tr_id], 16, false);

				fn = rp->path + "deltas_f_" + to_string(tr_id) + cp->fn_suffix;
				save_double_vector(fn, ed->deltas_f[tr_id], 16, false);
			}
		}
	}
}

void DimerCoreBehaviour::dump_std_evo(AllData * ad) const
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	ExpData * ed = ad->ed;

	int num_trajectories = cp->num_trajectories;
	int dump_num_total = ed->dump_num_total;

	int dump_obs = int(cp->params.find("dump_obs")->second);
	int jump = int(cp->params.find("jump")->second);

	if (dump_obs == 1)
	{
		int * dump_periods = ed->dump_periods;

		double * norm_evo = ed->norm_evo;
		double * mean_evo = ed->mean_evo;
		double * dispersion_evo = ed->dispersion_evo;
		double * m2_evo = ed->m2_evo;
		double * energy_evo = ed->energy_evo;

		string fn;

		fn = rp->path + "periods" + cp->fn_suffix;
		save_int_data(fn, dump_periods, dump_num_total, false);

		fn = rp->path + "norm_evo" + cp->fn_suffix;
		save_2d_inv_double_data(fn, norm_evo, dump_num_total, num_trajectories, 16, false);

		fn = rp->path + "mean_evo" + cp->fn_suffix;
		save_2d_inv_double_data(fn, mean_evo, dump_num_total, num_trajectories, 16, false);

		fn = rp->path + "dispersion_evo" + cp->fn_suffix;
		save_2d_inv_double_data(fn, dispersion_evo, dump_num_total, num_trajectories, 16, false);

		fn = rp->path + "m2_evo" + cp->fn_suffix;
		save_2d_inv_double_data(fn, m2_evo, dump_num_total, num_trajectories, 16, false);

		fn = rp->path + "energy_evo" + cp->fn_suffix;
		save_2d_inv_double_data(fn, energy_evo, dump_num_total, num_trajectories, 16, false);

		int num_random_obs = int(cp->params.find("num_random_obs")->second);
		int random_obs_seed = int(cp->params.find("random_obs_seed")->second);
		for (int obs_id = 0; obs_id < num_random_obs; obs_id++)
		{
			fn = rp->path + "random_obs_evo_" + to_string(random_obs_seed + obs_id) + cp->fn_suffix;
			std::vector<std::vector<std::complex<double>>> random_obs_evo_to_dump;
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				random_obs_evo_to_dump.push_back(ed->random_obs_evo[tr_id][obs_id]);
			}
			save_2d_vector(fn, random_obs_evo_to_dump, 16);
		}

		if (jump > 0)
		{
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				fn = rp->path + "jump_times_" + to_string(tr_id) + cp->fn_suffix;
				save_double_vector(fn, ed->jump_times[tr_id], 16, false);

				fn = rp->path + "diss_types_" + to_string(tr_id) + cp->fn_suffix;
				save_int_vector(fn, ed->diss_types[tr_id], false);

				if (jump > 1)
				{
					fn = rp->path + "jump_norms_" + to_string(tr_id) + cp->fn_suffix;
					save_double_vector(fn, ed->jump_norms[tr_id], 16, false);

					fn = rp->path + "jump_etas_" + to_string(tr_id) + cp->fn_suffix;
					save_double_vector(fn, ed->jump_etas[tr_id], 16, false);
				}
			}
		}
	}
}

void DimerCoreBehaviour::dump_lpn_evo(AllData * ad) const
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	ExpData * ed = ad->ed;

	int num_trajectories = cp->num_trajectories;
	int dump_num_total = ed->dump_num_total;

	int dump_obs = int(cp->params.find("dump_obs")->second);

	if (dump_obs == 1)
	{
		double * lambda_evo = ed->lambda_evo;
		double * mean_lpn_evo = ed->mean_lpn_evo;
		double * energy_lpn_evo = ed->energy_lpn_evo;

		string fn;

		fn = rp->path + "lambda_evo" + cp->fn_suffix;
		save_2d_inv_double_data(fn, lambda_evo, dump_num_total, num_trajectories, 16, false);

		fn = rp->path + "mean_lpn_evo" + cp->fn_suffix;
		save_2d_inv_double_data(fn, mean_lpn_evo, dump_num_total, num_trajectories, 16, false);

		fn = rp->path + "energy_lpn_evo" + cp->fn_suffix;
		save_2d_inv_double_data(fn, energy_lpn_evo, dump_num_total, num_trajectories, 16, false);

		int num_random_obs = int(cp->params.find("num_random_obs")->second);
		int random_obs_seed = int(cp->params.find("random_obs_seed")->second);
		for (int obs_id = 0; obs_id < num_random_obs; obs_id++)
		{
			fn = rp->path + "random_obs_lpn_evo_" + to_string(random_obs_seed + obs_id) + cp->fn_suffix;
			std::vector<std::vector<std::complex<double>>> random_obs_evo_to_dump;
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				random_obs_evo_to_dump.push_back(ed->random_obs_lpn_evo[tr_id][obs_id]);
			}
			save_2d_vector(fn, random_obs_evo_to_dump, 16);
		}
	}
}










void DimerSyncCoreBehaviour::init_splits(AllData* ad) const
{
	RunParam* rp = ad->rp;
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;

	int num_threads = rp->num_threads;

	md->structure = init_split_structure_dimersync(ad);
	md->splits = new Split[num_threads];
	for (int th_id = 0; th_id < num_threads; th_id++)
	{
		copy_struct_not_member(md->structure, &(md->splits)[th_id]);
	}

	cout << "Split initialization complete" << endl;
}

void DimerSyncCoreBehaviour::free_splits(AllData* ad) const
{
	RunParam* rp = ad->rp;
	MainData* md = ad->md;

	int num_threads = rp->num_threads;

	for (int i = 0; i < num_threads; i++)
	{
		delete_split_struct_not_member(&(md->splits[i]));
	}

	delete(md->splits);
	delete_split_struct(md->structure);
}

void DimerSyncCoreBehaviour::init_splits_deep(AllData* ad) const
{
	RunParam* rp = ad->rp;
	MainData* md = ad->md;

	int num_branches = md->num_ham_qj;
	int num_threads = rp->num_threads;

	int num_total = num_threads * num_branches;

	md->structure = init_split_structure_dimersync_deep(ad);
	md->splits = new Split[num_total];

	for (int b_id = 0; b_id < num_branches; b_id++)
	{
		for (int th_id = 0; th_id < num_threads; th_id++)
		{
			int index = b_id * num_threads + th_id;
			copy_struct_not_member(&(md->structure)[b_id], &(md->splits)[index]);
		}
	}
}

void DimerSyncCoreBehaviour::free_splits_deep(AllData* ad) const
{
	RunParam* rp = ad->rp;
	MainData* md = ad->md;

	int num_branches = md->num_ham_qj;
	int num_threads = rp->num_threads;

	for (int b_id = 0; b_id < num_branches; b_id++)
	{
		for (int th_id = 0; th_id < num_threads; th_id++)
		{
			int index = b_id * num_threads + th_id;
			delete_split_struct_not_member(&(md->splits[index]));
		}
	}

	delete(md->splits);
	delete_split_struct(md->structure);
}

void DimerSyncCoreBehaviour::ex_period(AllData* ad, int tr_id, int th_id, int period_id) const
{
	MainData* md = ad->md;

	int sys_size = md->sys_size;
	Split* head = &(md->splits)[th_id];

	for (int b_id = 0; b_id < head->counter; b_id++)
	{
		one_period_branch(ad, head, tr_id, &(head->next)[b_id]);
	}
}

void DimerSyncCoreBehaviour::ex_period_trp_deep(AllData* ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;

	int num_branches = md->num_ham_qj;
	int num_sub_steps = int(cp->params.find("deep_num_steps")->second);

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		for (int sub_step_id = 0; sub_step_id < num_sub_steps; sub_step_id++)
		{
			one_sub_period_deep(ad, tr_id, part_id, th_id);
		}
	}
}

void DimerSyncCoreBehaviour::ex_period_obs_deep(AllData* ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

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

			one_sub_period_deep(ad, tr_id, part_id, th_id);
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

void DimerSyncCoreBehaviour::ex_period_obs_deep_lpn(AllData* ad, int period_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);

	int num_trajectories = cp->num_trajectories;

	int num_branches = md->num_ham_qj;
	int num_sub_steps_per_part = int(cp->params.find("deep_num_steps")->second);
	int num_sub_steps = num_branches * int(cp->params.find("deep_num_steps")->second);

	int dump_point_id = 0;
	int global_point_id = 0;

	int step_id = 0;

	CoreBehavior* tmp = new DimerSyncCoreBehaviour;

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		for (int sub_step_id = 0; sub_step_id < num_sub_steps_per_part; sub_step_id++)
		{
			step_id = part_id * num_sub_steps_per_part + sub_step_id;

			global_point_id = period_id * num_sub_steps + dump_point_id;
			dump_point_id++;
			int dump_id = global_point_id + 1;

			ed->curr_time = double(dump_id) / double(num_sub_steps) * md->T;

			one_sub_period_deep(ad, 0, part_id, 0);
			calc_chars_std(ad, 0);
			calc_chars_lpn(ad, 0, 0);
			evo_chars_std(ad, 0, dump_id);
			evo_chars_lpn(ad, 0, dump_id);

#pragma omp parallel for
			for (int tr_id = 1; tr_id < num_trajectories; tr_id++)
			{
				int thread_id = omp_get_thread_num();
				one_sub_period_deep(ad, tr_id, part_id, thread_id);
				calc_chars_std(ad, tr_id);
				lambda_lpn(ad, tmp, tr_id, 0);
				evo_chars_std(ad, tr_id, dump_id);
				evo_chars_lpn(ad, tr_id, dump_id);
			}

#pragma omp parallel for
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				if (dump_evo_sep == 1)
				{
					dump_adr_single(ad, tr_id, true);
				}
			}
		}
	}

	delete tmp;
}

void DimerSyncCoreBehaviour::ex_period_obs_deep_mult_lpn(AllData* ad, int period_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);

	int num_trajectories = cp->num_trajectories;

	int num_branches = md->num_ham_qj;
	int num_sub_steps_per_part = int(cp->params.find("deep_num_steps")->second);
	int num_sub_steps = num_branches * int(cp->params.find("deep_num_steps")->second);

	int lambda_per_periods = int(cp->params.find("lambda_per_periods")->second);

	int dump_point_id = 0;
	int global_point_id = 0;

	int step_id = 0;

	CoreBehavior* tmp = new DimerSyncCoreBehaviour;

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		for (int sub_step_id = 0; sub_step_id < num_sub_steps_per_part; sub_step_id++)
		{
			step_id = part_id * num_sub_steps_per_part + sub_step_id;

			global_point_id = period_id * num_sub_steps + dump_point_id;
			dump_point_id++;
			int dump_id = global_point_id + 1;

			ed->curr_time = double(dump_id) / double(num_sub_steps) * md->T;

#pragma omp parallel for
			for (int tr_id = 0; tr_id < num_trajectories / 2; tr_id++)
			{
				int thread_id = omp_get_thread_num();
				one_sub_period_deep(ad, tr_id, part_id, thread_id);
				calc_chars_std(ad, tr_id);
				calc_chars_lpn(ad, tr_id, tr_id);
				evo_chars_std(ad, tr_id, dump_id);
				evo_chars_lpn(ad, tr_id, dump_id);
			}

#pragma omp parallel for
			for (int tr_id = num_trajectories / 2; tr_id < num_trajectories; tr_id++)
			{
				int thread_id = omp_get_thread_num();
				one_sub_period_deep(ad, tr_id, part_id, thread_id);
				calc_chars_std(ad, tr_id);
				if (lambda_per_periods == 0)
				{
					lambda_lpn(ad, tmp, tr_id, tr_id - num_trajectories / 2);
				}
				if (lambda_per_periods == 1)
				{
					if ((part_id == num_branches - 1) && (sub_step_id == num_sub_steps_per_part - 1))
					{
						lambda_lpn_now(ad, tmp, tr_id, tr_id - num_trajectories / 2);
					}
				}
				evo_chars_std(ad, tr_id, dump_id);
				evo_chars_lpn(ad, tr_id, dump_id);
			}

#pragma omp parallel for
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				if (dump_evo_sep == 1)
				{
					dump_adr_single(ad, tr_id, true);
				}
			}
		}
	}

	delete tmp;
}

void DimerSyncCoreBehaviour::ex_period_obs_deep_lpn_per_period(struct AllData* ad, int period_id, int num_periods) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);

	int num_trajectories = cp->num_trajectories;

	int num_branches = md->num_ham_qj;
	int num_sub_steps_per_part = int(cp->params.find("deep_num_steps")->second);
	int num_sub_steps = num_branches * int(cp->params.find("deep_num_steps")->second);

	int dump_point_id = 0;
	int global_point_id = 0;

	int step_id = 0;

	CoreBehavior* tmp = new DimerSyncCoreBehaviour;

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		for (int sub_step_id = 0; sub_step_id < num_sub_steps_per_part; sub_step_id++)
		{
			step_id = part_id * num_sub_steps_per_part + sub_step_id;

			global_point_id = period_id * num_sub_steps + dump_point_id;
			dump_point_id++;
			int dump_id = global_point_id + 1;

			ed->curr_time = double(dump_id) / double(num_sub_steps) * md->T;

			one_sub_period_deep(ad, 0, part_id, 0);
			calc_chars_std(ad, 0);
			calc_chars_lpn(ad, 0, 0);
			evo_chars_std(ad, 0, dump_id);
			evo_chars_lpn(ad, 0, dump_id);

#pragma omp parallel for
			for (int tr_id = 1; tr_id < num_trajectories; tr_id++)
			{
				int thread_id = omp_get_thread_num();
				one_sub_period_deep(ad, tr_id, part_id, thread_id);
				calc_chars_std(ad, tr_id);
				lambda_lpn_per_periods(ad, tmp, tr_id, 0, num_sub_steps, global_point_id, num_periods);
				evo_chars_std(ad, tr_id, dump_id);
				evo_chars_lpn(ad, tr_id, dump_id);
			}

#pragma omp parallel for
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				if (dump_evo_sep == 1)
				{
					dump_adr_single(ad, tr_id, true);
				}
			}
		}
	}

	delete tmp;
}

void DimerSyncCoreBehaviour::ex_period_obs_deep_cd(AllData* ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

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

			one_sub_period_deep(ad, tr_id, part_id, th_id);
			calc_chars_std(ad, tr_id);
		}
	}
}

void DimerSyncCoreBehaviour::ex_period_obs_deep_sigma(AllData* ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;

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

			one_sub_period_deep(ad, tr_id, part_id, th_id);
			calc_chars_std(ad, tr_id);
		}
	}
}

void DimerSyncCoreBehaviour::rk_period(AllData* ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	for (int step_id = 0; step_id < cp->rk_ns; step_id++)
	{
		ed->times_all[tr_id] = double(period_id) * md->T + double(step_id) * ed->rk_step;
		rk_step_dimersync(ad, tr_id, th_id, ed->rk_step);
	}
}

void DimerSyncCoreBehaviour::rk_period_trp_deep(AllData* ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int num_branches = md->num_ham_qj;
	int num_sub_steps_per_part = int(cp->params.find("deep_num_steps")->second);

	int step_id = 0;

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		for (int sub_step_id = 0; sub_step_id < num_sub_steps_per_part; sub_step_id++)
		{
			step_id = part_id * num_sub_steps_per_part + sub_step_id;

			for (int in_step_id = 0; in_step_id < cp->rk_ns; in_step_id++)
			{
				ed->times_all[tr_id] = double(period_id) * md->T + double(step_id * cp->rk_ns + in_step_id) * ed->rk_step;
				rk_step_dimersync(ad, tr_id, th_id, ed->rk_step);
			}
		}
	}
}

void DimerSyncCoreBehaviour::rk_period_obs_deep(AllData* ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

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

			for (int in_step_id = 0; in_step_id < cp->rk_ns; in_step_id++)
			{
				ed->times_all[tr_id] = double(period_id) * md->T + double(step_id * cp->rk_ns + in_step_id) * ed->rk_step;
				rk_step_dimersync(ad, tr_id, th_id, ed->rk_step);
			}

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

void DimerSyncCoreBehaviour::rk_period_obs_deep_lpn(AllData* ad, int period_id) const
{
	stringstream msg;
	msg << "Error: need to add functionality" << endl;
	Error(msg.str());
}

void DimerSyncCoreBehaviour::rk_period_obs_deep_cd(AllData* ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

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

			for (int in_step_id = 0; in_step_id < cp->rk_ns; in_step_id++)
			{
				ed->times_all[tr_id] = double(period_id) * md->T + double(step_id * cp->rk_ns + in_step_id) * ed->rk_step;
				rk_step_dimersync(ad, tr_id, th_id, ed->rk_step);
			}

			calc_chars_std(ad, tr_id);
		}
	}
}

void DimerSyncCoreBehaviour::rk_period_obs_deep_sigma(AllData* ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

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

			dump_point_id++;

			for (int in_step_id = 0; in_step_id < cp->rk_ns; in_step_id++)
			{
				ed->times_all[tr_id] = double(period_id) * md->T + double(step_id * cp->rk_ns + in_step_id) * ed->rk_step;
				rk_step_dimersync(ad, tr_id, th_id, ed->rk_step);
			}

			calc_chars_std(ad, tr_id);
		}
	}
}

void DimerSyncCoreBehaviour::calc_chars_std_start(AllData* ad, int tr_id) const
{
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int sys_size = md->sys_size;
	MKL_Complex16* phi = &(ed->phi_all[tr_id * sys_size]);
	double* adr = &(ed->abs_diag_rho_all[tr_id * sys_size]);

	double norm = norm_square(phi, sys_size);

	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		adr[st_id] = mult_scalar_double(mult_scalar_complex(&phi[st_id], &phi[st_id], 1), 1.0 / norm).real;
	}

	double mean = get_mean_simple(adr, sys_size);
	double dispersion = get_dispersion_simple(mean, mean);
	double m2 = get_m2(adr, sys_size, mean);
	double energy = get_energy_sync(ad, tr_id);

	ed->mean_start[tr_id] = mean;

	ed->norm[tr_id] = norm;
	ed->mean[tr_id] = mean;
	ed->dispersion[tr_id] = dispersion;
	ed->m2[tr_id] = m2;
	ed->energy[tr_id] = energy;

	calc_random_obs(ad, tr_id);
}

void DimerSyncCoreBehaviour::calc_chars_std(AllData* ad, int tr_id) const
{
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int sys_size = md->sys_size;
	MKL_Complex16* phi = &(ed->phi_all[tr_id * sys_size]);
	double* adr = &(ed->abs_diag_rho_all[tr_id * sys_size]);

	double norm = norm_square(phi, sys_size);

	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		adr[st_id] = mult_scalar_double(mult_scalar_complex(&phi[st_id], &phi[st_id], 1), 1.0 / norm).real;
	}

	double mean = get_mean_simple(adr, sys_size);
	double dispersion = get_dispersion_simple(mean, ed->mean_start[tr_id]);
	double m2 = get_m2(adr, sys_size, mean);
	double energy = get_energy_sync(ad, tr_id);

	ed->norm[tr_id] = norm;
	ed->mean[tr_id] = mean;
	ed->dispersion[tr_id] = dispersion;
	ed->m2[tr_id] = m2;
	ed->energy[tr_id] = energy;

	calc_random_obs(ad, tr_id);
}

void DimerSyncCoreBehaviour::calc_chars_lpn_start(AllData* ad, int tr_id, int base_tr_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int sys_size = md->sys_size;

	MKL_Complex16* phi = &(ed->phi_all[tr_id * sys_size]);
	double* adr = &(ed->abs_diag_rho_all[tr_id * sys_size]);

	double norm_2 = norm_square(phi, sys_size);

	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		adr[st_id] = mult_scalar_double(mult_scalar_complex(&phi[st_id], &phi[st_id], 1), 1.0 / norm_2).real;
	}

	double lambda = 0.0;
	double lambda_now = 0.0;
	double mean_lpn = ed->mean[tr_id];
	double energy_lpn = ed->energy[tr_id];

	double delta_s = this->calc_delta_f(ad, tr_id, base_tr_id); // Important! Here we use calc_delta_f not calc_delta_s

	ed->lambda[tr_id] = lambda;
	ed->lambda_now[tr_id] = lambda_now;
	ed->mean_lpn[tr_id] = mean_lpn;
	ed->energy_lpn[tr_id] = energy_lpn;
	ed->delta_s[tr_id] = delta_s;

	calc_random_obs_lpn_start(ad, tr_id);
}

void DimerSyncCoreBehaviour::calc_chars_lpn(AllData* ad, int tr_id, int base_tr_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int sys_size = md->sys_size;

	MKL_Complex16* phi = &(ed->phi_all[tr_id * sys_size]);
	double* adr = &(ed->abs_diag_rho_all[tr_id * sys_size]);

	double norm_2 = norm_square(phi, sys_size);

	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		adr[st_id] = mult_scalar_double(mult_scalar_complex(&phi[st_id], &phi[st_id], 1), 1.0 / norm_2).real;
	}

	double mean_lpn = get_mean_simple(adr, sys_size);
	double energy_lpn = get_energy_sync(ad, tr_id);

	ed->mean_lpn[tr_id] = mean_lpn;
	ed->energy_lpn[tr_id] = energy_lpn;

	calc_random_obs_lpn(ad, tr_id);
}

double DimerSyncCoreBehaviour::calc_T(AllData* ad) const
{
	MainData* md = ad->md;
	double T = md->T;
	return T;
}

void DimerSyncCoreBehaviour::evo_chars_std(AllData* ad, int tr_id, int dump_id) const
{
	ConfigParam* cp = ad->cp;
	ExpData* ed = ad->ed;

	int num_trajectories = cp->num_trajectories;
	int dump_num_total = ed->dump_num_total;

	int index = tr_id * dump_num_total + dump_id;

	ed->norm_evo[index] = ed->norm[tr_id];
	ed->mean_evo[index] = ed->mean[tr_id];
	ed->dispersion_evo[index] = ed->dispersion[tr_id];
	ed->m2_evo[index] = ed->m2[tr_id];
	ed->energy_evo[index] = ed->energy[tr_id];

	int num_random_obs = int(cp->params.find("num_random_obs")->second);
	for (int obs_id = 0; obs_id < num_random_obs; obs_id++)
	{
		ed->random_obs_evo[tr_id][obs_id].push_back(ed->random_obs[tr_id][obs_id]);
	}
}

void DimerSyncCoreBehaviour::evo_chars_lpn(AllData* ad, int tr_id, int dump_id) const
{
	ConfigParam* cp = ad->cp;
	ExpData* ed = ad->ed;

	int num_trajectories = cp->num_trajectories;
	int dump_num_total = ed->dump_num_total;

	int index = tr_id * dump_num_total + dump_id;

	ed->lambda_evo[index] = ed->lambda_now[tr_id];
	ed->mean_lpn_evo[index] = ed->mean_lpn[tr_id];
	ed->energy_lpn_evo[index] = ed->energy_lpn[tr_id];

	int num_random_obs = int(cp->params.find("num_random_obs")->second);
	for (int obs_id = 0; obs_id < num_random_obs; obs_id++)
	{
		ed->random_obs_lpn_evo[tr_id][obs_id].push_back(ed->random_obs_lpn[tr_id][obs_id]);
	}
}

double DimerSyncCoreBehaviour::calc_delta_s(AllData* ad, int tr_id, int base_tr_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int sys_size = md->sys_size;

	int lpn_type = int(cp->params.find("lpn_type")->second);
	int num_random_obs = int(ad->cp->params.find("num_random_obs")->second);

	double delta_s = 0.0;
	if (lpn_type == -1)
	{
		delta_s = fabs(ed->mean_lpn[tr_id] - ed->mean_lpn[base_tr_id]) / double(sys_size);
	}
	else if (lpn_type == -2)
	{
		delta_s = fabs(ed->energy_lpn[tr_id] - ed->energy_lpn[base_tr_id]) / ed->max_energy;
	}
	else if (lpn_type >= 0 && lpn_type < num_random_obs)
	{
		complex<double> base = ed->random_obs_lpn[base_tr_id][lpn_type];
		complex<double> var = ed->random_obs_lpn[tr_id][lpn_type];
		double tmp = pow((base.real() - var.real()), 2) + pow((base.imag() - var.imag()), 2);
		delta_s = sqrt(tmp);
	}
	else
	{
		stringstream msg;
		msg << "Error: Wrong lpn_type: " << lpn_type << endl;
		Error(msg.str());
	}

	return delta_s;
}

double DimerSyncCoreBehaviour::calc_delta_f(AllData* ad, int tr_id, int base_tr_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int sys_size = md->sys_size;

	int lpn_type = int(cp->params.find("lpn_type")->second);
	int num_random_obs = int(ad->cp->params.find("num_random_obs")->second);

	double delta_f = 0.0;
	if (lpn_type == -1)
	{
		delta_f = fabs(ed->mean[tr_id] - ed->mean[base_tr_id]) / double(sys_size);
	}
	else if (lpn_type == -2)
	{
		delta_f = fabs(ed->energy[tr_id] - ed->energy[base_tr_id]) / ed->max_energy;
	}
	else if (lpn_type >= 0 && lpn_type < num_random_obs)
	{
		complex<double> base = ed->random_obs[base_tr_id][lpn_type];
		complex<double> var = ed->random_obs[tr_id][lpn_type];
		double tmp = pow((base.real() - var.real()), 2) + pow((base.imag() - var.imag()), 2);
		delta_f = sqrt(tmp);
	}
	else
	{
		stringstream msg;
		msg << "Error: Wrong lpn_type: " << lpn_type << endl;
		Error(msg.str());
	}

	return delta_f;
}

void DimerSyncCoreBehaviour::calc_ci(AllData* ad, int tr_id) const
{
	calc_ci_double(ad, tr_id);
}

void DimerSyncCoreBehaviour::dump_std(AllData* ad) const
{
	RunParam* rp = ad->rp;
	ConfigParam* cp = ad->cp;
	ExpData* ed = ad->ed;

	int dump_obs = int(cp->params.find("dump_obs")->second);

	if (dump_obs == 1)
	{
		int num_trajectories = cp->num_trajectories;

		double* norm = ed->norm;
		double* mean_start = ed->mean_start;
		double* mean = ed->mean;
		double* dispersion = ed->dispersion;
		double* m2 = ed->m2;
		double* energy = ed->energy;

		string fn;

		fn = rp->path + "norm" + cp->fn_suffix;
		save_double_data(fn, norm, num_trajectories, 16, false);

		fn = rp->path + "mean_start" + cp->fn_suffix;
		save_double_data(fn, mean_start, num_trajectories, 16, false);

		fn = rp->path + "mean" + cp->fn_suffix;
		save_double_data(fn, mean, num_trajectories, 16, false);

		fn = rp->path + "dispersion" + cp->fn_suffix;
		save_double_data(fn, dispersion, num_trajectories, 16, false);

		fn = rp->path + "m2" + cp->fn_suffix;
		save_double_data(fn, m2, num_trajectories, 16, false);

		fn = rp->path + "energy" + cp->fn_suffix;
		save_double_data(fn, energy, num_trajectories, 16, false);
	}
}

void DimerSyncCoreBehaviour::dump_lpn(AllData* ad) const
{
	RunParam* rp = ad->rp;
	ConfigParam* cp = ad->cp;
	ExpData* ed = ad->ed;

	int num_trajectories = cp->num_trajectories;

	int dump_obs = int(cp->params.find("dump_obs")->second);

	string fn;

	double* lambda = ed->lambda_now;
	fn = rp->path + "lambda" + cp->fn_suffix;
	save_double_data(fn, lambda, num_trajectories, 16, false);

	if (dump_obs == 1)
	{
		int* num_renorms = ed->num_renorms;
		fn = rp->path + "num_renorms" + cp->fn_suffix;
		save_int_data(fn, num_renorms, num_trajectories, false);

		double* mean_lpn = ed->mean_lpn;
		fn = rp->path + "mean_lpn" + cp->fn_suffix;
		save_double_data(fn, mean_lpn, num_trajectories, 16, false);

		double* energy_lpn = ed->energy_lpn;
		fn = rp->path + "energy_lpn" + cp->fn_suffix;
		save_double_data(fn, energy_lpn, num_trajectories, 16, false);

		int save_lambdas = int(ad->cp->params.find("save_lambdas")->second);
		if (save_lambdas > 0)
		{
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				fn = rp->path + "delta_f_by_delta_s_" + to_string(tr_id) + cp->fn_suffix;
				save_double_vector(fn, ed->lambdas[tr_id], 16, false);

				fn = rp->path + "deltas_s_" + to_string(tr_id) + cp->fn_suffix;
				save_double_vector(fn, ed->deltas_s[tr_id], 16, false);

				fn = rp->path + "deltas_f_" + to_string(tr_id) + cp->fn_suffix;
				save_double_vector(fn, ed->deltas_f[tr_id], 16, false);
			}
		}
	}
}

void DimerSyncCoreBehaviour::dump_std_evo(AllData* ad) const
{
	RunParam* rp = ad->rp;
	ConfigParam* cp = ad->cp;
	ExpData* ed = ad->ed;

	int num_trajectories = cp->num_trajectories;
	int dump_num_total = ed->dump_num_total;

	int dump_obs = int(cp->params.find("dump_obs")->second);
	int jump = int(cp->params.find("jump")->second);

	if (dump_obs == 1)
	{
		int* dump_periods = ed->dump_periods;

		double* norm_evo = ed->norm_evo;
		double* mean_evo = ed->mean_evo;
		double* dispersion_evo = ed->dispersion_evo;
		double* m2_evo = ed->m2_evo;
		double* energy_evo = ed->energy_evo;

		string fn;

		fn = rp->path + "periods" + cp->fn_suffix;
		save_int_data(fn, dump_periods, dump_num_total, false);

		fn = rp->path + "norm_evo" + cp->fn_suffix;
		save_2d_inv_double_data(fn, norm_evo, dump_num_total, num_trajectories, 16, false);

		fn = rp->path + "mean_evo" + cp->fn_suffix;
		save_2d_inv_double_data(fn, mean_evo, dump_num_total, num_trajectories, 16, false);

		fn = rp->path + "dispersion_evo" + cp->fn_suffix;
		save_2d_inv_double_data(fn, dispersion_evo, dump_num_total, num_trajectories, 16, false);

		fn = rp->path + "m2_evo" + cp->fn_suffix;
		save_2d_inv_double_data(fn, m2_evo, dump_num_total, num_trajectories, 16, false);

		fn = rp->path + "energy_evo" + cp->fn_suffix;
		save_2d_inv_double_data(fn, energy_evo, dump_num_total, num_trajectories, 16, false);

		int num_random_obs = int(cp->params.find("num_random_obs")->second);
		int random_obs_seed = int(cp->params.find("random_obs_seed")->second);
		for (int obs_id = 0; obs_id < num_random_obs; obs_id++)
		{
			fn = rp->path + "random_obs_evo_" + to_string(random_obs_seed + obs_id) + cp->fn_suffix;
			std::vector<std::vector<std::complex<double>>> random_obs_evo_to_dump;
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				random_obs_evo_to_dump.push_back(ed->random_obs_evo[tr_id][obs_id]);
			}
			save_2d_vector(fn, random_obs_evo_to_dump, 16);
		}

		if (jump > 0)
		{
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				fn = rp->path + "jump_times_" + to_string(tr_id) + cp->fn_suffix;
				save_double_vector(fn, ed->jump_times[tr_id], 16, false);

				fn = rp->path + "diss_types_" + to_string(tr_id) + cp->fn_suffix;
				save_int_vector(fn, ed->diss_types[tr_id], false);

				if (jump > 1)
				{
					fn = rp->path + "jump_norms_" + to_string(tr_id) + cp->fn_suffix;
					save_double_vector(fn, ed->jump_norms[tr_id], 16, false);

					fn = rp->path + "jump_etas_" + to_string(tr_id) + cp->fn_suffix;
					save_double_vector(fn, ed->jump_etas[tr_id], 16, false);
				}
			}
		}
	}
}

void DimerSyncCoreBehaviour::dump_lpn_evo(AllData* ad) const
{
	RunParam* rp = ad->rp;
	ConfigParam* cp = ad->cp;
	ExpData* ed = ad->ed;

	int num_trajectories = cp->num_trajectories;
	int dump_num_total = ed->dump_num_total;

	int dump_obs = int(cp->params.find("dump_obs")->second);

	if (dump_obs == 1)
	{
		double* lambda_evo = ed->lambda_evo;
		double* mean_lpn_evo = ed->mean_lpn_evo;
		double* energy_lpn_evo = ed->energy_lpn_evo;

		string fn;

		fn = rp->path + "lambda_evo" + cp->fn_suffix;
		save_2d_inv_double_data(fn, lambda_evo, dump_num_total, num_trajectories, 16, false);

		fn = rp->path + "mean_lpn_evo" + cp->fn_suffix;
		save_2d_inv_double_data(fn, mean_lpn_evo, dump_num_total, num_trajectories, 16, false);

		fn = rp->path + "energy_lpn_evo" + cp->fn_suffix;
		save_2d_inv_double_data(fn, energy_lpn_evo, dump_num_total, num_trajectories, 16, false);

		int num_random_obs = int(cp->params.find("num_random_obs")->second);
		int random_obs_seed = int(cp->params.find("random_obs_seed")->second);
		for (int obs_id = 0; obs_id < num_random_obs; obs_id++)
		{
			fn = rp->path + "random_obs_lpn_evo_" + to_string(random_obs_seed + obs_id) + cp->fn_suffix;
			std::vector<std::vector<std::complex<double>>> random_obs_evo_to_dump;
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				random_obs_evo_to_dump.push_back(ed->random_obs_lpn_evo[tr_id][obs_id]);
			}
			save_2d_vector(fn, random_obs_evo_to_dump, 16);
		}
	}
}










void JCSCoreBehaviour::init_splits(AllData * ad) const
{
	RunParam * rp = ad->rp;
	MainData * md = ad->md;

	int num_branches = md->num_ham_qj;
	int num_threads = rp->num_threads;

	int num_total = num_threads * num_branches;

	md->structure = init_split_structure_jcs(ad);
	md->splits = new Split[num_total];

	for (int b_id = 0; b_id < num_branches; b_id++)
	{
		for (int th_id = 0; th_id < num_threads; th_id++)
		{
			int index = b_id * num_threads + th_id;
			copy_struct_not_member(&(md->structure)[b_id], &(md->splits)[index]);
		}
	}
}

void JCSCoreBehaviour::free_splits(AllData * ad) const
{
	RunParam * rp = ad->rp;
	MainData * md = ad->md;

	int num_branches = md->num_ham_qj;
	int num_threads = rp->num_threads;

	for (int b_id = 0; b_id < num_branches; b_id++)
	{
		for (int th_id = 0; th_id < num_threads; th_id++)
		{
			int index = b_id * num_threads + th_id;
			delete_split_struct_not_member(&(md->splits[index]));
		}
	}

	delete(md->splits);
	delete_split_struct(md->structure);
}

void JCSCoreBehaviour::init_splits_deep(AllData * ad) const
{
	RunParam * rp = ad->rp;
	MainData * md = ad->md;

	int num_branches = md->num_ham_qj;
	int num_threads = rp->num_threads;

	int num_total = num_threads * num_branches;

	md->structure = init_split_structure_jcs_deep(ad);
	md->splits = new Split[num_total];

	for (int b_id = 0; b_id < num_branches; b_id++)
	{
		for (int th_id = 0; th_id < num_threads; th_id++)
		{
			int index = b_id * num_threads + th_id;
			copy_struct_not_member(&(md->structure)[b_id], &(md->splits)[index]);
		}
	}
}

void JCSCoreBehaviour::free_splits_deep(AllData * ad) const
{
	RunParam * rp = ad->rp;
	MainData * md = ad->md;

	int num_branches = md->num_ham_qj;
	int num_threads = rp->num_threads;

	for (int b_id = 0; b_id < num_branches; b_id++)
	{
		for (int th_id = 0; th_id < num_threads; th_id++)
		{
			int index = b_id * num_threads + th_id;
			delete_split_struct_not_member(&(md->splits[index]));
		}
	}

	delete(md->splits);
	delete_split_struct(md->structure);
}

void JCSCoreBehaviour::ex_period(AllData * ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;

	int num_branches = md->num_ham_qj;

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		one_sub_period_deep(ad, tr_id, part_id, th_id);
	}
}

void JCSCoreBehaviour::ex_period_trp_deep(AllData * ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;

	int num_branches = md->num_ham_qj;
	int num_sub_steps = int(cp->params.find("deep_num_steps")->second);

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		for (int sub_step_id = 0; sub_step_id < num_sub_steps; sub_step_id++)
		{
			one_sub_period_deep(ad, tr_id, part_id, th_id);
		}
	}
}

void JCSCoreBehaviour::ex_period_obs_deep(AllData * ad, int tr_id, int th_id, int period_id) const
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

			one_sub_period_deep(ad, tr_id, part_id, th_id);
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

void JCSCoreBehaviour::ex_period_obs_deep_lpn(AllData * ad, int period_id) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	double jcs_drv_part_1 = double(cp->params.find("jcs_drv_part_1")->second);
	double jcs_drv_part_2 = double(cp->params.find("jcs_drv_part_2")->second);
	double T_1 = jcs_drv_part_1 * md->T;
	double T_2 = jcs_drv_part_2 * md->T;
	double T = T_1 + T_2;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);

	int num_trajectories = cp->num_trajectories;

	int num_branches = md->num_ham_qj;
	int num_sub_steps_per_part = int(cp->params.find("deep_num_steps")->second);
	int num_sub_steps = num_branches * int(cp->params.find("deep_num_steps")->second);

	int dump_point_id = 0;
	int global_point_id = 0;

	int step_id = 0;

	CoreBehavior * tmp = new JCSCoreBehaviour;

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		for (int sub_step_id = 0; sub_step_id < num_sub_steps_per_part; sub_step_id++)
		{
			step_id = part_id * num_sub_steps_per_part + sub_step_id;

			global_point_id = period_id * num_sub_steps + dump_point_id;
			dump_point_id++;
			int dump_id = global_point_id + 1;

			ed->curr_time = double(dump_id) / double(num_sub_steps) * T;

			one_sub_period_deep(ad, 0, part_id, 0);
			calc_chars_std(ad, 0);
			calc_chars_lpn(ad, 0, 0);
			evo_chars_std(ad, 0, dump_id);
			evo_chars_lpn(ad, 0, dump_id);

#pragma omp parallel for
			for (int tr_id = 1; tr_id < num_trajectories; tr_id++)
			{
				int thread_id = omp_get_thread_num();
				one_sub_period_deep(ad, tr_id, part_id, thread_id);
				calc_chars_std(ad, tr_id);
				lambda_lpn(ad, tmp, tr_id,0);
				evo_chars_std(ad, tr_id, dump_id);
				evo_chars_lpn(ad, tr_id, dump_id);
			}

#pragma omp parallel for
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				if (dump_evo_sep == 1)
				{
					dump_adr_single(ad, tr_id, true);
				}
			}
		}
	}

	delete tmp;
}

void JCSCoreBehaviour::ex_period_obs_deep_mult_lpn(AllData * ad, int period_id) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	double jcs_drv_part_1 = double(cp->params.find("jcs_drv_part_1")->second);
	double jcs_drv_part_2 = double(cp->params.find("jcs_drv_part_2")->second);
	double T_1 = jcs_drv_part_1 * md->T;
	double T_2 = jcs_drv_part_2 * md->T;
	double T = T_1 + T_2;

	int lambda_per_periods = int(cp->params.find("lambda_per_periods")->second);

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);

	int num_trajectories = cp->num_trajectories;

	int num_branches = md->num_ham_qj;
	int num_sub_steps_per_part = int(cp->params.find("deep_num_steps")->second);
	int num_sub_steps = num_branches * int(cp->params.find("deep_num_steps")->second);

	int dump_point_id = 0;
	int global_point_id = 0;

	int step_id = 0;

	CoreBehavior * tmp = new JCSCoreBehaviour;

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		for (int sub_step_id = 0; sub_step_id < num_sub_steps_per_part; sub_step_id++)
		{
			step_id = part_id * num_sub_steps_per_part + sub_step_id;

			global_point_id = period_id * num_sub_steps + dump_point_id;
			dump_point_id++;
			int dump_id = global_point_id + 1;

			ed->curr_time = double(dump_id) / double(num_sub_steps) * T;

#pragma omp parallel for
			for (int tr_id = 0; tr_id < num_trajectories / 2; tr_id++)
			{
				int thread_id = omp_get_thread_num();
				one_sub_period_deep(ad, tr_id, part_id, thread_id);
				calc_chars_std(ad, tr_id);
				calc_chars_lpn(ad, tr_id, tr_id);
				evo_chars_std(ad, tr_id, dump_id);
				evo_chars_lpn(ad, tr_id, dump_id);
			}

#pragma omp parallel for
			for (int tr_id = num_trajectories / 2; tr_id < num_trajectories; tr_id++)
			{
				int thread_id = omp_get_thread_num();
				one_sub_period_deep(ad, tr_id, part_id, thread_id);
				calc_chars_std(ad, tr_id);
				if (lambda_per_periods == 0)
				{
					lambda_lpn(ad, tmp, tr_id, tr_id - num_trajectories / 2);
				}
				if (lambda_per_periods == 1)
				{
					if ((part_id == num_branches - 1) && (sub_step_id == num_sub_steps_per_part - 1))
					{
						lambda_lpn_now(ad, tmp, tr_id, tr_id - num_trajectories / 2);
					}
				}
				evo_chars_std(ad, tr_id, dump_id);
				evo_chars_lpn(ad, tr_id, dump_id);
			}

#pragma omp parallel for
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				if (dump_evo_sep == 1)
				{
					dump_adr_single(ad, tr_id, true);
				}
			}
		}
	}

	delete tmp;
}

void JCSCoreBehaviour::ex_period_obs_deep_lpn_per_period(struct AllData* ad, int period_id, int num_periods) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	double jcs_drv_part_1 = double(cp->params.find("jcs_drv_part_1")->second);
	double jcs_drv_part_2 = double(cp->params.find("jcs_drv_part_2")->second);
	double T_1 = jcs_drv_part_1 * md->T;
	double T_2 = jcs_drv_part_2 * md->T;
	double T = T_1 + T_2;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);

	int num_trajectories = cp->num_trajectories;

	int num_branches = md->num_ham_qj;
	int num_sub_steps_per_part = int(cp->params.find("deep_num_steps")->second);
	int num_sub_steps = num_branches * int(cp->params.find("deep_num_steps")->second);

	int dump_point_id = 0;
	int global_point_id = 0;

	int step_id = 0;

	CoreBehavior * tmp = new JCSCoreBehaviour;

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		for (int sub_step_id = 0; sub_step_id < num_sub_steps_per_part; sub_step_id++)
		{
			step_id = part_id * num_sub_steps_per_part + sub_step_id;

			global_point_id = period_id * num_sub_steps + dump_point_id;
			dump_point_id++;
			int dump_id = global_point_id + 1;

			ed->curr_time = double(dump_id) / double(num_sub_steps) * T;

			one_sub_period_deep(ad, 0, part_id, 0);
			calc_chars_std(ad, 0);
			calc_chars_lpn(ad, 0, 0);
			evo_chars_std(ad, 0, dump_id);
			evo_chars_lpn(ad, 0, dump_id);

#pragma omp parallel for
			for (int tr_id = 1; tr_id < num_trajectories; tr_id++)
			{
				int thread_id = omp_get_thread_num();
				one_sub_period_deep(ad, tr_id, part_id, thread_id);
				calc_chars_std(ad, tr_id);
				lambda_lpn_per_periods(ad, tmp, tr_id, 0, num_sub_steps, global_point_id, num_periods);
				evo_chars_std(ad, tr_id, dump_id);
				evo_chars_lpn(ad, tr_id, dump_id);
			}

#pragma omp parallel for
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				if (dump_evo_sep == 1)
				{
					dump_adr_single(ad, tr_id, true);
				}
			}
		}
	}

	delete tmp;
}

void JCSCoreBehaviour::ex_period_obs_deep_cd(AllData * ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

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

			double observable = ed->spec[tr_id].real;

			for (int cd_st_id = 0; cd_st_id < ed->cd_dim; cd_st_id++)
			{
				curr_point_id = global_point_id - cd_st_id;

				if (curr_point_id >= 0 && curr_point_id < ed->cd_num_points)
				{
					ed->cd_rec_data[tr_id][curr_point_id][cd_st_id] = observable;
				}
			}

			one_sub_period_deep(ad, tr_id, part_id, th_id);
			calc_chars_std(ad, tr_id);
		}
	}
}

void JCSCoreBehaviour::ex_period_obs_deep_sigma(AllData * ad, int tr_id, int th_id, int period_id) const
{
	stringstream msg;
	msg << "Error: need to add functionality" << endl;
	Error(msg.str());
}

void JCSCoreBehaviour::rk_period(AllData * ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int num_branches = md->num_ham_qj;

	double jcs_drv_part_1 = double(cp->params.find("jcs_drv_part_1")->second);
	double jcs_drv_part_2 = double(cp->params.find("jcs_drv_part_2")->second);

	double T_1 = jcs_drv_part_1 * md->T;
	double T_2 = jcs_drv_part_2 * md->T;

	int step_id = 0;

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		double time = 0.0;
		double step = 0.0;
		if (part_id == 0)
		{
			time = double(period_id) * (T_1 + T_2);
			step = T_1 / double(cp->rk_ns);
		}
		else if (part_id == 1)
		{
			time = double(period_id) * (T_1 + T_2) + T_1;
			step = T_2 / double(cp->rk_ns);
		}
		else
		{
			stringstream msg;
			msg << "Error: Wrong part_id" << endl;
			Error(msg.str());
		}

		for (int in_step_id = 0; in_step_id < cp->rk_ns; in_step_id++)
		{
			ed->times_all[tr_id] = time + double(in_step_id) * step;
			rk_step_jcs(ad, tr_id, th_id, step);
		}
	}
}

void JCSCoreBehaviour::rk_period_trp_deep(AllData * ad, int tr_id, int th_id, int period_id) const
{
	stringstream msg;
	msg << "Error: need to add functionality" << endl;
	Error(msg.str());
}

void JCSCoreBehaviour::rk_period_obs_deep(AllData * ad, int tr_id, int th_id, int period_id) const
{
	stringstream msg;
	msg << "Error: need to add functionality" << endl;
	Error(msg.str());
}

void JCSCoreBehaviour::rk_period_obs_deep_lpn(AllData * ad, int period_id) const
{
	stringstream msg;
	msg << "Error: need to add functionality" << endl;
	Error(msg.str());
}

void JCSCoreBehaviour::rk_period_obs_deep_cd(AllData * ad, int tr_id, int th_id, int period_id) const
{
	stringstream msg;
	msg << "Error: need to add functionality" << endl;
	Error(msg.str());
}

void JCSCoreBehaviour::rk_period_obs_deep_sigma(AllData * ad, int tr_id, int th_id, int period_id) const
{
	stringstream msg;
	msg << "Error: need to add functionality" << endl;
	Error(msg.str());
}

void JCSCoreBehaviour::calc_chars_std_start(AllData * ad, int tr_id) const
{
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;
	MKL_Complex16 * phi = &(ed->phi_all[tr_id * sys_size]);
	double * adr = &(ed->abs_diag_rho_all[tr_id * sys_size]);

	double norm = norm_square(phi, sys_size);

	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		adr[st_id] = mult_scalar_double(mult_scalar_complex(&phi[st_id], &phi[st_id], 1), 1.0 / norm).real;
	}

	MKL_Complex16 spec = get_spec_jcs(ad, tr_id);
	MKL_Complex16 num_photons = get_num_photons_jcs(ad, tr_id);

	ed->norm[tr_id] = norm;
	ed->spec[tr_id] = spec;
	ed->mean[tr_id] = num_photons.real;

	calc_random_obs(ad, tr_id);
}

void JCSCoreBehaviour::calc_chars_std(AllData * ad, int tr_id) const
{
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;
	MKL_Complex16 * phi = &(ed->phi_all[tr_id * sys_size]);
	double * adr = &(ed->abs_diag_rho_all[tr_id * sys_size]);

	double norm = norm_square(phi, sys_size);

	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		adr[st_id] = mult_scalar_double(mult_scalar_complex(&phi[st_id], &phi[st_id], 1), 1.0 / norm).real;
	}

	MKL_Complex16 spec = get_spec_jcs(ad, tr_id);
	MKL_Complex16 num_photons = get_num_photons_jcs(ad, tr_id);

	ed->norm[tr_id] = norm;
	ed->spec[tr_id] = spec;
	ed->mean[tr_id] = num_photons.real;

	calc_random_obs(ad, tr_id);
}

void JCSCoreBehaviour::calc_chars_lpn_start(AllData * ad, int tr_id, int base_tr_id) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;

	double lambda = 0.0;
	double lambda_now = 0.0;
	MKL_Complex16 spec_lpn = ed->spec[tr_id];
	double mean_lpn = ed->mean[tr_id];

	double delta_s = this->calc_delta_f(ad, tr_id, base_tr_id); // Important! Here we use calc_delta_f not calc_delta_s

	ed->lambda[tr_id] = lambda;
	ed->lambda_now[tr_id] = lambda_now;
	ed->delta_s[tr_id] = delta_s;

	ed->spec_lpn[tr_id] = spec_lpn;
	ed->mean_lpn[tr_id] = mean_lpn;

	calc_random_obs_lpn_start(ad, tr_id);
}

void JCSCoreBehaviour::calc_chars_lpn(AllData * ad, int tr_id, int base_tr_id) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;

	MKL_Complex16 spec_lpn = get_spec_jcs(ad, tr_id);
	double mean_lpn = get_num_photons_jcs(ad, tr_id).real;

	ed->spec_lpn[tr_id] = spec_lpn;
	ed->mean_lpn[tr_id] = mean_lpn;

	calc_random_obs_lpn(ad, tr_id);
}

double JCSCoreBehaviour::calc_T(AllData * ad) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;

	double jcs_drv_part_1 = double(cp->params.find("jcs_drv_part_1")->second);
	double jcs_drv_part_2 = double(cp->params.find("jcs_drv_part_2")->second);
	double T_1 = jcs_drv_part_1 * md->T;
	double T_2 = jcs_drv_part_2 * md->T;
	double T = T_1 + T_2;

	return T;
}

void JCSCoreBehaviour::evo_chars_std(AllData * ad, int tr_id, int dump_id) const
{
	ConfigParam * cp = ad->cp;
	ExpData * ed = ad->ed;

	int dump_num_total = ed->dump_num_total;

	int index = tr_id * dump_num_total + dump_id;

	ed->norm_evo[index] = ed->norm[tr_id];
	ed->spec_evo[index] = ed->spec[tr_id];
	ed->mean_evo[index] = ed->mean[tr_id];

	int num_random_obs = int(cp->params.find("num_random_obs")->second);
	for (int obs_id = 0; obs_id < num_random_obs; obs_id++)
	{
		ed->random_obs_evo[tr_id][obs_id].push_back(ed->random_obs[tr_id][obs_id]);
	}
}

void JCSCoreBehaviour::evo_chars_lpn(AllData * ad, int tr_id, int dump_id) const
{
	ConfigParam * cp = ad->cp;
	ExpData * ed = ad->ed;

	int dump_num_total = ed->dump_num_total;

	int index = tr_id * dump_num_total + dump_id;

	ed->lambda_evo[index] = ed->lambda_now[tr_id];
	ed->spec_lpn_evo[index] = ed->spec_lpn[tr_id];
	ed->mean_lpn_evo[index] = ed->mean_lpn[tr_id];

	int num_random_obs = int(cp->params.find("num_random_obs")->second);
	for (int obs_id = 0; obs_id < num_random_obs; obs_id++)
	{
		ed->random_obs_lpn_evo[tr_id][obs_id].push_back(ed->random_obs_lpn[tr_id][obs_id]);
	}
}

double JCSCoreBehaviour::calc_delta_s(AllData * ad, int tr_id, int base_tr_id) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;

	int lpn_type = int(cp->params.find("lpn_type")->second);
	int num_random_obs = int(ad->cp->params.find("num_random_obs")->second);

	double delta_s = 0.0;
	if (lpn_type == -1)
	{
		MKL_Complex16 base = ed->spec_lpn[base_tr_id];
		MKL_Complex16 var = ed->spec_lpn[tr_id];
		double tmp = pow((base.real - var.real), 2) + pow((base.imag - var.imag), 2);
		delta_s = sqrt(tmp);
	}
	else if (lpn_type == -2)
	{
		double base = ed->mean_lpn[base_tr_id];
		double var = ed->mean_lpn[tr_id];
		delta_s = fabs(var - base) / double(sys_size);
	}
	else if (lpn_type >= 0 && lpn_type < num_random_obs)
	{
		complex<double> base = ed->random_obs_lpn[base_tr_id][lpn_type];
		complex<double> var = ed->random_obs_lpn[tr_id][lpn_type];
		double tmp = pow((base.real() - var.real()), 2) + pow((base.imag() - var.imag()), 2);
		delta_s = sqrt(tmp);
	}
	else
	{
		stringstream msg;
		msg << "Error: Wrong lpn_type: " << lpn_type << endl;
		Error(msg.str());
	}

	return delta_s;
}

double JCSCoreBehaviour::calc_delta_f(AllData * ad, int tr_id, int base_tr_id) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;

	int lpn_type = int(cp->params.find("lpn_type")->second);
	int num_random_obs = int(ad->cp->params.find("num_random_obs")->second);

	double delta_f = 0.0;
	if (lpn_type == -1)
	{
		MKL_Complex16 base = ed->spec[base_tr_id];
		MKL_Complex16 var = ed->spec[tr_id];
		double tmp = pow((base.real - var.real), 2) + pow((base.imag - var.imag), 2);
		delta_f = sqrt(tmp);
	}
	else if (lpn_type == -2)
	{
		double base = ed->mean[base_tr_id];
		double var = ed->mean[tr_id];
		delta_f = fabs(var - base) / double(sys_size);
	}
	else if (lpn_type >= 0 && lpn_type < num_random_obs)
	{
		complex<double> base = ed->random_obs[base_tr_id][lpn_type];
		complex<double> var = ed->random_obs[tr_id][lpn_type];
		double tmp = pow((base.real() - var.real()), 2) + pow((base.imag() - var.imag()), 2);
		delta_f = sqrt(tmp);
	}
	else
	{
		stringstream msg;
		msg << "Error: Wrong lpn_type: " << lpn_type << endl;
		Error(msg.str());
	}

	return delta_f;
}

void JCSCoreBehaviour::calc_ci(AllData * ad, int tr_id) const
{
	calc_ci_double(ad, tr_id);
}

void JCSCoreBehaviour::dump_std(AllData * ad) const
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	ExpData * ed = ad->ed;

	int dump_obs = int(cp->params.find("dump_obs")->second);

	if (dump_obs == 1)
	{
		int num_trajectories = cp->num_trajectories;

		double * norm = ed->norm;
		MKL_Complex16 * spec = ed->spec;

		string fn;

		fn = rp->path + "norm" + cp->fn_suffix;
		save_double_data(fn, norm, num_trajectories, 16, false);

		fn = rp->path + "spec" + cp->fn_suffix;
		save_complex_data(fn, spec, num_trajectories, 16, false);
	}
}

void JCSCoreBehaviour::dump_lpn(AllData * ad) const
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	ExpData * ed = ad->ed;

	int num_trajectories = cp->num_trajectories;

	int dump_obs = int(cp->params.find("dump_obs")->second);

	if (dump_obs == 1)
	{
		double * lambda = ed->lambda_now;
		MKL_Complex16 * spec_lpn = ed->spec_lpn;
		double * mean = ed->mean_lpn;

		string fn;

		fn = rp->path + "lambda" + cp->fn_suffix;
		save_double_data(fn, lambda, num_trajectories, 16, false);

		int* num_renorms = ed->num_renorms;
		fn = rp->path + "num_renorms" + cp->fn_suffix;
		save_int_data(fn, num_renorms, num_trajectories, false);

		fn = rp->path + "spec_lpn" + cp->fn_suffix;
		save_complex_data(fn, spec_lpn, num_trajectories, 16, false);

		fn = rp->path + "mean" + cp->fn_suffix;
		save_double_data(fn, mean, num_trajectories, 16, false);
	}
}

void JCSCoreBehaviour::dump_std_evo(AllData * ad) const
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	ExpData * ed = ad->ed;

	int num_trajectories = cp->num_trajectories;
	int dump_num_total = ed->dump_num_total;

	int dump_obs = int(cp->params.find("dump_obs")->second);
	int jump = int(cp->params.find("jump")->second);

	if (dump_obs == 1)
	{
		int * dump_periods = ed->dump_periods;

		double * norm_evo = ed->norm_evo;
		MKL_Complex16 * spec_evo = ed->spec_evo;
		double * mean_evo = ed->mean_evo;

		string fn;

		fn = rp->path + "norm_evo" + cp->fn_suffix;
		save_2d_inv_double_data(fn, norm_evo, dump_num_total, num_trajectories, 16, false);

		fn = rp->path + "periods" + cp->fn_suffix;
		save_int_data(fn, dump_periods, dump_num_total, false);

		fn = rp->path + "spec_evo" + cp->fn_suffix;
		save_2d_inv_complex_data(fn, spec_evo, dump_num_total, num_trajectories, 16, false);

		fn = rp->path + "mean_evo" + cp->fn_suffix;
		save_2d_inv_double_data(fn, mean_evo, dump_num_total, num_trajectories, 16, false);

		int num_random_obs = int(cp->params.find("num_random_obs")->second);
		int random_obs_seed = int(cp->params.find("random_obs_seed")->second);
		for (int obs_id = 0; obs_id < num_random_obs; obs_id++)
		{
			fn = rp->path + "random_obs_evo_" + to_string(random_obs_seed + obs_id) + cp->fn_suffix;
			std::vector<std::vector<std::complex<double>>> random_obs_evo_to_dump;
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				random_obs_evo_to_dump.push_back(ed->random_obs_evo[tr_id][obs_id]);
			}
			save_2d_vector(fn, random_obs_evo_to_dump, 16);
		}

		if (jump > 0)
		{
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				fn = rp->path + "jump_times_" + to_string(tr_id) + cp->fn_suffix;
				save_double_vector(fn, ed->jump_times[tr_id], 16, false);

				fn = rp->path + "diss_types_" + to_string(tr_id) + cp->fn_suffix;
				save_int_vector(fn, ed->diss_types[tr_id], false);

				if (jump > 1)
				{
					fn = rp->path + "jump_norms_" + to_string(tr_id) + cp->fn_suffix;
					save_double_vector(fn, ed->jump_norms[tr_id], 16, false);

					fn = rp->path + "jump_etas_" + to_string(tr_id) + cp->fn_suffix;
					save_double_vector(fn, ed->jump_etas[tr_id], 16, false);
				}
			}
		}
	}
}

void JCSCoreBehaviour::dump_lpn_evo(AllData * ad) const
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	ExpData * ed = ad->ed;

	int num_trajectories = cp->num_trajectories;
	int dump_num_total = ed->dump_num_total;

	int dump_obs = int(cp->params.find("dump_obs")->second);

	if (dump_obs == 1)
	{
		double * lambda_evo = ed->lambda_evo;
		MKL_Complex16 * spec_lpn_evo = ed->spec_lpn_evo;
		double * mean_lpn_evo = ed->mean_lpn_evo;

		string fn;

		fn = rp->path + "lambda_evo" + cp->fn_suffix;
		save_2d_inv_double_data(fn, lambda_evo, dump_num_total, num_trajectories, 16, false);

		fn = rp->path + "spec_lpn_evo" + cp->fn_suffix;
		save_2d_inv_complex_data(fn, spec_lpn_evo, dump_num_total, num_trajectories, 16, false);

		fn = rp->path + "mean_lpn_evo" + cp->fn_suffix;
		save_2d_inv_double_data(fn, mean_lpn_evo, dump_num_total, num_trajectories, 16, false);

		int num_random_obs = int(cp->params.find("num_random_obs")->second);
		int random_obs_seed = int(cp->params.find("random_obs_seed")->second);
		for (int obs_id = 0; obs_id < num_random_obs; obs_id++)
		{
			fn = rp->path + "random_obs_lpn_evo_" + to_string(random_obs_seed + obs_id) + cp->fn_suffix;
			std::vector<std::vector<std::complex<double>>> random_obs_evo_to_dump;
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				random_obs_evo_to_dump.push_back(ed->random_obs_lpn_evo[tr_id][obs_id]);
			}
			save_2d_vector(fn, random_obs_evo_to_dump, 16);
		}
	}
}










void PSCoreBehaviour::init_splits(AllData * ad) const
{
	RunParam * rp = ad->rp;
	MainData * md = ad->md;

	int num_branches = md->num_ham_qj;
	int num_threads = rp->num_threads;

	int num_total = num_threads * num_branches;

	md->structure = init_split_structure_ps(ad);
	cout << "init_split_structure_ps" << endl;
	md->splits = new Split[num_total];

	for (int b_id = 0; b_id < num_branches; b_id++)
	{
		for (int th_id = 0; th_id < num_threads; th_id++)
		{
			int index = b_id * num_threads + th_id;
			copy_struct_not_member(&(md->structure)[b_id], &(md->splits)[index]);
		}
	}
}

void PSCoreBehaviour::free_splits(AllData * ad) const
{
	RunParam * rp = ad->rp;
	MainData * md = ad->md;

	int num_branches = md->num_ham_qj;
	int num_threads = rp->num_threads;

	for (int b_id = 0; b_id < num_branches; b_id++)
	{
		for (int th_id = 0; th_id < num_threads; th_id++)
		{
			int index = b_id * num_threads + th_id;
			delete_split_struct_not_member(&(md->splits[index]));
		}
	}

	delete(md->splits);
	delete_split_struct(md->structure);
}

void PSCoreBehaviour::init_splits_deep(AllData * ad) const
{
	RunParam * rp = ad->rp;
	MainData * md = ad->md;

	int num_branches = md->num_ham_qj;
	int num_threads = rp->num_threads;

	int num_total = num_threads * num_branches;

	md->structure = init_split_structure_ps_deep(ad);
	md->splits = new Split[num_total];

	for (int b_id = 0; b_id < num_branches; b_id++)
	{
		for (int th_id = 0; th_id < num_threads; th_id++)
		{
			int index = b_id * num_threads + th_id;
			copy_struct_not_member(&(md->structure)[b_id], &(md->splits)[index]);
		}
	}
}

void PSCoreBehaviour::free_splits_deep(AllData * ad) const
{
	RunParam * rp = ad->rp;
	MainData * md = ad->md;

	int num_branches = md->num_ham_qj;
	int num_threads = rp->num_threads;

	for (int b_id = 0; b_id < num_branches; b_id++)
	{
		for (int th_id = 0; th_id < num_threads; th_id++)
		{
			int index = b_id * num_threads + th_id;
			delete_split_struct_not_member(&(md->splits[index]));
		}
	}

	delete(md->splits);
	delete_split_struct(md->structure);
}

void PSCoreBehaviour::ex_period(AllData * ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;

	int num_branches = md->num_ham_qj;

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		one_sub_period_deep(ad, tr_id, part_id, th_id);
	}
}

void PSCoreBehaviour::ex_period_trp_deep(AllData * ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;

	int num_branches = md->num_ham_qj;
	int num_sub_steps = int(cp->params.find("deep_num_steps")->second);

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		for (int sub_step_id = 0; sub_step_id < num_sub_steps; sub_step_id++)
		{
			one_sub_period_deep(ad, tr_id, part_id, th_id);
		}
	}
}

void PSCoreBehaviour::ex_period_obs_deep(AllData * ad, int tr_id, int th_id, int period_id) const
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

			one_sub_period_deep(ad, tr_id, part_id, th_id);
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

void PSCoreBehaviour::ex_period_obs_deep_lpn(AllData * ad, int period_id) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	double ps_drv_part_1 = double(cp->params.find("ps_drv_part_1")->second);
	double ps_drv_part_2 = double(cp->params.find("ps_drv_part_2")->second);
	double T_1 = ps_drv_part_1 * md->T;
	double T_2 = ps_drv_part_2 * md->T;
	double T = T_1 + T_2;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);

	int num_trajectories = cp->num_trajectories;

	int num_branches = md->num_ham_qj;
	int num_sub_steps_per_part = int(cp->params.find("deep_num_steps")->second);
	int num_sub_steps = num_branches * int(cp->params.find("deep_num_steps")->second);

	int dump_point_id = 0;
	int global_point_id = 0;

	int step_id = 0;

	CoreBehavior * tmp = new PSCoreBehaviour;

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		for (int sub_step_id = 0; sub_step_id < num_sub_steps_per_part; sub_step_id++)
		{
			step_id = part_id * num_sub_steps_per_part + sub_step_id;

			global_point_id = period_id * num_sub_steps + dump_point_id;
			dump_point_id++;
			int dump_id = global_point_id + 1;

			ed->curr_time = double(dump_id) / double(num_sub_steps) * T;

			one_sub_period_deep(ad, 0, part_id, 0);
			calc_chars_std(ad, 0);
			calc_chars_lpn(ad, 0, 0);
			evo_chars_std(ad, 0, dump_id);
			evo_chars_lpn(ad, 0, dump_id);

#pragma omp parallel for
			for (int tr_id = 1; tr_id < num_trajectories; tr_id++)
			{
				int thread_id = omp_get_thread_num();
				one_sub_period_deep(ad, tr_id, part_id, thread_id);
				calc_chars_std(ad, tr_id);
				lambda_lpn(ad, tmp, tr_id, 0);
				evo_chars_std(ad, tr_id, dump_id);
				evo_chars_lpn(ad, tr_id, dump_id);
			}

#pragma omp parallel for
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				if (dump_evo_sep == 1)
				{
					dump_adr_single(ad, tr_id, true);
				}
			}
		}
	}

	delete tmp;
}

void PSCoreBehaviour::ex_period_obs_deep_mult_lpn(AllData * ad, int period_id) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	double ps_drv_part_1 = double(cp->params.find("ps_drv_part_1")->second);
	double ps_drv_part_2 = double(cp->params.find("ps_drv_part_2")->second);
	double T_1 = ps_drv_part_1 * md->T;
	double T_2 = ps_drv_part_2 * md->T;
	double T = T_1 + T_2;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);
	int lambda_per_periods = int(cp->params.find("lambda_per_periods")->second);

	int num_trajectories = cp->num_trajectories;

	int num_branches = md->num_ham_qj;
	int num_sub_steps_per_part = int(cp->params.find("deep_num_steps")->second);
	int num_sub_steps = num_branches * int(cp->params.find("deep_num_steps")->second);

	int dump_point_id = 0;
	int global_point_id = 0;

	int step_id = 0;

	CoreBehavior * tmp = new PSCoreBehaviour;

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		for (int sub_step_id = 0; sub_step_id < num_sub_steps_per_part; sub_step_id++)
		{
			step_id = part_id * num_sub_steps_per_part + sub_step_id;

			global_point_id = period_id * num_sub_steps + dump_point_id;
			dump_point_id++;
			int dump_id = global_point_id + 1;

			ed->curr_time = double(dump_id) / double(num_sub_steps) * T;

#pragma omp parallel for
			for (int tr_id = 0; tr_id < num_trajectories / 2; tr_id++)
			{
				int thread_id = omp_get_thread_num();
				one_sub_period_deep(ad, tr_id, part_id, thread_id);
				calc_chars_std(ad, tr_id);
				calc_chars_lpn(ad, tr_id, tr_id);
				evo_chars_std(ad, tr_id, dump_id);
				evo_chars_lpn(ad, tr_id, dump_id);
			}

#pragma omp parallel for
			for (int tr_id = num_trajectories / 2; tr_id < num_trajectories; tr_id++)
			{
				int thread_id = omp_get_thread_num();
				one_sub_period_deep(ad, tr_id, part_id, thread_id);
				calc_chars_std(ad, tr_id);
				if (lambda_per_periods == 0)
				{
					lambda_lpn(ad, tmp, tr_id, tr_id - num_trajectories / 2);
				}
				if (lambda_per_periods == 1)
				{
					if ((part_id == num_branches - 1) && (sub_step_id == num_sub_steps_per_part - 1))
					{
						lambda_lpn_now(ad, tmp, tr_id, tr_id - num_trajectories / 2);
					}
				}
				evo_chars_std(ad, tr_id, dump_id);
				evo_chars_lpn(ad, tr_id, dump_id);
			}

#pragma omp parallel for
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				if (dump_evo_sep == 1)
				{
					dump_adr_single(ad, tr_id, true);
				}
			}
		}
	}

	delete tmp;
}

void PSCoreBehaviour::ex_period_obs_deep_lpn_per_period(struct AllData* ad, int period_id, int num_periods) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	double ps_drv_part_1 = double(cp->params.find("ps_drv_part_1")->second);
	double ps_drv_part_2 = double(cp->params.find("ps_drv_part_2")->second);
	double T_1 = ps_drv_part_1 * md->T;
	double T_2 = ps_drv_part_2 * md->T;
	double T = T_1 + T_2;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);

	int num_trajectories = cp->num_trajectories;

	int num_branches = md->num_ham_qj;
	int num_sub_steps_per_part = int(cp->params.find("deep_num_steps")->second);
	int num_sub_steps = num_branches * int(cp->params.find("deep_num_steps")->second);

	int dump_point_id = 0;
	int global_point_id = 0;

	int step_id = 0;

	CoreBehavior * tmp = new PSCoreBehaviour;

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		for (int sub_step_id = 0; sub_step_id < num_sub_steps_per_part; sub_step_id++)
		{
			step_id = part_id * num_sub_steps_per_part + sub_step_id;

			global_point_id = period_id * num_sub_steps + dump_point_id;
			dump_point_id++;
			int dump_id = global_point_id + 1;

			ed->curr_time = double(dump_id) / double(num_sub_steps) * T;

			one_sub_period_deep(ad, 0, part_id, 0);
			calc_chars_std(ad, 0);
			calc_chars_lpn(ad, 0, 0);
			evo_chars_std(ad, 0, dump_id);
			evo_chars_lpn(ad, 0, dump_id);

#pragma omp parallel for
			for (int tr_id = 1; tr_id < num_trajectories; tr_id++)
			{
				int thread_id = omp_get_thread_num();
				one_sub_period_deep(ad, tr_id, part_id, thread_id);
				calc_chars_std(ad, tr_id);
				lambda_lpn_per_periods(ad, tmp, tr_id, 0, num_sub_steps, global_point_id, num_periods);
				evo_chars_std(ad, tr_id, dump_id);
				evo_chars_lpn(ad, tr_id, dump_id);
			}

#pragma omp parallel for
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				if (dump_evo_sep == 1)
				{
					dump_adr_single(ad, tr_id, true);
				}
			}
		}
	}

	delete tmp;
}

void PSCoreBehaviour::ex_period_obs_deep_cd(AllData * ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

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

			double observable = ed->spec[tr_id].real;

			for (int cd_st_id = 0; cd_st_id < ed->cd_dim; cd_st_id++)
			{
				curr_point_id = global_point_id - cd_st_id;

				if (curr_point_id >= 0 && curr_point_id < ed->cd_num_points)
				{
					ed->cd_rec_data[tr_id][curr_point_id][cd_st_id] = observable;
				}
			}

			one_sub_period_deep(ad, tr_id, part_id, th_id);
			calc_chars_std(ad, tr_id);
		}
	}
}

void PSCoreBehaviour::ex_period_obs_deep_sigma(AllData * ad, int tr_id, int th_id, int period_id) const
{
	stringstream msg;
	msg << "Error: need to add functionality" << endl;
	Error(msg.str());
}

void PSCoreBehaviour::rk_period(AllData * ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int num_branches = md->num_ham_qj;

	double ps_drv_part_1 = double(cp->params.find("ps_drv_part_1")->second);
	double ps_drv_part_2 = double(cp->params.find("ps_drv_part_2")->second);

	double T_1 = ps_drv_part_1 * md->T;
	double T_2 = ps_drv_part_2 * md->T;

	int step_id = 0;

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		double time = 0.0;
		double step = 0.0;
		if (part_id == 0)
		{
			time = double(period_id) * (T_1 + T_2);
			step = T_1 / double(cp->rk_ns);
		}
		else if (part_id == 1)
		{
			time = double(period_id) * (T_1 + T_2) + T_1;
			step = T_2 / double(cp->rk_ns);
		}
		else
		{
			stringstream msg;
			msg << "Error: Wrong part_id" << endl;
			Error(msg.str());
		}

		for (int in_step_id = 0; in_step_id < cp->rk_ns; in_step_id++)
		{
			ed->times_all[tr_id] = time + double(in_step_id) * step;
			rk_step_ps(ad, tr_id, th_id, step);
		}
	}
}

void PSCoreBehaviour::rk_period_trp_deep(AllData * ad, int tr_id, int th_id, int period_id) const
{
	stringstream msg;
	msg << "Error: need to add functionality" << endl;
	Error(msg.str());
}

void PSCoreBehaviour::rk_period_obs_deep(AllData * ad, int tr_id, int th_id, int period_id) const
{
	stringstream msg;
	msg << "Error: need to add functionality" << endl;
	Error(msg.str());
}

void PSCoreBehaviour::rk_period_obs_deep_lpn(AllData * ad, int period_id) const
{
	stringstream msg;
	msg << "Error: need to add functionality" << endl;
	Error(msg.str());
}

void PSCoreBehaviour::rk_period_obs_deep_cd(AllData * ad, int tr_id, int th_id, int period_id) const
{
	stringstream msg;
	msg << "Error: need to add functionality" << endl;
	Error(msg.str());
}

void PSCoreBehaviour::rk_period_obs_deep_sigma(AllData * ad, int tr_id, int th_id, int period_id) const
{
	stringstream msg;
	msg << "Error: need to add functionality" << endl;
	Error(msg.str());
}

void PSCoreBehaviour::calc_chars_std_start(AllData * ad, int tr_id) const
{
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;
	MKL_Complex16 * phi = &(ed->phi_all[tr_id * sys_size]);
	double * adr = &(ed->abs_diag_rho_all[tr_id * sys_size]);

	double norm = norm_square(phi, sys_size);

	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		adr[st_id] = mult_scalar_double(mult_scalar_complex(&phi[st_id], &phi[st_id], 1), 1.0 / norm).real;
	}

	MKL_Complex16 spec = get_spec_ps(ad, tr_id);
	MKL_Complex16 spec_2 = get_spec_2_ps(ad, tr_id);
	MKL_Complex16 spec_3 = get_spec_3_ps(ad, tr_id);
	MKL_Complex16 num_photons = get_num_photons_ps(ad, tr_id);

	ed->norm[tr_id] = norm;
	ed->spec[tr_id] = spec;
	ed->spec_2[tr_id] = spec_2;
	ed->spec_3[tr_id] = spec_3;
	ed->mean[tr_id] = num_photons.real;

	calc_random_obs(ad, tr_id);
}

void PSCoreBehaviour::calc_chars_std(AllData * ad, int tr_id) const
{
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;
	MKL_Complex16 * phi = &(ed->phi_all[tr_id * sys_size]);
	double * adr = &(ed->abs_diag_rho_all[tr_id * sys_size]);

	double norm = norm_square(phi, sys_size);

	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		adr[st_id] = mult_scalar_double(mult_scalar_complex(&phi[st_id], &phi[st_id], 1), 1.0 / norm).real;
	}

	MKL_Complex16 spec = get_spec_ps(ad, tr_id);
	MKL_Complex16 spec_2 = get_spec_2_ps(ad, tr_id);
	MKL_Complex16 spec_3 = get_spec_3_ps(ad, tr_id);
	MKL_Complex16 num_photons = get_num_photons_ps(ad, tr_id);

	ed->norm[tr_id] = norm;
	ed->spec[tr_id] = spec;
	ed->spec_2[tr_id] = spec_2;
	ed->spec_3[tr_id] = spec_3;
	ed->mean[tr_id] = num_photons.real;

	calc_random_obs(ad, tr_id);
}

void PSCoreBehaviour::calc_chars_lpn_start(AllData * ad, int tr_id, int base_tr_id) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;

	double lambda = 0.0;
	double lambda_now = 0.0;
	MKL_Complex16 spec_lpn = ed->spec[tr_id];
	MKL_Complex16 spec_2_lpn = ed->spec_2[tr_id];
	MKL_Complex16 spec_3_lpn = ed->spec_3[tr_id];
	double mean_lpn = ed->mean[tr_id];

	double delta_s = this->calc_delta_f(ad, tr_id, base_tr_id); // Important! Here we use calc_delta_f not calc_delta_s

	ed->lambda[tr_id] = lambda;
	ed->lambda_now[tr_id] = lambda_now;
	ed->delta_s[tr_id] = delta_s;

	ed->spec_lpn[tr_id] = spec_lpn;
	ed->spec_2_lpn[tr_id] = spec_2_lpn;
	ed->spec_3_lpn[tr_id] = spec_3_lpn;
	ed->mean_lpn[tr_id] = mean_lpn;

	calc_random_obs_lpn_start(ad, tr_id);
}

void PSCoreBehaviour::calc_chars_lpn(AllData * ad, int tr_id, int base_tr_id) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;

	MKL_Complex16 spec_lpn = get_spec_ps(ad, tr_id);
	MKL_Complex16 spec_2_lpn = get_spec_2_ps(ad, tr_id);
	MKL_Complex16 spec_3_lpn = get_spec_3_ps(ad, tr_id);
	double mean_lpn = get_num_photons_ps(ad, tr_id).real;

	ed->spec_lpn[tr_id] = spec_lpn;
	ed->spec_2_lpn[tr_id] = spec_2_lpn;
	ed->spec_3_lpn[tr_id] = spec_3_lpn;
	ed->mean_lpn[tr_id] = mean_lpn;

	calc_random_obs_lpn(ad, tr_id);
}

double PSCoreBehaviour::calc_T(AllData * ad) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;

	double ps_drv_part_1 = double(cp->params.find("ps_drv_part_1")->second);
	double ps_drv_part_2 = double(cp->params.find("ps_drv_part_2")->second);
	double T_1 = ps_drv_part_1 * md->T;
	double T_2 = ps_drv_part_2 * md->T;
	double T = T_1 + T_2;

	return T;
}

void PSCoreBehaviour::evo_chars_std(AllData * ad, int tr_id, int dump_id) const
{
	ConfigParam * cp = ad->cp;
	ExpData * ed = ad->ed;

	int dump_num_total = ed->dump_num_total;

	int index = tr_id * dump_num_total + dump_id;

	ed->norm_evo[index] = ed->norm[tr_id];
	ed->spec_evo[index] = ed->spec[tr_id];
	ed->spec_2_evo[index] = ed->spec_2[tr_id];
	ed->spec_3_evo[index] = ed->spec_3[tr_id];
	ed->mean_evo[index] = ed->mean[tr_id];

	int num_random_obs = int(cp->params.find("num_random_obs")->second);
	for (int obs_id = 0; obs_id < num_random_obs; obs_id++)
	{
		ed->random_obs_evo[tr_id][obs_id].push_back(ed->random_obs[tr_id][obs_id]);
	}
}

void PSCoreBehaviour::evo_chars_lpn(AllData * ad, int tr_id, int dump_id) const
{
	ConfigParam * cp = ad->cp;
	ExpData * ed = ad->ed;

	int dump_num_total = ed->dump_num_total;

	int index = tr_id * dump_num_total + dump_id;

	ed->lambda_evo[index] = ed->lambda_now[tr_id];
	ed->spec_lpn_evo[index] = ed->spec_lpn[tr_id];
	ed->spec_2_lpn_evo[index] = ed->spec_2_lpn[tr_id];
	ed->spec_3_lpn_evo[index] = ed->spec_3_lpn[tr_id];
	ed->mean_lpn_evo[index] = ed->mean_lpn[tr_id];

	int num_random_obs = int(cp->params.find("num_random_obs")->second);
	for (int obs_id = 0; obs_id < num_random_obs; obs_id++)
	{
		ed->random_obs_lpn_evo[tr_id][obs_id].push_back(ed->random_obs_lpn[tr_id][obs_id]);
	}
}

double PSCoreBehaviour::calc_delta_s(AllData * ad, int tr_id, int base_tr_id) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;

	int lpn_type = int(cp->params.find("lpn_type")->second);
	int num_random_obs = int(ad->cp->params.find("num_random_obs")->second);

	double delta_s = 0.0;
	if (lpn_type == -1)
	{
		MKL_Complex16 base = ed->spec_lpn[base_tr_id];
		MKL_Complex16 var = ed->spec_lpn[tr_id];
		double tmp = pow((base.real - var.real), 2) + pow((base.imag - var.imag), 2);
		delta_s = sqrt(tmp);
	}
	else if (lpn_type == -2)
	{
		MKL_Complex16 base = ed->spec_2_lpn[base_tr_id];
		MKL_Complex16 var = ed->spec_2_lpn[tr_id];
		double tmp = pow((base.real - var.real), 2) + pow((base.imag - var.imag), 2);
		delta_s = sqrt(tmp);
	}
	else if (lpn_type == -3)
	{
		MKL_Complex16 base = ed->spec_3_lpn[base_tr_id];
		MKL_Complex16 var = ed->spec_3_lpn[tr_id];
		double tmp = pow((base.real - var.real), 2) + pow((base.imag - var.imag), 2);
		delta_s = sqrt(tmp);
	}
	else if (lpn_type == -4)
	{
		double base = ed->mean_lpn[base_tr_id];
		double var = ed->mean_lpn[tr_id];
		delta_s = fabs(var - base) / double(sys_size);
	}
	else if (lpn_type >= 0 && lpn_type < num_random_obs)
	{
		complex<double> base = ed->random_obs_lpn[base_tr_id][lpn_type];
		complex<double> var = ed->random_obs_lpn[tr_id][lpn_type];
		double tmp = pow((base.real() - var.real()), 2) + pow((base.imag() - var.imag()), 2);
		delta_s = sqrt(tmp);
	}
	else
	{
		stringstream msg;
		msg << "Error: Wrong lpn_type: " << lpn_type << endl;
		Error(msg.str());
	}

	return delta_s;
}

double PSCoreBehaviour::calc_delta_f(AllData * ad, int tr_id, int base_tr_id) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;

	int lpn_type = int(cp->params.find("lpn_type")->second);
	int num_random_obs = int(ad->cp->params.find("num_random_obs")->second);

	double delta_f = 0.0;
	if (lpn_type == -1)
	{
		MKL_Complex16 base = ed->spec[base_tr_id];
		MKL_Complex16 var = ed->spec[tr_id];
		double tmp = pow((base.real - var.real), 2) + pow((base.imag - var.imag), 2);
		delta_f = sqrt(tmp);
	}
	else if (lpn_type == -2)
	{
		MKL_Complex16 base = ed->spec_2[base_tr_id];
		MKL_Complex16 var = ed->spec_2[tr_id];
		double tmp = pow((base.real - var.real), 2) + pow((base.imag - var.imag), 2);
		delta_f = sqrt(tmp);
	}
	else if (lpn_type == -3)
	{
		MKL_Complex16 base = ed->spec_3[base_tr_id];
		MKL_Complex16 var = ed->spec_3[tr_id];
		double tmp = pow((base.real - var.real), 2) + pow((base.imag - var.imag), 2);
		delta_f = sqrt(tmp);
	}
	else if (lpn_type == -4)
	{
		double base = ed->mean[base_tr_id];
		double var = ed->mean[tr_id];
		delta_f = fabs(var - base) / double(sys_size);
	}
	else if (lpn_type >= 0 && lpn_type < num_random_obs)
	{
		complex<double> base = ed->random_obs[base_tr_id][lpn_type];
		complex<double> var = ed->random_obs[tr_id][lpn_type];
		double tmp = pow((base.real() - var.real()), 2) + pow((base.imag() - var.imag()), 2);
		delta_f = sqrt(tmp);
	}
	else
	{
		stringstream msg;
		msg << "Error: Wrong lpn_type: " << lpn_type << endl;
		Error(msg.str());
	}

	return delta_f;
}

void PSCoreBehaviour::calc_ci(AllData * ad, int tr_id) const
{
	calc_ci_double(ad, tr_id);
}

void PSCoreBehaviour::dump_std(AllData * ad) const
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	ExpData * ed = ad->ed;

	int dump_obs = int(cp->params.find("dump_obs")->second);

	if (dump_obs == 1)
	{
		int num_trajectories = cp->num_trajectories;

		double * norm = ed->norm;
		double * mean = ed->mean;
		MKL_Complex16 * spec = ed->spec;
		MKL_Complex16 * spec_2 = ed->spec_2;
		MKL_Complex16 * spec_3 = ed->spec_3;

		string fn;

		fn = rp->path + "norm" + cp->fn_suffix;
		save_double_data(fn, norm, num_trajectories, 16, false);

		fn = rp->path + "spec" + cp->fn_suffix;
		save_complex_data(fn, spec, num_trajectories, 16, false);

		fn = rp->path + "spec_2" + cp->fn_suffix;
		save_complex_data(fn, spec_2, num_trajectories, 16, false);

		fn = rp->path + "spec_3" + cp->fn_suffix;
		save_complex_data(fn, spec_3, num_trajectories, 16, false);

		fn = rp->path + "mean" + cp->fn_suffix;
		save_double_data(fn, mean, num_trajectories, 16, false);
	}
}

void PSCoreBehaviour::dump_lpn(AllData * ad) const
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	ExpData * ed = ad->ed;

	int num_trajectories = cp->num_trajectories;

	int dump_obs = int(cp->params.find("dump_obs")->second);

	if (dump_obs == 1)
	{
		double * lambda = ed->lambda_now;
		MKL_Complex16 * spec_lpn = ed->spec_lpn;
		MKL_Complex16 * spec_2_lpn = ed->spec_2_lpn;
		MKL_Complex16 * spec_3_lpn = ed->spec_3_lpn;
		double * mean = ed->mean_lpn;

		string fn;

		fn = rp->path + "lambda" + cp->fn_suffix;
		save_double_data(fn, lambda, num_trajectories, 16, false);

		int* num_renorms = ed->num_renorms;
		fn = rp->path + "num_renorms" + cp->fn_suffix;
		save_int_data(fn, num_renorms, num_trajectories, false);

		fn = rp->path + "spec_lpn" + cp->fn_suffix;
		save_complex_data(fn, spec_lpn, num_trajectories, 16, false);

		fn = rp->path + "spec_2_lpn" + cp->fn_suffix;
		save_complex_data(fn, spec_2_lpn, num_trajectories, 16, false);

		fn = rp->path + "spec_3_lpn" + cp->fn_suffix;
		save_complex_data(fn, spec_3_lpn, num_trajectories, 16, false);

		fn = rp->path + "mean_lpn" + cp->fn_suffix;
		save_double_data(fn, mean, num_trajectories, 16, false);

		int save_lambdas = int(ad->cp->params.find("save_lambdas")->second);
		if (save_lambdas > 0)
		{
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				fn = rp->path + "delta_f_by_delta_s_" + to_string(tr_id) + cp->fn_suffix;
				save_double_vector(fn, ed->lambdas[tr_id], 16, false);

				fn = rp->path + "deltas_s_" + to_string(tr_id) + cp->fn_suffix;
				save_double_vector(fn, ed->deltas_s[tr_id], 16, false);

				fn = rp->path + "deltas_f_" + to_string(tr_id) + cp->fn_suffix;
				save_double_vector(fn, ed->deltas_f[tr_id], 16, false);
			}
		}
	}
}

void PSCoreBehaviour::dump_std_evo(AllData * ad) const
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	ExpData * ed = ad->ed;

	int num_trajectories = cp->num_trajectories;
	int dump_num_total = ed->dump_num_total;

	int dump_obs = int(cp->params.find("dump_obs")->second);
	int jump = int(cp->params.find("jump")->second);

	if (dump_obs == 1)
	{
		int * dump_periods = ed->dump_periods;

		double * norm_evo = ed->norm_evo;
		MKL_Complex16 * spec_evo = ed->spec_evo;
		MKL_Complex16 * spec_2_evo = ed->spec_2_evo;
		MKL_Complex16 * spec_3_evo = ed->spec_3_evo;
		double * mean_evo = ed->mean_evo;

		string fn;

		fn = rp->path + "norm_evo" + cp->fn_suffix;
		save_2d_inv_double_data(fn, norm_evo, dump_num_total, num_trajectories, 16, false);

		fn = rp->path + "periods" + cp->fn_suffix;
		save_int_data(fn, dump_periods, dump_num_total, false);

		fn = rp->path + "spec_evo" + cp->fn_suffix;
		save_2d_inv_complex_data(fn, spec_evo, dump_num_total, num_trajectories, 16, false);

		fn = rp->path + "spec_2_evo" + cp->fn_suffix;
		save_2d_inv_complex_data(fn, spec_2_evo, dump_num_total, num_trajectories, 16, false);

		fn = rp->path + "spec_3_evo" + cp->fn_suffix;
		save_2d_inv_complex_data(fn, spec_3_evo, dump_num_total, num_trajectories, 16, false);

		fn = rp->path + "mean_evo" + cp->fn_suffix;
		save_2d_inv_double_data(fn, mean_evo, dump_num_total, num_trajectories, 16, false);

		int num_random_obs = int(cp->params.find("num_random_obs")->second);
		int random_obs_seed = int(cp->params.find("random_obs_seed")->second);
		for (int obs_id = 0; obs_id < num_random_obs; obs_id++)
		{
			fn = rp->path + "random_obs_evo_" + to_string(random_obs_seed + obs_id) + cp->fn_suffix;
			std::vector<std::vector<std::complex<double>>> random_obs_evo_to_dump;
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				random_obs_evo_to_dump.push_back(ed->random_obs_evo[tr_id][obs_id]);
			}
			save_2d_vector(fn, random_obs_evo_to_dump, 16);
		}

		if (jump > 0)
		{
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				fn = rp->path + "jump_times_" + to_string(tr_id) + cp->fn_suffix;
				save_double_vector(fn, ed->jump_times[tr_id], 16, false);

				fn = rp->path + "diss_types_" + to_string(tr_id) + cp->fn_suffix;
				save_int_vector(fn, ed->diss_types[tr_id], false);

				if (jump > 1)
				{
					fn = rp->path + "jump_norms_" + to_string(tr_id) + cp->fn_suffix;
					save_double_vector(fn, ed->jump_norms[tr_id], 16, false);

					fn = rp->path + "jump_etas_" + to_string(tr_id) + cp->fn_suffix;
					save_double_vector(fn, ed->jump_etas[tr_id], 16, false);
				}
			}
		}
	}
}

void PSCoreBehaviour::dump_lpn_evo(AllData * ad) const
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	ExpData * ed = ad->ed;

	int num_trajectories = cp->num_trajectories;
	int dump_num_total = ed->dump_num_total;

	int dump_obs = int(cp->params.find("dump_obs")->second);

	if (dump_obs == 1)
	{
		double * lambda_evo = ed->lambda_evo;
		MKL_Complex16 * spec_lpn_evo = ed->spec_lpn_evo;
		MKL_Complex16 * spec_2_lpn_evo = ed->spec_2_lpn_evo;
		MKL_Complex16 * spec_3_lpn_evo = ed->spec_3_lpn_evo;
		double * mean_lpn_evo = ed->mean_lpn_evo;

		string fn;

		fn = rp->path + "lambda_evo" + cp->fn_suffix;
		save_2d_inv_double_data(fn, lambda_evo, dump_num_total, num_trajectories, 16, false);

		fn = rp->path + "spec_lpn_evo" + cp->fn_suffix;
		save_2d_inv_complex_data(fn, spec_lpn_evo, dump_num_total, num_trajectories, 16, false);

		fn = rp->path + "spec_2_lpn_evo" + cp->fn_suffix;
		save_2d_inv_complex_data(fn, spec_2_lpn_evo, dump_num_total, num_trajectories, 16, false);

		fn = rp->path + "spec_3_lpn_evo" + cp->fn_suffix;
		save_2d_inv_complex_data(fn, spec_3_lpn_evo, dump_num_total, num_trajectories, 16, false);

		fn = rp->path + "mean_lpn_evo" + cp->fn_suffix;
		save_2d_inv_double_data(fn, mean_lpn_evo, dump_num_total, num_trajectories, 16, false);

		int num_random_obs = int(cp->params.find("num_random_obs")->second);
		int random_obs_seed = int(cp->params.find("random_obs_seed")->second);
		for (int obs_id = 0; obs_id < num_random_obs; obs_id++)
		{
			fn = rp->path + "random_obs_lpn_evo_" + to_string(random_obs_seed + obs_id) + cp->fn_suffix;
			std::vector<std::vector<std::complex<double>>> random_obs_evo_to_dump;
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				random_obs_evo_to_dump.push_back(ed->random_obs_lpn_evo[tr_id][obs_id]);
			}
			save_2d_vector(fn, random_obs_evo_to_dump, 16);
		}
	}
}










void MBLCoreBehaviour::init_splits(AllData * ad) const
{
	RunParam * rp = ad->rp;
	MainData * md = ad->md;

	int num_branches = md->num_ham_qj;
	int num_threads = rp->num_threads;

	int num_total = num_threads * num_branches;

	md->structure = init_split_structure_mbl(ad);
	cout << "init_split_structure_mbl" << endl;
	md->splits = new Split[num_total];

	for (int b_id = 0; b_id < num_branches; b_id++)
	{
		for (int th_id = 0; th_id < num_threads; th_id++)
		{
			int index = b_id * num_threads + th_id;
			copy_struct_not_member(&(md->structure)[b_id], &(md->splits)[index]);
		}
	}
}

void MBLCoreBehaviour::free_splits(AllData * ad) const
{
	RunParam * rp = ad->rp;
	MainData * md = ad->md;

	int num_branches = md->num_ham_qj;
	int num_threads = rp->num_threads;

	for (int b_id = 0; b_id < num_branches; b_id++)
	{
		for (int th_id = 0; th_id < num_threads; th_id++)
		{
			int index = b_id * num_threads + th_id;
			delete_split_struct_not_member(&(md->splits[index]));
		}
	}

	delete(md->splits);
	delete_split_struct(md->structure);
}

void MBLCoreBehaviour::init_splits_deep(AllData * ad) const
{
	RunParam * rp = ad->rp;
	MainData * md = ad->md;

	int num_branches = md->num_ham_qj;
	int num_threads = rp->num_threads;

	int num_total = num_threads * num_branches;

	md->structure = init_split_structure_mbl_deep(ad);
	md->splits = new Split[num_total];

	for (int b_id = 0; b_id < num_branches; b_id++)
	{
		for (int th_id = 0; th_id < num_threads; th_id++)
		{
			int index = b_id * num_threads + th_id;
			copy_struct_not_member(&(md->structure)[b_id], &(md->splits)[index]);
		}
	}
}

void MBLCoreBehaviour::free_splits_deep(AllData * ad) const
{
	RunParam * rp = ad->rp;
	MainData * md = ad->md;

	int num_branches = md->num_ham_qj;
	int num_threads = rp->num_threads;

	for (int b_id = 0; b_id < num_branches; b_id++)
	{
		for (int th_id = 0; th_id < num_threads; th_id++)
		{
			int index = b_id * num_threads + th_id;
			delete_split_struct_not_member(&(md->splits[index]));
		}
	}

	delete(md->splits);
	delete_split_struct(md->structure);
}

void MBLCoreBehaviour::ex_period(AllData * ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;

	int num_branches = md->num_ham_qj;

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		one_sub_period_deep(ad, tr_id, part_id, th_id);
	}
}

void MBLCoreBehaviour::ex_period_trp_deep(AllData * ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;

	int num_branches = md->num_ham_qj;
	int num_sub_steps = int(cp->params.find("deep_num_steps")->second);

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		for (int sub_step_id = 0; sub_step_id < num_sub_steps; sub_step_id++)
		{
			one_sub_period_deep(ad, tr_id, part_id, th_id);
		}
	}
}

void MBLCoreBehaviour::ex_period_obs_deep(AllData * ad, int tr_id, int th_id, int period_id) const
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

			one_sub_period_deep(ad, tr_id, part_id, th_id);
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

void MBLCoreBehaviour::ex_period_obs_deep_lpn(AllData * ad, int period_id) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	double T = md->T;
	
	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);

	int num_trajectories = cp->num_trajectories;

	int num_branches = md->num_ham_qj;
	int num_sub_steps_per_part = int(cp->params.find("deep_num_steps")->second);
	int num_sub_steps = num_branches * int(cp->params.find("deep_num_steps")->second);

	int dump_point_id = 0;
	int global_point_id = 0;

	int step_id = 0;

	CoreBehavior * tmp = new MBLCoreBehaviour;

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		for (int sub_step_id = 0; sub_step_id < num_sub_steps_per_part; sub_step_id++)
		{
			step_id = part_id * num_sub_steps_per_part + sub_step_id;

			global_point_id = period_id * num_sub_steps + dump_point_id;
			dump_point_id++;
			int dump_id = global_point_id + 1;

			ed->curr_time = double(dump_id) / double(num_sub_steps) * T;

			one_sub_period_deep(ad, 0, part_id, 0);
			calc_chars_std(ad, 0);
			calc_chars_lpn(ad, 0, 0);
			evo_chars_std(ad, 0, dump_id);
			evo_chars_lpn(ad, 0, dump_id);

#pragma omp parallel for
			for (int tr_id = 1; tr_id < num_trajectories; tr_id++)
			{
				int thread_id = omp_get_thread_num();
				one_sub_period_deep(ad, tr_id, part_id, thread_id);
				calc_chars_std(ad, tr_id);
				lambda_lpn(ad, tmp, tr_id, 0);
				evo_chars_std(ad, tr_id, dump_id);
				evo_chars_lpn(ad, tr_id, dump_id);
			}

#pragma omp parallel for
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				if (dump_evo_sep == 1)
				{
					dump_adr_single(ad, tr_id, true);
				}
			}
		}
	}

	delete tmp;
}

void MBLCoreBehaviour::ex_period_obs_deep_mult_lpn(AllData * ad, int period_id) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	double T = md->T;

	int lambda_per_periods = int(cp->params.find("lambda_per_periods")->second);

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);

	int num_trajectories = cp->num_trajectories;

	int num_branches = md->num_ham_qj;
	int num_sub_steps_per_part = int(cp->params.find("deep_num_steps")->second);
	int num_sub_steps = num_branches * int(cp->params.find("deep_num_steps")->second);

	int dump_point_id = 0;
	int global_point_id = 0;

	int step_id = 0;

	CoreBehavior * tmp = new MBLCoreBehaviour;

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		for (int sub_step_id = 0; sub_step_id < num_sub_steps_per_part; sub_step_id++)
		{
			step_id = part_id * num_sub_steps_per_part + sub_step_id;

			global_point_id = period_id * num_sub_steps + dump_point_id;
			dump_point_id++;
			int dump_id = global_point_id + 1;

			ed->curr_time = double(dump_id) / double(num_sub_steps) * T;

#pragma omp parallel for
			for (int tr_id = 0; tr_id < num_trajectories / 2; tr_id++)
			{
				int thread_id = omp_get_thread_num();
				one_sub_period_deep(ad, tr_id, part_id, thread_id);
				calc_chars_std(ad, tr_id);
				calc_chars_lpn(ad, tr_id, tr_id);
				evo_chars_std(ad, tr_id, dump_id);
				evo_chars_lpn(ad, tr_id, dump_id);
			}

#pragma omp parallel for
			for (int tr_id = num_trajectories / 2; tr_id < num_trajectories; tr_id++)
			{
				int thread_id = omp_get_thread_num();
				one_sub_period_deep(ad, tr_id, part_id, thread_id);
				calc_chars_std(ad, tr_id);
				if (lambda_per_periods == 0)
				{
					lambda_lpn(ad, tmp, tr_id, tr_id - num_trajectories / 2);
				}
				if (lambda_per_periods == 1)
				{
					if ((part_id == num_branches - 1) && (sub_step_id == num_sub_steps_per_part - 1))
					{
						lambda_lpn_now(ad, tmp, tr_id, tr_id - num_trajectories / 2);
					}
					else
					{
						lambda_lpn_now_wo_renorm_just_to_save(ad, tmp, tr_id, tr_id - num_trajectories / 2);
					}
				}
				evo_chars_std(ad, tr_id, dump_id);
				evo_chars_lpn(ad, tr_id, dump_id);
			}

#pragma omp parallel for
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				if (dump_evo_sep == 1)
				{
					dump_adr_single(ad, tr_id, true);
				}
			}
		}
	}

	delete tmp;
}

void MBLCoreBehaviour::ex_period_obs_deep_lpn_per_period(struct AllData* ad, int period_id, int num_periods) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	double T = md->T;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);

	int num_trajectories = cp->num_trajectories;

	int num_branches = md->num_ham_qj;
	int num_sub_steps_per_part = int(cp->params.find("deep_num_steps")->second);
	int num_sub_steps = num_branches * int(cp->params.find("deep_num_steps")->second);

	int dump_point_id = 0;
	int global_point_id = 0;

	int step_id = 0;

	CoreBehavior * tmp = new MBLCoreBehaviour;

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		for (int sub_step_id = 0; sub_step_id < num_sub_steps_per_part; sub_step_id++)
		{
			step_id = part_id * num_sub_steps_per_part + sub_step_id;

			global_point_id = period_id * num_sub_steps + dump_point_id;
			dump_point_id++;
			int dump_id = global_point_id + 1;

			ed->curr_time = double(dump_id) / double(num_sub_steps) * T;

			one_sub_period_deep(ad, 0, part_id, 0);
			calc_chars_std(ad, 0);
			calc_chars_lpn(ad, 0, 0);
			evo_chars_std(ad, 0, dump_id);
			evo_chars_lpn(ad, 0, dump_id);

#pragma omp parallel for
			for (int tr_id = 1; tr_id < num_trajectories; tr_id++)
			{
				int thread_id = omp_get_thread_num();
				one_sub_period_deep(ad, tr_id, part_id, thread_id);
				calc_chars_std(ad, tr_id);
				lambda_lpn_per_periods(ad, tmp, tr_id, 0, num_sub_steps, global_point_id, num_periods);
				evo_chars_std(ad, tr_id, dump_id);
				evo_chars_lpn(ad, tr_id, dump_id);
			}

#pragma omp parallel for
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				if (dump_evo_sep == 1)
				{
					dump_adr_single(ad, tr_id, true);
				}
			}
		}
	}

	delete tmp;
}

void MBLCoreBehaviour::ex_period_obs_deep_cd(AllData * ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

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

			double observable = ed->spec[tr_id].real;

			for (int cd_st_id = 0; cd_st_id < ed->cd_dim; cd_st_id++)
			{
				curr_point_id = global_point_id - cd_st_id;

				if (curr_point_id >= 0 && curr_point_id < ed->cd_num_points)
				{
					ed->cd_rec_data[tr_id][curr_point_id][cd_st_id] = observable;
				}
			}

			one_sub_period_deep(ad, tr_id, part_id, th_id);
			calc_chars_std(ad, tr_id);
		}
	}
}

void MBLCoreBehaviour::ex_period_obs_deep_sigma(AllData * ad, int tr_id, int th_id, int period_id) const
{
	stringstream msg;
	msg << "Error: need to add functionality" << endl;
	Error(msg.str());
}

void MBLCoreBehaviour::rk_period(AllData * ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int num_branches = md->num_ham_qj;

	double T = md->T;

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		double time = double(period_id) * T;
		
		for (int in_step_id = 0; in_step_id < cp->rk_ns; in_step_id++)
		{
			ed->times_all[tr_id] = time + double(in_step_id) * ed->rk_step;
			rk_step_mbl(ad, tr_id, th_id, ed->rk_step);
		}
	}
}

void MBLCoreBehaviour::rk_period_trp_deep(AllData * ad, int tr_id, int th_id, int period_id) const
{
	stringstream msg;
	msg << "Error: need to add functionality" << endl;
	Error(msg.str());
}

void MBLCoreBehaviour::rk_period_obs_deep(AllData * ad, int tr_id, int th_id, int period_id) const
{
	stringstream msg;
	msg << "Error: need to add functionality" << endl;
	Error(msg.str());
}

void MBLCoreBehaviour::rk_period_obs_deep_lpn(AllData * ad, int period_id) const
{
	stringstream msg;
	msg << "Error: need to add functionality" << endl;
	Error(msg.str());
}

void MBLCoreBehaviour::rk_period_obs_deep_cd(AllData * ad, int tr_id, int th_id, int period_id) const
{
	stringstream msg;
	msg << "Error: need to add functionality" << endl;
	Error(msg.str());
}

void MBLCoreBehaviour::rk_period_obs_deep_sigma(AllData * ad, int tr_id, int th_id, int period_id) const
{
	stringstream msg;
	msg << "Error: need to add functionality" << endl;
	Error(msg.str());
}

void MBLCoreBehaviour::calc_chars_std_start(AllData * ad, int tr_id) const
{
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;
	MKL_Complex16 * phi = &(ed->phi_all[tr_id * sys_size]);
	double * adr = &(ed->abs_diag_rho_all[tr_id * sys_size]);

	double norm = norm_square(phi, sys_size);

	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		adr[st_id] = mult_scalar_double(mult_scalar_complex(&phi[st_id], &phi[st_id], 1), 1.0 / norm).real;
	}

	MKL_Complex16 spec = get_spec_mbl(ad, tr_id);
	double imbalance = get_imbalance_mbl(ad, adr);

	ed->norm[tr_id] = norm;
	ed->spec[tr_id] = spec;
	ed->mean[tr_id] = imbalance;

	calc_random_obs(ad, tr_id);
}

void MBLCoreBehaviour::calc_chars_std(AllData * ad, int tr_id) const
{
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;
	MKL_Complex16 * phi = &(ed->phi_all[tr_id * sys_size]);
	double * adr = &(ed->abs_diag_rho_all[tr_id * sys_size]);

	double norm = norm_square(phi, sys_size);

	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		adr[st_id] = mult_scalar_double(mult_scalar_complex(&phi[st_id], &phi[st_id], 1), 1.0 / norm).real;
	}

	MKL_Complex16 spec = get_spec_mbl(ad, tr_id);
	double imbalance = get_imbalance_mbl(ad, adr);

	ed->norm[tr_id] = norm;
	ed->spec[tr_id] = spec;
	ed->mean[tr_id] = imbalance;

	calc_random_obs(ad, tr_id);
}

void MBLCoreBehaviour::calc_chars_lpn_start(AllData * ad, int tr_id, int base_tr_id) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;


	double lambda = 0.0;
	double lambda_now = 0.0;
	MKL_Complex16 spec_lpn = ed->spec[tr_id];
	double imbalance_lpn = ed->mean[tr_id];

	double delta_s = this->calc_delta_f(ad, tr_id, base_tr_id); // Important! Here we use calc_delta_f not calc_delta_s

	ed->lambda[tr_id] = lambda;
	ed->lambda_now[tr_id] = lambda_now;
	ed->delta_s[tr_id] = delta_s;

	ed->spec_lpn[tr_id] = spec_lpn;
	ed->mean_lpn[tr_id] = imbalance_lpn;

	calc_random_obs_lpn_start(ad, tr_id);
}

void MBLCoreBehaviour::calc_chars_lpn(AllData * ad, int tr_id, int base_tr_id) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;
	MKL_Complex16 * phi = &(ed->phi_all[tr_id * sys_size]);
	double * adr = &(ed->abs_diag_rho_all[tr_id * sys_size]);
	double norm = norm_square(phi, sys_size);
	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		adr[st_id] = mult_scalar_double(mult_scalar_complex(&phi[st_id], &phi[st_id], 1), 1.0 / norm).real;
	}

	MKL_Complex16 spec_lpn = get_spec_mbl(ad, tr_id);
	double imbalance_lpn = get_imbalance_mbl(ad, adr);

	ed->spec_lpn[tr_id] = spec_lpn;
	ed->mean_lpn[tr_id] = imbalance_lpn;

	calc_random_obs_lpn(ad, tr_id);
}

double MBLCoreBehaviour::calc_T(AllData * ad) const
{
	MainData * md = ad->md;
	return md->T;
}

void MBLCoreBehaviour::evo_chars_std(AllData * ad, int tr_id, int dump_id) const
{
	ConfigParam * cp = ad->cp;
	ExpData * ed = ad->ed;

	int dump_num_total = ed->dump_num_total;

	int index = tr_id * dump_num_total + dump_id;

	ed->norm_evo[index] = ed->norm[tr_id];
	ed->spec_evo[index] = ed->spec[tr_id];
	ed->mean_evo[index] = ed->mean[tr_id];

	int num_random_obs = int(cp->params.find("num_random_obs")->second);
	for (int obs_id = 0; obs_id < num_random_obs; obs_id++)
	{
		ed->random_obs_evo[tr_id][obs_id].push_back(ed->random_obs[tr_id][obs_id]);
	}
}

void MBLCoreBehaviour::evo_chars_lpn(AllData * ad, int tr_id, int dump_id) const
{
	ConfigParam * cp = ad->cp;
	ExpData * ed = ad->ed;

	int dump_num_total = ed->dump_num_total;

	int index = tr_id * dump_num_total + dump_id;

	ed->lambda_evo[index] = ed->lambda_now[tr_id];
	ed->spec_lpn_evo[index] = ed->spec_lpn[tr_id];
	ed->mean_lpn_evo[index] = ed->mean_lpn[tr_id];

	int num_random_obs = int(cp->params.find("num_random_obs")->second);
	for (int obs_id = 0; obs_id < num_random_obs; obs_id++)
	{
		ed->random_obs_lpn_evo[tr_id][obs_id].push_back(ed->random_obs_lpn[tr_id][obs_id]);
	}
}

double MBLCoreBehaviour::calc_delta_s(AllData * ad, int tr_id, int base_tr_id) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;

	int num_random_obs = int(ad->cp->params.find("num_random_obs")->second);
	int lpn_type = int(cp->params.find("lpn_type")->second);

	double delta_s = 0.0;
	if (lpn_type == -1)
	{
		MKL_Complex16 base = ed->spec_lpn[base_tr_id];
		MKL_Complex16 var = ed->spec_lpn[tr_id];
		double tmp = pow((base.real - var.real), 2) + pow((base.imag - var.imag), 2);
		delta_s = sqrt(tmp);
	}
	else if (lpn_type == -2)
	{
		delta_s = fabs(ed->mean_lpn[tr_id] - ed->mean_lpn[base_tr_id]);
	}
	else if (lpn_type >= 0 && lpn_type < num_random_obs)
	{
		complex<double> base = ed->random_obs_lpn[base_tr_id][lpn_type];
		complex<double> var = ed->random_obs_lpn[tr_id][lpn_type];
		double tmp = pow((base.real() - var.real()), 2) + pow((base.imag() - var.imag()), 2);
		delta_s = sqrt(tmp);
	}
	else
	{
		stringstream msg;
		msg << "Error: Wrong lpn_type: " << lpn_type << endl;
		Error(msg.str());
	}

	return delta_s;
}

double MBLCoreBehaviour::calc_delta_f(AllData * ad, int tr_id, int base_tr_id) const
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;

	int num_random_obs = int(ad->cp->params.find("num_random_obs")->second);
	int lpn_type = int(cp->params.find("lpn_type")->second);

	double delta_f = 0.0;
	if (lpn_type == -1)
	{
		MKL_Complex16 base = ed->spec[base_tr_id];
		MKL_Complex16 var = ed->spec[tr_id];
		double tmp = pow((base.real - var.real), 2) + pow((base.imag - var.imag), 2);
		delta_f = sqrt(tmp);
	}
	else if (lpn_type == -2)
	{
		delta_f = fabs(ed->mean[tr_id] - ed->mean[base_tr_id]);
	}
	else if (lpn_type >= 0 && lpn_type < num_random_obs)
	{
		complex<double> base = ed->random_obs[base_tr_id][lpn_type];
		complex<double> var = ed->random_obs[tr_id][lpn_type];
		double tmp = pow((base.real() - var.real()), 2) + pow((base.imag() - var.imag()), 2);
		delta_f = sqrt(tmp);
	}
	else
	{
		stringstream msg;
		msg << "Error: Wrong lpn_type: " << lpn_type << endl;
		Error(msg.str());
	}

	return delta_f;
}

void MBLCoreBehaviour::calc_ci(AllData * ad, int tr_id) const
{
	stringstream msg;
	msg << "Error: need to add functionality" << endl;
	Error(msg.str());
}

void MBLCoreBehaviour::dump_std(AllData * ad) const
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	ExpData * ed = ad->ed;

	int dump_obs = int(cp->params.find("dump_obs")->second);

	if (dump_obs == 1)
	{
		int num_trajectories = cp->num_trajectories;

		double * norm = ed->norm;
		MKL_Complex16 * spec = ed->spec;
		double * imbalance = ed->mean;

		string fn;

		fn = rp->path + "norm" + cp->fn_suffix;
		save_double_data(fn, norm, num_trajectories, 16, false);

		fn = rp->path + "spec" + cp->fn_suffix;
		save_complex_data(fn, spec, num_trajectories, 16, false);

		fn = rp->path + "imbalance" + cp->fn_suffix;
		save_double_data(fn, imbalance, num_trajectories, 16, false);
	}
}

void MBLCoreBehaviour::dump_lpn(AllData * ad) const
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	ExpData * ed = ad->ed;

	int num_trajectories = cp->num_trajectories;

	int dump_obs = int(cp->params.find("dump_obs")->second);

	if (dump_obs == 1)
	{
		double * lambda = ed->lambda_now;
		MKL_Complex16 * spec_lpn = ed->spec_lpn;
		double * imbalance_lpn = ed->mean_lpn;

		string fn;

		fn = rp->path + "lambda" + cp->fn_suffix;
		save_double_data(fn, lambda, num_trajectories, 16, false);

		int* num_renorms = ed->num_renorms;
		fn = rp->path + "num_renorms" + cp->fn_suffix;
		save_int_data(fn, num_renorms, num_trajectories, false);

		fn = rp->path + "spec_lpn" + cp->fn_suffix;
		save_complex_data(fn, spec_lpn, num_trajectories, 16, false);

		fn = rp->path + "imbalance_lpn" + cp->fn_suffix;
		save_double_data(fn, imbalance_lpn, num_trajectories, 16, false);

		int save_lambdas = int(ad->cp->params.find("save_lambdas")->second);
		if (save_lambdas > 0)
		{
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				fn = rp->path + "delta_f_by_delta_s_" + to_string(tr_id) + cp->fn_suffix;
				save_double_vector(fn, ed->lambdas[tr_id], 16, false);

				fn = rp->path + "deltas_s_" + to_string(tr_id) + cp->fn_suffix;
				save_double_vector(fn, ed->deltas_s[tr_id], 16, false);

				fn = rp->path + "deltas_f_" + to_string(tr_id) + cp->fn_suffix;
				save_double_vector(fn, ed->deltas_f[tr_id], 16, false);
			}
		}
	}
}

void MBLCoreBehaviour::dump_std_evo(AllData * ad) const
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	ExpData * ed = ad->ed;

	int num_trajectories = cp->num_trajectories;
	int dump_num_total = ed->dump_num_total;

	int dump_obs = int(cp->params.find("dump_obs")->second);
	int jump = int(cp->params.find("jump")->second);

	if (dump_obs == 1)
	{
		int * dump_periods = ed->dump_periods;

		double * norm_evo = ed->norm_evo;
		MKL_Complex16 * spec_evo = ed->spec_evo;
		double * imbalance_evo = ed->mean_evo;

		string fn;

		fn = rp->path + "norm_evo" + cp->fn_suffix;
		save_2d_inv_double_data(fn, norm_evo, dump_num_total, num_trajectories, 16, false);

		fn = rp->path + "periods" + cp->fn_suffix;
		save_int_data(fn, dump_periods, dump_num_total, false);

		fn = rp->path + "spec_evo" + cp->fn_suffix;
		save_2d_inv_complex_data(fn, spec_evo, dump_num_total, num_trajectories, 16, false);

		fn = rp->path + "imbalance_evo" + cp->fn_suffix;
		save_2d_inv_double_data(fn, imbalance_evo, dump_num_total, num_trajectories, 16, false);

		int num_random_obs = int(cp->params.find("num_random_obs")->second);
		int random_obs_seed = int(cp->params.find("random_obs_seed")->second);
		for (int obs_id = 0; obs_id < num_random_obs; obs_id++)
		{
			fn = rp->path + "random_obs_evo_" + to_string(random_obs_seed + obs_id) + cp->fn_suffix;
			std::vector<std::vector<std::complex<double>>> random_obs_evo_to_dump;
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				random_obs_evo_to_dump.push_back(ed->random_obs_evo[tr_id][obs_id]);
			}
			save_2d_vector(fn, random_obs_evo_to_dump, 16);
		}

		if (jump > 0)
		{
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				fn = rp->path + "jump_times_" + to_string(tr_id) + cp->fn_suffix;
				save_double_vector(fn, ed->jump_times[tr_id], 16, false);

				fn = rp->path + "diss_types_" + to_string(tr_id) + cp->fn_suffix;
				save_int_vector(fn, ed->diss_types[tr_id], false);

				if (jump > 1)
				{
					fn = rp->path + "jump_norms_" + to_string(tr_id) + cp->fn_suffix;
					save_double_vector(fn, ed->jump_norms[tr_id], 16, false);

					fn = rp->path + "jump_etas_" + to_string(tr_id) + cp->fn_suffix;
					save_double_vector(fn, ed->jump_etas[tr_id], 16, false);
				}
			}
		}
	}
}

void MBLCoreBehaviour::dump_lpn_evo(AllData * ad) const
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	ExpData * ed = ad->ed;

	int num_trajectories = cp->num_trajectories;
	int dump_num_total = ed->dump_num_total;

	int dump_obs = int(cp->params.find("dump_obs")->second);

	if (dump_obs == 1)
	{
		double * lambda_evo = ed->lambda_evo;
		MKL_Complex16 * spec_lpn_evo = ed->spec_lpn_evo;
		double * imbalance_lpn_evo = ed->mean_lpn_evo;

		string fn;

		fn = rp->path + "lambda_evo" + cp->fn_suffix;
		save_2d_inv_double_data(fn, lambda_evo, dump_num_total, num_trajectories, 16, false);

		fn = rp->path + "spec_lpn_evo" + cp->fn_suffix;
		save_2d_inv_complex_data(fn, spec_lpn_evo, dump_num_total, num_trajectories, 16, false);

		fn = rp->path + "imbalance_lpn_evo" + cp->fn_suffix;
		save_2d_inv_double_data(fn, imbalance_lpn_evo, dump_num_total, num_trajectories, 16, false);

		int num_random_obs = int(cp->params.find("num_random_obs")->second);
		int random_obs_seed = int(cp->params.find("random_obs_seed")->second);
		for (int obs_id = 0; obs_id < num_random_obs; obs_id++)
		{
			fn = rp->path + "random_obs_lpn_evo_" + to_string(random_obs_seed + obs_id) + cp->fn_suffix;
			std::vector<std::vector<std::complex<double>>> random_obs_evo_to_dump;
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				random_obs_evo_to_dump.push_back(ed->random_obs_lpn_evo[tr_id][obs_id]);
			}
			save_2d_vector(fn, random_obs_evo_to_dump, 16);
		}
	}
}










void LndHamCoreBehaviour::init_splits(AllData* ad) const
{
	RunParam* rp = ad->rp;
	MainData* md = ad->md;

	int num_branches = md->num_ham_qj;
	int num_threads = rp->num_threads;

	int num_total = num_threads * num_branches;

	md->structure = init_split_structure_lndham(ad);
	cout << "init_split_structure_lndham" << endl;
	md->splits = new Split[num_total];

	for (int b_id = 0; b_id < num_branches; b_id++)
	{
		for (int th_id = 0; th_id < num_threads; th_id++)
		{
			int index = b_id * num_threads + th_id;
			copy_struct_not_member(&(md->structure)[b_id], &(md->splits)[index]);
		}
	}
}

void LndHamCoreBehaviour::free_splits(AllData* ad) const
{
	RunParam* rp = ad->rp;
	MainData* md = ad->md;

	int num_branches = md->num_ham_qj;
	int num_threads = rp->num_threads;

	for (int b_id = 0; b_id < num_branches; b_id++)
	{
		for (int th_id = 0; th_id < num_threads; th_id++)
		{
			int index = b_id * num_threads + th_id;
			delete_split_struct_not_member(&(md->splits[index]));
		}
	}

	delete(md->splits);
	delete_split_struct(md->structure);
}

void LndHamCoreBehaviour::init_splits_deep(AllData* ad) const
{
	RunParam* rp = ad->rp;
	MainData* md = ad->md;

	int num_branches = md->num_ham_qj;
	int num_threads = rp->num_threads;

	int num_total = num_threads * num_branches;

	md->structure = init_split_structure_lndham_deep(ad);
	md->splits = new Split[num_total];

	for (int b_id = 0; b_id < num_branches; b_id++)
	{
		for (int th_id = 0; th_id < num_threads; th_id++)
		{
			int index = b_id * num_threads + th_id;
			copy_struct_not_member(&(md->structure)[b_id], &(md->splits)[index]);
		}
	}
}

void LndHamCoreBehaviour::free_splits_deep(AllData* ad) const
{
	RunParam* rp = ad->rp;
	MainData* md = ad->md;

	int num_branches = md->num_ham_qj;
	int num_threads = rp->num_threads;

	for (int b_id = 0; b_id < num_branches; b_id++)
	{
		for (int th_id = 0; th_id < num_threads; th_id++)
		{
			int index = b_id * num_threads + th_id;
			delete_split_struct_not_member(&(md->splits[index]));
		}
	}

	delete(md->splits);
	delete_split_struct(md->structure);
}

void LndHamCoreBehaviour::ex_period(AllData* ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;

	int num_branches = md->num_ham_qj;

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		one_sub_period_deep(ad, tr_id, part_id, th_id);
	}
}

void LndHamCoreBehaviour::ex_period_trp_deep(AllData* ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;

	int num_branches = md->num_ham_qj;
	int num_sub_steps = int(cp->params.find("deep_num_steps")->second);

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		for (int sub_step_id = 0; sub_step_id < num_sub_steps; sub_step_id++)
		{
			one_sub_period_deep(ad, tr_id, part_id, th_id);
		}
	}
}

void LndHamCoreBehaviour::ex_period_obs_deep(AllData* ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

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

			one_sub_period_deep(ad, tr_id, part_id, th_id);
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

void LndHamCoreBehaviour::ex_period_obs_deep_lpn(AllData* ad, int period_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	double T = md->T;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);

	int num_trajectories = cp->num_trajectories;

	int num_branches = md->num_ham_qj;
	int num_sub_steps_per_part = int(cp->params.find("deep_num_steps")->second);
	int num_sub_steps = num_branches * int(cp->params.find("deep_num_steps")->second);

	int dump_point_id = 0;
	int global_point_id = 0;

	int step_id = 0;

	CoreBehavior* tmp = new LndHamCoreBehaviour;

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		for (int sub_step_id = 0; sub_step_id < num_sub_steps_per_part; sub_step_id++)
		{
			step_id = part_id * num_sub_steps_per_part + sub_step_id;

			global_point_id = period_id * num_sub_steps + dump_point_id;
			dump_point_id++;
			int dump_id = global_point_id + 1;

			ed->curr_time = double(dump_id) / double(num_sub_steps) * T;

			one_sub_period_deep(ad, 0, part_id, 0);
			calc_chars_std(ad, 0);
			calc_chars_lpn(ad, 0, 0);
			evo_chars_std(ad, 0, dump_id);
			evo_chars_lpn(ad, 0, dump_id);

#pragma omp parallel for
			for (int tr_id = 1; tr_id < num_trajectories; tr_id++)
			{
				int thread_id = omp_get_thread_num();
				one_sub_period_deep(ad, tr_id, part_id, thread_id);
				calc_chars_std(ad, tr_id);
				lambda_lpn(ad, tmp, tr_id, 0);
				evo_chars_std(ad, tr_id, dump_id);
				evo_chars_lpn(ad, tr_id, dump_id);
			}

#pragma omp parallel for
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				if (dump_evo_sep == 1)
				{
					dump_adr_single(ad, tr_id, true);
				}
			}
		}
	}

	delete tmp;
}

void LndHamCoreBehaviour::ex_period_obs_deep_mult_lpn(AllData* ad, int period_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	double T = md->T;

	int lambda_per_periods = int(cp->params.find("lambda_per_periods")->second);

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);

	int num_trajectories = cp->num_trajectories;

	int num_branches = md->num_ham_qj;
	int num_sub_steps_per_part = int(cp->params.find("deep_num_steps")->second);
	int num_sub_steps = num_branches * int(cp->params.find("deep_num_steps")->second);

	int dump_point_id = 0;
	int global_point_id = 0;

	int step_id = 0;

	CoreBehavior* tmp = new LndHamCoreBehaviour;

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		for (int sub_step_id = 0; sub_step_id < num_sub_steps_per_part; sub_step_id++)
		{
			step_id = part_id * num_sub_steps_per_part + sub_step_id;

			global_point_id = period_id * num_sub_steps + dump_point_id;
			dump_point_id++;
			int dump_id = global_point_id + 1;

			ed->curr_time = double(dump_id) / double(num_sub_steps) * T;

#pragma omp parallel for
			for (int tr_id = 0; tr_id < num_trajectories / 2; tr_id++)
			{
				int thread_id = omp_get_thread_num();
				one_sub_period_deep(ad, tr_id, part_id, thread_id);
				calc_chars_std(ad, tr_id);
				calc_chars_lpn(ad, tr_id, tr_id);
				evo_chars_std(ad, tr_id, dump_id);
				evo_chars_lpn(ad, tr_id, dump_id);
			}

#pragma omp parallel for
			for (int tr_id = num_trajectories / 2; tr_id < num_trajectories; tr_id++)
			{
				int thread_id = omp_get_thread_num();
				one_sub_period_deep(ad, tr_id, part_id, thread_id);
				calc_chars_std(ad, tr_id);
				if (lambda_per_periods == 0)
				{
					lambda_lpn(ad, tmp, tr_id, tr_id - num_trajectories / 2);
				}
				if (lambda_per_periods == 1)
				{
					if ((part_id == num_branches - 1) && (sub_step_id == num_sub_steps_per_part - 1))
					{
						lambda_lpn_now(ad, tmp, tr_id, tr_id - num_trajectories / 2);
					}
					else
					{
						lambda_lpn_now_wo_renorm_just_to_save(ad, tmp, tr_id, tr_id - num_trajectories / 2);
					}
				}
				evo_chars_std(ad, tr_id, dump_id);
				evo_chars_lpn(ad, tr_id, dump_id);
			}

#pragma omp parallel for
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				if (dump_evo_sep == 1)
				{
					dump_adr_single(ad, tr_id, true);
				}
			}
		}
	}

	delete tmp;
}

void LndHamCoreBehaviour::ex_period_obs_deep_lpn_per_period(struct AllData* ad, int period_id, int num_periods) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	double T = md->T;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);

	int num_trajectories = cp->num_trajectories;

	int num_branches = md->num_ham_qj;
	int num_sub_steps_per_part = int(cp->params.find("deep_num_steps")->second);
	int num_sub_steps = num_branches * int(cp->params.find("deep_num_steps")->second);

	int dump_point_id = 0;
	int global_point_id = 0;

	int step_id = 0;

	CoreBehavior* tmp = new LndHamCoreBehaviour;

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		for (int sub_step_id = 0; sub_step_id < num_sub_steps_per_part; sub_step_id++)
		{
			step_id = part_id * num_sub_steps_per_part + sub_step_id;

			global_point_id = period_id * num_sub_steps + dump_point_id;
			dump_point_id++;
			int dump_id = global_point_id + 1;

			ed->curr_time = double(dump_id) / double(num_sub_steps) * T;

			one_sub_period_deep(ad, 0, part_id, 0);
			calc_chars_std(ad, 0);
			calc_chars_lpn(ad, 0, 0);
			evo_chars_std(ad, 0, dump_id);
			evo_chars_lpn(ad, 0, dump_id);

#pragma omp parallel for
			for (int tr_id = 1; tr_id < num_trajectories; tr_id++)
			{
				int thread_id = omp_get_thread_num();
				one_sub_period_deep(ad, tr_id, part_id, thread_id);
				calc_chars_std(ad, tr_id);
				lambda_lpn_per_periods(ad, tmp, tr_id, 0, num_sub_steps, global_point_id, num_periods);
				evo_chars_std(ad, tr_id, dump_id);
				evo_chars_lpn(ad, tr_id, dump_id);
			}

#pragma omp parallel for
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				if (dump_evo_sep == 1)
				{
					dump_adr_single(ad, tr_id, true);
				}
			}
		}
	}

	delete tmp;
}

void LndHamCoreBehaviour::ex_period_obs_deep_cd(AllData* ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

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

			double observable = ed->spec[tr_id].real;

			for (int cd_st_id = 0; cd_st_id < ed->cd_dim; cd_st_id++)
			{
				curr_point_id = global_point_id - cd_st_id;

				if (curr_point_id >= 0 && curr_point_id < ed->cd_num_points)
				{
					ed->cd_rec_data[tr_id][curr_point_id][cd_st_id] = observable;
				}
			}

			one_sub_period_deep(ad, tr_id, part_id, th_id);
			calc_chars_std(ad, tr_id);
		}
	}
}

void LndHamCoreBehaviour::ex_period_obs_deep_sigma(AllData* ad, int tr_id, int th_id, int period_id) const
{
	stringstream msg;
	msg << "Error: need to add functionality" << endl;
	Error(msg.str());
}

void LndHamCoreBehaviour::rk_period(AllData* ad, int tr_id, int th_id, int period_id) const
{
	stringstream msg;
	msg << "Error: need to add functionality" << endl;
	Error(msg.str());
}

void LndHamCoreBehaviour::rk_period_trp_deep(AllData* ad, int tr_id, int th_id, int period_id) const
{
	stringstream msg;
	msg << "Error: need to add functionality" << endl;
	Error(msg.str());
}

void LndHamCoreBehaviour::rk_period_obs_deep(AllData* ad, int tr_id, int th_id, int period_id) const
{
	stringstream msg;
	msg << "Error: need to add functionality" << endl;
	Error(msg.str());
}

void LndHamCoreBehaviour::rk_period_obs_deep_lpn(AllData* ad, int period_id) const
{
	stringstream msg;
	msg << "Error: need to add functionality" << endl;
	Error(msg.str());
}

void LndHamCoreBehaviour::rk_period_obs_deep_cd(AllData* ad, int tr_id, int th_id, int period_id) const
{
	stringstream msg;
	msg << "Error: need to add functionality" << endl;
	Error(msg.str());
}

void LndHamCoreBehaviour::rk_period_obs_deep_sigma(AllData* ad, int tr_id, int th_id, int period_id) const
{
	stringstream msg;
	msg << "Error: need to add functionality" << endl;
	Error(msg.str());
}

void LndHamCoreBehaviour::calc_chars_std_start(AllData* ad, int tr_id) const
{
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int sys_size = md->sys_size;
	MKL_Complex16* phi = &(ed->phi_all[tr_id * sys_size]);
	double* adr = &(ed->abs_diag_rho_all[tr_id * sys_size]);

	double norm = norm_square(phi, sys_size);

	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		adr[st_id] = mult_scalar_double(mult_scalar_complex(&phi[st_id], &phi[st_id], 1), 1.0 / norm).real;
	}

	MKL_Complex16 spec = get_spec_lndham(ad, tr_id);

	ed->norm[tr_id] = norm;
	ed->spec[tr_id] = spec;

	calc_random_obs(ad, tr_id);
}

void LndHamCoreBehaviour::calc_chars_std(AllData* ad, int tr_id) const
{
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int sys_size = md->sys_size;
	MKL_Complex16* phi = &(ed->phi_all[tr_id * sys_size]);
	double* adr = &(ed->abs_diag_rho_all[tr_id * sys_size]);

	double norm = norm_square(phi, sys_size);

	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		adr[st_id] = mult_scalar_double(mult_scalar_complex(&phi[st_id], &phi[st_id], 1), 1.0 / norm).real;
	}

	MKL_Complex16 spec = get_spec_lndham(ad, tr_id);

	ed->norm[tr_id] = norm;
	ed->spec[tr_id] = spec;

	calc_random_obs(ad, tr_id);
}

void LndHamCoreBehaviour::calc_chars_lpn_start(AllData* ad, int tr_id, int base_tr_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;


	double lambda = 0.0;
	double lambda_now = 0.0;
	MKL_Complex16 spec_lpn = ed->spec[tr_id];

	double delta_s = this->calc_delta_f(ad, tr_id, base_tr_id); // Important! Here we use calc_delta_f not calc_delta_s

	ed->lambda[tr_id] = lambda;
	ed->lambda_now[tr_id] = lambda_now;
	ed->delta_s[tr_id] = delta_s;

	ed->spec_lpn[tr_id] = spec_lpn;

	calc_random_obs_lpn_start(ad, tr_id);
}

void LndHamCoreBehaviour::calc_chars_lpn(AllData* ad, int tr_id, int base_tr_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int sys_size = md->sys_size;
	MKL_Complex16* phi = &(ed->phi_all[tr_id * sys_size]);
	double* adr = &(ed->abs_diag_rho_all[tr_id * sys_size]);
	double norm = norm_square(phi, sys_size);
	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		adr[st_id] = mult_scalar_double(mult_scalar_complex(&phi[st_id], &phi[st_id], 1), 1.0 / norm).real;
	}

	MKL_Complex16 spec_lpn = get_spec_lndham(ad, tr_id);

	ed->spec_lpn[tr_id] = spec_lpn;

	calc_random_obs_lpn(ad, tr_id);
}

double LndHamCoreBehaviour::calc_T(AllData* ad) const
{
	MainData* md = ad->md;
	return md->T;
}

void LndHamCoreBehaviour::evo_chars_std(AllData* ad, int tr_id, int dump_id) const
{
	ConfigParam* cp = ad->cp;
	ExpData* ed = ad->ed;

	int dump_num_total = ed->dump_num_total;

	int index = tr_id * dump_num_total + dump_id;

	ed->norm_evo[index] = ed->norm[tr_id];
	ed->spec_evo[index] = ed->spec[tr_id];

	int num_random_obs = int(cp->params.find("num_random_obs")->second);
	for (int obs_id = 0; obs_id < num_random_obs; obs_id++)
	{
		ed->random_obs_evo[tr_id][obs_id].push_back(ed->random_obs[tr_id][obs_id]);
	}
}

void LndHamCoreBehaviour::evo_chars_lpn(AllData* ad, int tr_id, int dump_id) const
{
	ConfigParam* cp = ad->cp;
	ExpData* ed = ad->ed;

	int dump_num_total = ed->dump_num_total;

	int index = tr_id * dump_num_total + dump_id;

	ed->lambda_evo[index] = ed->lambda_now[tr_id];
	ed->spec_lpn_evo[index] = ed->spec_lpn[tr_id];

	int num_random_obs = int(cp->params.find("num_random_obs")->second);
	for (int obs_id = 0; obs_id < num_random_obs; obs_id++)
	{
		ed->random_obs_lpn_evo[tr_id][obs_id].push_back(ed->random_obs_lpn[tr_id][obs_id]);
	}
}

double LndHamCoreBehaviour::calc_delta_s(AllData* ad, int tr_id, int base_tr_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int sys_size = md->sys_size;

	int num_random_obs = int(ad->cp->params.find("num_random_obs")->second);
	int lpn_type = int(cp->params.find("lpn_type")->second);

	double delta_s = 0.0;
	if (lpn_type == -1)
	{
		MKL_Complex16 base = ed->spec_lpn[base_tr_id];
		MKL_Complex16 var = ed->spec_lpn[tr_id];
		double tmp = pow((base.real - var.real), 2) + pow((base.imag - var.imag), 2);
		delta_s = sqrt(tmp);
	}
	else if (lpn_type >= 0 && lpn_type < num_random_obs)
	{
		complex<double> base = ed->random_obs_lpn[base_tr_id][lpn_type];
		complex<double> var = ed->random_obs_lpn[tr_id][lpn_type];
		double tmp = pow((base.real() - var.real()), 2) + pow((base.imag() - var.imag()), 2);
		delta_s = sqrt(tmp);
	}
	else
	{
		stringstream msg;
		msg << "Error: Wrong lpn_type: " << lpn_type << endl;
		Error(msg.str());
	}

	return delta_s;
}

double LndHamCoreBehaviour::calc_delta_f(AllData* ad, int tr_id, int base_tr_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int sys_size = md->sys_size;

	int num_random_obs = int(ad->cp->params.find("num_random_obs")->second);
	int lpn_type = int(cp->params.find("lpn_type")->second);

	double delta_f = 0.0;
	if (lpn_type == -1)
	{
		MKL_Complex16 base = ed->spec[base_tr_id];
		MKL_Complex16 var = ed->spec[tr_id];
		double tmp = pow((base.real - var.real), 2) + pow((base.imag - var.imag), 2);
		delta_f = sqrt(tmp);
	}
	else if (lpn_type >= 0 && lpn_type < num_random_obs)
	{
		complex<double> base = ed->random_obs[base_tr_id][lpn_type];
		complex<double> var = ed->random_obs[tr_id][lpn_type];
		double tmp = pow((base.real() - var.real()), 2) + pow((base.imag() - var.imag()), 2);
		delta_f = sqrt(tmp);
	}
	else
	{
		stringstream msg;
		msg << "Error: Wrong lpn_type: " << lpn_type << endl;
		Error(msg.str());
	}

	return delta_f;
}

void LndHamCoreBehaviour::calc_ci(AllData* ad, int tr_id) const
{
	stringstream msg;
	msg << "Error: need to add functionality" << endl;
	Error(msg.str());
}

void LndHamCoreBehaviour::dump_std(AllData* ad) const
{
	RunParam* rp = ad->rp;
	ConfigParam* cp = ad->cp;
	ExpData* ed = ad->ed;

	int dump_obs = int(cp->params.find("dump_obs")->second);

	if (dump_obs == 1)
	{
		int num_trajectories = cp->num_trajectories;

		double* norm = ed->norm;
		MKL_Complex16* spec = ed->spec;

		string fn;

		fn = rp->path + "norm" + cp->fn_suffix;
		save_double_data(fn, norm, num_trajectories, 16, false);

		fn = rp->path + "spec" + cp->fn_suffix;
		save_complex_data(fn, spec, num_trajectories, 16, false);
	}
}

void LndHamCoreBehaviour::dump_lpn(AllData* ad) const
{
	RunParam* rp = ad->rp;
	ConfigParam* cp = ad->cp;
	ExpData* ed = ad->ed;

	int num_trajectories = cp->num_trajectories;

	int dump_obs = int(cp->params.find("dump_obs")->second);

	if (dump_obs == 1)
	{
		double* lambda = ed->lambda_now;
		MKL_Complex16* spec_lpn = ed->spec_lpn;

		string fn;

		fn = rp->path + "lambda" + cp->fn_suffix;
		save_double_data(fn, lambda, num_trajectories, 16, false);

		int* num_renorms = ed->num_renorms;
		fn = rp->path + "num_renorms" + cp->fn_suffix;
		save_int_data(fn, num_renorms, num_trajectories, false);

		fn = rp->path + "spec_lpn" + cp->fn_suffix;
		save_complex_data(fn, spec_lpn, num_trajectories, 16, false);

		int save_lambdas = int(ad->cp->params.find("save_lambdas")->second);
		if (save_lambdas > 0)
		{
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				fn = rp->path + "delta_f_by_delta_s_" + to_string(tr_id) + cp->fn_suffix;
				save_double_vector(fn, ed->lambdas[tr_id], 16, false);

				fn = rp->path + "deltas_s_" + to_string(tr_id) + cp->fn_suffix;
				save_double_vector(fn, ed->deltas_s[tr_id], 16, false);

				fn = rp->path + "deltas_f_" + to_string(tr_id) + cp->fn_suffix;
				save_double_vector(fn, ed->deltas_f[tr_id], 16, false);
			}
		}
	}
}

void LndHamCoreBehaviour::dump_std_evo(AllData* ad) const
{
	RunParam* rp = ad->rp;
	ConfigParam* cp = ad->cp;
	ExpData* ed = ad->ed;

	int num_trajectories = cp->num_trajectories;
	int dump_num_total = ed->dump_num_total;

	int dump_obs = int(cp->params.find("dump_obs")->second);
	int jump = int(cp->params.find("jump")->second);

	if (dump_obs == 1)
	{
		int* dump_periods = ed->dump_periods;

		double* norm_evo = ed->norm_evo;
		MKL_Complex16* spec_evo = ed->spec_evo;

		string fn;

		fn = rp->path + "norm_evo" + cp->fn_suffix;
		save_2d_inv_double_data(fn, norm_evo, dump_num_total, num_trajectories, 16, false);

		fn = rp->path + "periods" + cp->fn_suffix;
		save_int_data(fn, dump_periods, dump_num_total, false);

		fn = rp->path + "spec_evo" + cp->fn_suffix;
		save_2d_inv_complex_data(fn, spec_evo, dump_num_total, num_trajectories, 16, false);

		int num_random_obs = int(cp->params.find("num_random_obs")->second);
		int random_obs_seed = int(cp->params.find("random_obs_seed")->second);
		for (int obs_id = 0; obs_id < num_random_obs; obs_id++)
		{
			fn = rp->path + "random_obs_evo_" + to_string(random_obs_seed + obs_id) + cp->fn_suffix;
			std::vector<std::vector<std::complex<double>>> random_obs_evo_to_dump;
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				random_obs_evo_to_dump.push_back(ed->random_obs_evo[tr_id][obs_id]);
			}
			save_2d_vector(fn, random_obs_evo_to_dump, 16);
		}

		if (jump > 0)
		{
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				fn = rp->path + "jump_times_" + to_string(tr_id) + cp->fn_suffix;
				save_double_vector(fn, ed->jump_times[tr_id], 16, false);

				fn = rp->path + "diss_types_" + to_string(tr_id) + cp->fn_suffix;
				save_int_vector(fn, ed->diss_types[tr_id], false);

				if (jump > 1)
				{
					fn = rp->path + "jump_norms_" + to_string(tr_id) + cp->fn_suffix;
					save_double_vector(fn, ed->jump_norms[tr_id], 16, false);

					fn = rp->path + "jump_etas_" + to_string(tr_id) + cp->fn_suffix;
					save_double_vector(fn, ed->jump_etas[tr_id], 16, false);
				}
			}
		}
	}
}

void LndHamCoreBehaviour::dump_lpn_evo(AllData* ad) const
{
	RunParam* rp = ad->rp;
	ConfigParam* cp = ad->cp;
	ExpData* ed = ad->ed;

	int num_trajectories = cp->num_trajectories;
	int dump_num_total = ed->dump_num_total;

	int dump_obs = int(cp->params.find("dump_obs")->second);

	if (dump_obs == 1)
	{
		double* lambda_evo = ed->lambda_evo;
		MKL_Complex16* spec_lpn_evo = ed->spec_lpn_evo;

		string fn;

		fn = rp->path + "lambda_evo" + cp->fn_suffix;
		save_2d_inv_double_data(fn, lambda_evo, dump_num_total, num_trajectories, 16, false);

		fn = rp->path + "spec_lpn_evo" + cp->fn_suffix;
		save_2d_inv_complex_data(fn, spec_lpn_evo, dump_num_total, num_trajectories, 16, false);

		int num_random_obs = int(cp->params.find("num_random_obs")->second);
		int random_obs_seed = int(cp->params.find("random_obs_seed")->second);
		for (int obs_id = 0; obs_id < num_random_obs; obs_id++)
		{
			fn = rp->path + "random_obs_lpn_evo_" + to_string(random_obs_seed + obs_id) + cp->fn_suffix;
			std::vector<std::vector<std::complex<double>>> random_obs_evo_to_dump;
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				random_obs_evo_to_dump.push_back(ed->random_obs_lpn_evo[tr_id][obs_id]);
			}
			save_2d_vector(fn, random_obs_evo_to_dump, 16);
		}
	}
}










void IntegrableCoreBehaviour::init_splits(AllData* ad) const
{
	RunParam* rp = ad->rp;
	MainData* md = ad->md;

	int num_branches = md->num_ham_qj;
	int num_threads = rp->num_threads;

	int num_total = num_threads * num_branches;

	md->structure = init_split_structure_integrable(ad);
	cout << "init_split_structure_integrable" << endl;
	md->splits = new Split[num_total];

	for (int b_id = 0; b_id < num_branches; b_id++)
	{
		for (int th_id = 0; th_id < num_threads; th_id++)
		{
			int index = b_id * num_threads + th_id;
			copy_struct_not_member(&(md->structure)[b_id], &(md->splits)[index]);
		}
	}
}

void IntegrableCoreBehaviour::free_splits(AllData* ad) const
{
	RunParam* rp = ad->rp;
	MainData* md = ad->md;

	int num_branches = md->num_ham_qj;
	int num_threads = rp->num_threads;

	for (int b_id = 0; b_id < num_branches; b_id++)
	{
		for (int th_id = 0; th_id < num_threads; th_id++)
		{
			int index = b_id * num_threads + th_id;
			delete_split_struct_not_member(&(md->splits[index]));
		}
	}

	delete(md->splits);
	delete_split_struct(md->structure);
}

void IntegrableCoreBehaviour::init_splits_deep(AllData* ad) const
{
	RunParam* rp = ad->rp;
	MainData* md = ad->md;

	int num_branches = md->num_ham_qj;
	int num_threads = rp->num_threads;

	int num_total = num_threads * num_branches;

	md->structure = init_split_structure_integrable_deep(ad);
	md->splits = new Split[num_total];

	for (int b_id = 0; b_id < num_branches; b_id++)
	{
		for (int th_id = 0; th_id < num_threads; th_id++)
		{
			int index = b_id * num_threads + th_id;
			copy_struct_not_member(&(md->structure)[b_id], &(md->splits)[index]);
		}
	}
}

void IntegrableCoreBehaviour::free_splits_deep(AllData* ad) const
{
	RunParam* rp = ad->rp;
	MainData* md = ad->md;

	int num_branches = md->num_ham_qj;
	int num_threads = rp->num_threads;

	for (int b_id = 0; b_id < num_branches; b_id++)
	{
		for (int th_id = 0; th_id < num_threads; th_id++)
		{
			int index = b_id * num_threads + th_id;
			delete_split_struct_not_member(&(md->splits[index]));
		}
	}

	delete(md->splits);
	delete_split_struct(md->structure);
}

void IntegrableCoreBehaviour::ex_period(AllData* ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;

	int num_branches = md->num_ham_qj;

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		one_sub_period_deep(ad, tr_id, part_id, th_id);
	}
}

void IntegrableCoreBehaviour::ex_period_trp_deep(AllData* ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;

	int num_branches = md->num_ham_qj;
	int num_sub_steps = int(cp->params.find("deep_num_steps")->second);

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		for (int sub_step_id = 0; sub_step_id < num_sub_steps; sub_step_id++)
		{
			one_sub_period_deep(ad, tr_id, part_id, th_id);
		}
	}
}

void IntegrableCoreBehaviour::ex_period_obs_deep(AllData* ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

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

			one_sub_period_deep(ad, tr_id, part_id, th_id);
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

void IntegrableCoreBehaviour::ex_period_obs_deep_lpn(AllData* ad, int period_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	double T = md->T;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);

	int num_trajectories = cp->num_trajectories;

	int num_branches = md->num_ham_qj;
	int num_sub_steps_per_part = int(cp->params.find("deep_num_steps")->second);
	int num_sub_steps = num_branches * int(cp->params.find("deep_num_steps")->second);

	int dump_point_id = 0;
	int global_point_id = 0;

	int step_id = 0;

	CoreBehavior* tmp = new IntegrableCoreBehaviour;

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		for (int sub_step_id = 0; sub_step_id < num_sub_steps_per_part; sub_step_id++)
		{
			step_id = part_id * num_sub_steps_per_part + sub_step_id;

			global_point_id = period_id * num_sub_steps + dump_point_id;
			dump_point_id++;
			int dump_id = global_point_id + 1;

			ed->curr_time = double(dump_id) / double(num_sub_steps) * T;

			one_sub_period_deep(ad, 0, part_id, 0);
			calc_chars_std(ad, 0);
			calc_chars_lpn(ad, 0, 0);
			evo_chars_std(ad, 0, dump_id);
			evo_chars_lpn(ad, 0, dump_id);

#pragma omp parallel for
			for (int tr_id = 1; tr_id < num_trajectories; tr_id++)
			{
				int thread_id = omp_get_thread_num();
				one_sub_period_deep(ad, tr_id, part_id, thread_id);
				calc_chars_std(ad, tr_id);
				lambda_lpn(ad, tmp, tr_id, 0);
				evo_chars_std(ad, tr_id, dump_id);
				evo_chars_lpn(ad, tr_id, dump_id);
			}

#pragma omp parallel for
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				if (dump_evo_sep == 1)
				{
					dump_adr_single(ad, tr_id, true);
				}
			}
		}
	}

	delete tmp;
}

void IntegrableCoreBehaviour::ex_period_obs_deep_mult_lpn(AllData* ad, int period_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	double T = md->T;

	int lambda_per_periods = int(cp->params.find("lambda_per_periods")->second);

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);

	int num_trajectories = cp->num_trajectories;

	int num_branches = md->num_ham_qj;
	int num_sub_steps_per_part = int(cp->params.find("deep_num_steps")->second);
	int num_sub_steps = num_branches * int(cp->params.find("deep_num_steps")->second);

	int dump_point_id = 0;
	int global_point_id = 0;

	int step_id = 0;

	CoreBehavior* tmp = new IntegrableCoreBehaviour;

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		for (int sub_step_id = 0; sub_step_id < num_sub_steps_per_part; sub_step_id++)
		{
			step_id = part_id * num_sub_steps_per_part + sub_step_id;

			global_point_id = period_id * num_sub_steps + dump_point_id;
			dump_point_id++;
			int dump_id = global_point_id + 1;

			ed->curr_time = double(dump_id) / double(num_sub_steps) * T;

#pragma omp parallel for
			for (int tr_id = 0; tr_id < num_trajectories / 2; tr_id++)
			{
				int thread_id = omp_get_thread_num();
				one_sub_period_deep(ad, tr_id, part_id, thread_id);
				calc_chars_std(ad, tr_id);
				calc_chars_lpn(ad, tr_id, tr_id);
				evo_chars_std(ad, tr_id, dump_id);
				evo_chars_lpn(ad, tr_id, dump_id);
			}

#pragma omp parallel for
			for (int tr_id = num_trajectories / 2; tr_id < num_trajectories; tr_id++)
			{
				int thread_id = omp_get_thread_num();
				one_sub_period_deep(ad, tr_id, part_id, thread_id);
				calc_chars_std(ad, tr_id);
				if (lambda_per_periods == 0)
				{
					lambda_lpn(ad, tmp, tr_id, tr_id - num_trajectories / 2);
				}
				if (lambda_per_periods == 1)
				{
					if ((part_id == num_branches - 1) && (sub_step_id == num_sub_steps_per_part - 1))
					{
						lambda_lpn_now(ad, tmp, tr_id, tr_id - num_trajectories / 2);
					}
					else
					{
						lambda_lpn_now_wo_renorm_just_to_save(ad, tmp, tr_id, tr_id - num_trajectories / 2);
					}
				}
				evo_chars_std(ad, tr_id, dump_id);
				evo_chars_lpn(ad, tr_id, dump_id);
			}

			dump_phi_evo(ad, true);

#pragma omp parallel for
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				if (dump_evo_sep == 1)
				{
					dump_adr_single(ad, tr_id, true);
				}
			}
		}
	}

	delete tmp;
}

void IntegrableCoreBehaviour::ex_period_obs_deep_lpn_per_period(struct AllData* ad, int period_id, int num_periods) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	double T = md->T;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);

	int num_trajectories = cp->num_trajectories;

	int num_branches = md->num_ham_qj;
	int num_sub_steps_per_part = int(cp->params.find("deep_num_steps")->second);
	int num_sub_steps = num_branches * int(cp->params.find("deep_num_steps")->second);

	int dump_point_id = 0;
	int global_point_id = 0;

	int step_id = 0;

	CoreBehavior* tmp = new IntegrableCoreBehaviour;

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		for (int sub_step_id = 0; sub_step_id < num_sub_steps_per_part; sub_step_id++)
		{
			step_id = part_id * num_sub_steps_per_part + sub_step_id;

			global_point_id = period_id * num_sub_steps + dump_point_id;
			dump_point_id++;
			int dump_id = global_point_id + 1;

			ed->curr_time = double(dump_id) / double(num_sub_steps) * T;

			one_sub_period_deep(ad, 0, part_id, 0);
			calc_chars_std(ad, 0);
			calc_chars_lpn(ad, 0, 0);
			evo_chars_std(ad, 0, dump_id);
			evo_chars_lpn(ad, 0, dump_id);

#pragma omp parallel for
			for (int tr_id = 1; tr_id < num_trajectories; tr_id++)
			{
				int thread_id = omp_get_thread_num();
				one_sub_period_deep(ad, tr_id, part_id, thread_id);
				calc_chars_std(ad, tr_id);
				lambda_lpn_per_periods(ad, tmp, tr_id, 0, num_sub_steps, global_point_id, num_periods);
				evo_chars_std(ad, tr_id, dump_id);
				evo_chars_lpn(ad, tr_id, dump_id);
			}

#pragma omp parallel for
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				if (dump_evo_sep == 1)
				{
					dump_adr_single(ad, tr_id, true);
				}
			}
		}
	}

	delete tmp;
}

void IntegrableCoreBehaviour::ex_period_obs_deep_cd(AllData* ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

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

			double observable = ed->spec[tr_id].real;

			for (int cd_st_id = 0; cd_st_id < ed->cd_dim; cd_st_id++)
			{
				curr_point_id = global_point_id - cd_st_id;

				if (curr_point_id >= 0 && curr_point_id < ed->cd_num_points)
				{
					ed->cd_rec_data[tr_id][curr_point_id][cd_st_id] = observable;
				}
			}

			one_sub_period_deep(ad, tr_id, part_id, th_id);
			calc_chars_std(ad, tr_id);
		}
	}
}

void IntegrableCoreBehaviour::ex_period_obs_deep_sigma(AllData* ad, int tr_id, int th_id, int period_id) const
{
	stringstream msg;
	msg << "Error: need to add functionality" << endl;
	Error(msg.str());
}

void IntegrableCoreBehaviour::rk_period(AllData* ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int num_branches = md->num_ham_qj;

	double T = md->T;

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		double time = double(period_id) * T;

		for (int in_step_id = 0; in_step_id < cp->rk_ns; in_step_id++)
		{
			ed->times_all[tr_id] = time + double(in_step_id) * ed->rk_step;
			rk_step_integrable(ad, tr_id, th_id, ed->rk_step);
		}
	}
}

void IntegrableCoreBehaviour::rk_period_trp_deep(AllData* ad, int tr_id, int th_id, int period_id) const
{
	stringstream msg;
	msg << "Error: need to add functionality" << endl;
	Error(msg.str());
}

void IntegrableCoreBehaviour::rk_period_obs_deep(AllData* ad, int tr_id, int th_id, int period_id) const
{
	stringstream msg;
	msg << "Error: need to add functionality" << endl;
	Error(msg.str());
}

void IntegrableCoreBehaviour::rk_period_obs_deep_lpn(AllData* ad, int period_id) const
{
	stringstream msg;
	msg << "Error: need to add functionality" << endl;
	Error(msg.str());
}

void IntegrableCoreBehaviour::rk_period_obs_deep_cd(AllData* ad, int tr_id, int th_id, int period_id) const
{
	stringstream msg;
	msg << "Error: need to add functionality" << endl;
	Error(msg.str());
}

void IntegrableCoreBehaviour::rk_period_obs_deep_sigma(AllData* ad, int tr_id, int th_id, int period_id) const
{
	stringstream msg;
	msg << "Error: need to add functionality" << endl;
	Error(msg.str());
}

void IntegrableCoreBehaviour::calc_chars_std_start(AllData* ad, int tr_id) const
{
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int sys_size = md->sys_size;
	MKL_Complex16* phi = &(ed->phi_all[tr_id * sys_size]);
	double* adr = &(ed->abs_diag_rho_all[tr_id * sys_size]);

	double norm = norm_square(phi, sys_size);

	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		adr[st_id] = mult_scalar_double(mult_scalar_complex(&phi[st_id], &phi[st_id], 1), 1.0 / norm).real;
	}

	MKL_Complex16 spec = get_spec_integrable(ad, tr_id);

	ed->norm[tr_id] = norm;
	ed->spec[tr_id] = spec;

	calc_random_obs(ad, tr_id);
}

void IntegrableCoreBehaviour::calc_chars_std(AllData* ad, int tr_id) const
{
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int sys_size = md->sys_size;
	MKL_Complex16* phi = &(ed->phi_all[tr_id * sys_size]);
	double* adr = &(ed->abs_diag_rho_all[tr_id * sys_size]);

	double norm = norm_square(phi, sys_size);

	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		adr[st_id] = mult_scalar_double(mult_scalar_complex(&phi[st_id], &phi[st_id], 1), 1.0 / norm).real;
	}

	MKL_Complex16 spec = get_spec_integrable(ad, tr_id);

	ed->norm[tr_id] = norm;
	ed->spec[tr_id] = spec;

	calc_random_obs(ad, tr_id);
}

void IntegrableCoreBehaviour::calc_chars_lpn_start(AllData* ad, int tr_id, int base_tr_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;


	double lambda = 0.0;
	double lambda_now = 0.0;
	MKL_Complex16 spec_lpn = ed->spec[tr_id];

	double delta_s = this->calc_delta_f(ad, tr_id, base_tr_id); // Important! Here we use calc_delta_f not calc_delta_s

	ed->lambda[tr_id] = lambda;
	ed->lambda_now[tr_id] = lambda_now;
	ed->delta_s[tr_id] = delta_s;

	ed->spec_lpn[tr_id] = spec_lpn;

	calc_random_obs_lpn_start(ad, tr_id);
}

void IntegrableCoreBehaviour::calc_chars_lpn(AllData* ad, int tr_id, int base_tr_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int sys_size = md->sys_size;
	MKL_Complex16* phi = &(ed->phi_all[tr_id * sys_size]);
	double* adr = &(ed->abs_diag_rho_all[tr_id * sys_size]);
	double norm = norm_square(phi, sys_size);
	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		adr[st_id] = mult_scalar_double(mult_scalar_complex(&phi[st_id], &phi[st_id], 1), 1.0 / norm).real;
	}

	MKL_Complex16 spec_lpn = get_spec_integrable(ad, tr_id);

	ed->spec_lpn[tr_id] = spec_lpn;

	calc_random_obs_lpn(ad, tr_id);
}

double IntegrableCoreBehaviour::calc_T(AllData* ad) const
{
	MainData* md = ad->md;
	return md->T;
}

void IntegrableCoreBehaviour::evo_chars_std(AllData* ad, int tr_id, int dump_id) const
{
	ConfigParam* cp = ad->cp;
	ExpData* ed = ad->ed;

	int dump_num_total = ed->dump_num_total;

	int index = tr_id * dump_num_total + dump_id;

	ed->norm_evo[index] = ed->norm[tr_id];
	ed->spec_evo[index] = ed->spec[tr_id];

	int num_random_obs = int(cp->params.find("num_random_obs")->second);
	for (int obs_id = 0; obs_id < num_random_obs; obs_id++)
	{
		ed->random_obs_evo[tr_id][obs_id].push_back(ed->random_obs[tr_id][obs_id]);
	}
}

void IntegrableCoreBehaviour::evo_chars_lpn(AllData* ad, int tr_id, int dump_id) const
{
	ConfigParam* cp = ad->cp;
	ExpData* ed = ad->ed;

	int dump_num_total = ed->dump_num_total;

	int index = tr_id * dump_num_total + dump_id;

	ed->lambda_evo[index] = ed->lambda_now[tr_id];
	ed->spec_lpn_evo[index] = ed->spec_lpn[tr_id];

	int num_random_obs = int(cp->params.find("num_random_obs")->second);
	for (int obs_id = 0; obs_id < num_random_obs; obs_id++)
	{
		ed->random_obs_lpn_evo[tr_id][obs_id].push_back(ed->random_obs_lpn[tr_id][obs_id]);
	}
}

double IntegrableCoreBehaviour::calc_delta_s(AllData* ad, int tr_id, int base_tr_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int sys_size = md->sys_size;

	int num_random_obs = int(ad->cp->params.find("num_random_obs")->second);
	int lpn_type = int(cp->params.find("lpn_type")->second);

	double delta_s = 0.0;
	if (lpn_type == -1)
	{
		MKL_Complex16 base = ed->spec_lpn[base_tr_id];
		MKL_Complex16 var = ed->spec_lpn[tr_id];
		double tmp = pow((base.real - var.real), 2) + pow((base.imag - var.imag), 2);
		delta_s = sqrt(tmp);
	}
	else if (lpn_type >= 0 && lpn_type < num_random_obs)
	{
		complex<double> base = ed->random_obs_lpn[base_tr_id][lpn_type];
		complex<double> var = ed->random_obs_lpn[tr_id][lpn_type];
		double tmp = pow((base.real() - var.real()), 2) + pow((base.imag() - var.imag()), 2);
		delta_s = sqrt(tmp);
	}
	else
	{
		stringstream msg;
		msg << "Error: Wrong lpn_type: " << lpn_type << endl;
		Error(msg.str());
	}

	return delta_s;
}

double IntegrableCoreBehaviour::calc_delta_f(AllData* ad, int tr_id, int base_tr_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int sys_size = md->sys_size;

	int num_random_obs = int(ad->cp->params.find("num_random_obs")->second);
	int lpn_type = int(cp->params.find("lpn_type")->second);

	double delta_f = 0.0;
	if (lpn_type == -1)
	{
		MKL_Complex16 base = ed->spec[base_tr_id];
		MKL_Complex16 var = ed->spec[tr_id];
		double tmp = pow((base.real - var.real), 2) + pow((base.imag - var.imag), 2);
		delta_f = sqrt(tmp);
	}
	else if (lpn_type >= 0 && lpn_type < num_random_obs)
	{
		complex<double> base = ed->random_obs[base_tr_id][lpn_type];
		complex<double> var = ed->random_obs[tr_id][lpn_type];
		double tmp = pow((base.real() - var.real()), 2) + pow((base.imag() - var.imag()), 2);
		delta_f = sqrt(tmp);
	}
	else
	{
		stringstream msg;
		msg << "Error: Wrong lpn_type: " << lpn_type << endl;
		Error(msg.str());
	}

	return delta_f;
}

void IntegrableCoreBehaviour::calc_ci(AllData* ad, int tr_id) const
{
	stringstream msg;
	msg << "Error: need to add functionality" << endl;
	Error(msg.str());
}

void IntegrableCoreBehaviour::dump_std(AllData* ad) const
{
	RunParam* rp = ad->rp;
	ConfigParam* cp = ad->cp;
	ExpData* ed = ad->ed;

	int dump_obs = int(cp->params.find("dump_obs")->second);

	if (dump_obs == 1)
	{
		int num_trajectories = cp->num_trajectories;

		double* norm = ed->norm;
		MKL_Complex16* spec = ed->spec;

		string fn;

		fn = rp->path + "norm" + cp->fn_suffix;
		save_double_data(fn, norm, num_trajectories, 16, false);

		fn = rp->path + "spec" + cp->fn_suffix;
		save_complex_data(fn, spec, num_trajectories, 16, false);
	}
}

void IntegrableCoreBehaviour::dump_lpn(AllData* ad) const
{
	RunParam* rp = ad->rp;
	ConfigParam* cp = ad->cp;
	ExpData* ed = ad->ed;

	int num_trajectories = cp->num_trajectories;

	int dump_obs = int(cp->params.find("dump_obs")->second);

	if (dump_obs == 1)
	{
		double* lambda = ed->lambda_now;
		MKL_Complex16* spec_lpn = ed->spec_lpn;

		string fn;

		fn = rp->path + "lambda" + cp->fn_suffix;
		save_double_data(fn, lambda, num_trajectories, 16, false);

		int* num_renorms = ed->num_renorms;
		fn = rp->path + "num_renorms" + cp->fn_suffix;
		save_int_data(fn, num_renorms, num_trajectories, false);

		fn = rp->path + "spec_lpn" + cp->fn_suffix;
		save_complex_data(fn, spec_lpn, num_trajectories, 16, false);

		int save_lambdas = int(ad->cp->params.find("save_lambdas")->second);
		if (save_lambdas > 0)
		{
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				fn = rp->path + "delta_f_by_delta_s_" + to_string(tr_id) + cp->fn_suffix;
				save_double_vector(fn, ed->lambdas[tr_id], 16, false);

				fn = rp->path + "deltas_s_" + to_string(tr_id) + cp->fn_suffix;
				save_double_vector(fn, ed->deltas_s[tr_id], 16, false);

				fn = rp->path + "deltas_f_" + to_string(tr_id) + cp->fn_suffix;
				save_double_vector(fn, ed->deltas_f[tr_id], 16, false);
			}
		}
	}
}

void IntegrableCoreBehaviour::dump_std_evo(AllData* ad) const
{
	RunParam* rp = ad->rp;
	ConfigParam* cp = ad->cp;
	ExpData* ed = ad->ed;

	int num_trajectories = cp->num_trajectories;
	int dump_num_total = ed->dump_num_total;

	int dump_obs = int(cp->params.find("dump_obs")->second);
	int jump = int(cp->params.find("jump")->second);

	if (dump_obs == 1)
	{
		int* dump_periods = ed->dump_periods;

		double* norm_evo = ed->norm_evo;
		MKL_Complex16* spec_evo = ed->spec_evo;

		string fn;

		fn = rp->path + "norm_evo" + cp->fn_suffix;
		save_2d_inv_double_data(fn, norm_evo, dump_num_total, num_trajectories, 16, false);

		fn = rp->path + "periods" + cp->fn_suffix;
		save_int_data(fn, dump_periods, dump_num_total, false);

		fn = rp->path + "spec_evo" + cp->fn_suffix;
		save_2d_inv_complex_data(fn, spec_evo, dump_num_total, num_trajectories, 16, false);

		int num_random_obs = int(cp->params.find("num_random_obs")->second);
		int random_obs_seed = int(cp->params.find("random_obs_seed")->second);
		for (int obs_id = 0; obs_id < num_random_obs; obs_id++)
		{
			fn = rp->path + "random_obs_evo_" + to_string(random_obs_seed + obs_id) + cp->fn_suffix;
			std::vector<std::vector<std::complex<double>>> random_obs_evo_to_dump;
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				random_obs_evo_to_dump.push_back(ed->random_obs_evo[tr_id][obs_id]);
			}
			save_2d_vector(fn, random_obs_evo_to_dump, 16);
		}

		if (jump > 0)
		{
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				fn = rp->path + "jump_times_" + to_string(tr_id) + cp->fn_suffix;
				save_double_vector(fn, ed->jump_times[tr_id], 16, false);

				fn = rp->path + "diss_types_" + to_string(tr_id) + cp->fn_suffix;
				save_int_vector(fn, ed->diss_types[tr_id], false);

				if (jump > 1)
				{
					fn = rp->path + "jump_norms_" + to_string(tr_id) + cp->fn_suffix;
					save_double_vector(fn, ed->jump_norms[tr_id], 16, false);

					fn = rp->path + "jump_etas_" + to_string(tr_id) + cp->fn_suffix;
					save_double_vector(fn, ed->jump_etas[tr_id], 16, false);
				}
			}
		}
	}
}

void IntegrableCoreBehaviour::dump_lpn_evo(AllData* ad) const
{
	RunParam* rp = ad->rp;
	ConfigParam* cp = ad->cp;
	ExpData* ed = ad->ed;

	int num_trajectories = cp->num_trajectories;
	int dump_num_total = ed->dump_num_total;

	int dump_obs = int(cp->params.find("dump_obs")->second);

	if (dump_obs == 1)
	{
		double* lambda_evo = ed->lambda_evo;
		MKL_Complex16* spec_lpn_evo = ed->spec_lpn_evo;

		string fn;

		fn = rp->path + "lambda_evo" + cp->fn_suffix;
		save_2d_inv_double_data(fn, lambda_evo, dump_num_total, num_trajectories, 16, false);

		fn = rp->path + "spec_lpn_evo" + cp->fn_suffix;
		save_2d_inv_complex_data(fn, spec_lpn_evo, dump_num_total, num_trajectories, 16, false);

		int num_random_obs = int(cp->params.find("num_random_obs")->second);
		int random_obs_seed = int(cp->params.find("random_obs_seed")->second);
		for (int obs_id = 0; obs_id < num_random_obs; obs_id++)
		{
			fn = rp->path + "random_obs_lpn_evo_" + to_string(random_obs_seed + obs_id) + cp->fn_suffix;
			std::vector<std::vector<std::complex<double>>> random_obs_evo_to_dump;
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				random_obs_evo_to_dump.push_back(ed->random_obs_lpn_evo[tr_id][obs_id]);
			}
			save_2d_vector(fn, random_obs_evo_to_dump, 16);
		}
	}
}






void Floq2SpinsCoreBehaviour::init_splits(AllData* ad) const
{
	RunParam* rp = ad->rp;
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;

	int num_threads = rp->num_threads;

	md->structure = init_split_structure_floq2spins(ad);
	md->splits = new Split[num_threads];
	for (int th_id = 0; th_id < num_threads; th_id++)
	{
		copy_struct_not_member(md->structure, &(md->splits)[th_id]);
	}

	cout << "Split initialization complete" << endl;
}

void Floq2SpinsCoreBehaviour::free_splits(AllData* ad) const
{
	RunParam* rp = ad->rp;
	MainData* md = ad->md;

	int num_threads = rp->num_threads;

	for (int i = 0; i < num_threads; i++)
	{
		delete_split_struct_not_member(&(md->splits[i]));
	}

	delete(md->splits);
	delete_split_struct(md->structure);
}

void Floq2SpinsCoreBehaviour::init_splits_deep(AllData* ad) const
{
	RunParam* rp = ad->rp;
	MainData* md = ad->md;

	int num_branches = md->num_ham_qj;
	int num_threads = rp->num_threads;

	int num_total = num_threads * num_branches;

	md->structure = init_split_structure_floq2spins_deep(ad);
	md->splits = new Split[num_total];

	for (int b_id = 0; b_id < num_branches; b_id++)
	{
		for (int th_id = 0; th_id < num_threads; th_id++)
		{
			int index = b_id * num_threads + th_id;
			copy_struct_not_member(&(md->structure)[b_id], &(md->splits)[index]);
		}
	}
}

void Floq2SpinsCoreBehaviour::free_splits_deep(AllData* ad) const
{
	RunParam* rp = ad->rp;
	MainData* md = ad->md;

	int num_branches = md->num_ham_qj;
	int num_threads = rp->num_threads;

	for (int b_id = 0; b_id < num_branches; b_id++)
	{
		for (int th_id = 0; th_id < num_threads; th_id++)
		{
			int index = b_id * num_threads + th_id;
			delete_split_struct_not_member(&(md->splits[index]));
		}
	}

	delete(md->splits);
	delete_split_struct(md->structure);
}

void Floq2SpinsCoreBehaviour::ex_period(AllData* ad, int tr_id, int th_id, int period_id) const
{
	MainData* md = ad->md;

	int sys_size = md->sys_size;
	Split* head = &(md->splits)[th_id];

	for (int b_id = 0; b_id < head->counter; b_id++)
	{
		one_period_branch(ad, head, tr_id, &(head->next)[b_id]);
	}
}

void Floq2SpinsCoreBehaviour::ex_period_trp_deep(AllData* ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;

	int num_branches = md->num_ham_qj;
	int num_sub_steps = int(cp->params.find("deep_num_steps")->second);

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		for (int sub_step_id = 0; sub_step_id < num_sub_steps; sub_step_id++)
		{
			one_sub_period_deep(ad, tr_id, part_id, th_id);
		}
	}
}

void Floq2SpinsCoreBehaviour::ex_period_obs_deep(AllData* ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

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

			one_sub_period_deep(ad, tr_id, part_id, th_id);
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

void Floq2SpinsCoreBehaviour::ex_period_obs_deep_lpn(AllData* ad, int period_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);

	int num_trajectories = cp->num_trajectories;

	int num_branches = md->num_ham_qj;
	int num_sub_steps_per_part = int(cp->params.find("deep_num_steps")->second);
	int num_sub_steps = num_branches * int(cp->params.find("deep_num_steps")->second);

	int dump_point_id = 0;
	int global_point_id = 0;

	int step_id = 0;

	CoreBehavior* tmp = new Floq2SpinsCoreBehaviour;

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		for (int sub_step_id = 0; sub_step_id < num_sub_steps_per_part; sub_step_id++)
		{
			step_id = part_id * num_sub_steps_per_part + sub_step_id;

			global_point_id = period_id * num_sub_steps + dump_point_id;
			dump_point_id++;
			int dump_id = global_point_id + 1;

			ed->curr_time = double(dump_id) / double(num_sub_steps) * md->T;

			one_sub_period_deep(ad, 0, part_id, 0);
			calc_chars_std(ad, 0);
			calc_chars_lpn(ad, 0, 0);
			evo_chars_std(ad, 0, dump_id);
			evo_chars_lpn(ad, 0, dump_id);

#pragma omp parallel for
			for (int tr_id = 1; tr_id < num_trajectories; tr_id++)
			{
				int thread_id = omp_get_thread_num();
				one_sub_period_deep(ad, tr_id, part_id, thread_id);
				calc_chars_std(ad, tr_id);
				lambda_lpn(ad, tmp, tr_id, 0);
				evo_chars_std(ad, tr_id, dump_id);
				evo_chars_lpn(ad, tr_id, dump_id);
			}

#pragma omp parallel for
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				if (dump_evo_sep == 1)
				{
					dump_adr_single(ad, tr_id, true);
				}
			}
		}
	}

	delete tmp;
}

void Floq2SpinsCoreBehaviour::ex_period_obs_deep_mult_lpn(AllData* ad, int period_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);

	int num_trajectories = cp->num_trajectories;

	int num_branches = md->num_ham_qj;
	int num_sub_steps_per_part = int(cp->params.find("deep_num_steps")->second);
	int num_sub_steps = num_branches * int(cp->params.find("deep_num_steps")->second);

	int lambda_per_periods = int(cp->params.find("lambda_per_periods")->second);

	int dump_point_id = 0;
	int global_point_id = 0;

	int step_id = 0;

	CoreBehavior* tmp = new Floq2SpinsCoreBehaviour;

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		for (int sub_step_id = 0; sub_step_id < num_sub_steps_per_part; sub_step_id++)
		{
			step_id = part_id * num_sub_steps_per_part + sub_step_id;

			global_point_id = period_id * num_sub_steps + dump_point_id;
			dump_point_id++;
			int dump_id = global_point_id + 1;

			ed->curr_time = double(dump_id) / double(num_sub_steps) * md->T;

#pragma omp parallel for
			for (int tr_id = 0; tr_id < num_trajectories / 2; tr_id++)
			{
				int thread_id = omp_get_thread_num();
				one_sub_period_deep(ad, tr_id, part_id, thread_id);
				calc_chars_std(ad, tr_id);
				calc_chars_lpn(ad, tr_id, tr_id);
				evo_chars_std(ad, tr_id, dump_id);
				evo_chars_lpn(ad, tr_id, dump_id);
			}

#pragma omp parallel for
			for (int tr_id = num_trajectories / 2; tr_id < num_trajectories; tr_id++)
			{
				int thread_id = omp_get_thread_num();
				one_sub_period_deep(ad, tr_id, part_id, thread_id);
				calc_chars_std(ad, tr_id);
				if (lambda_per_periods == 0)
				{
					lambda_lpn(ad, tmp, tr_id, tr_id - num_trajectories / 2);
				}
				if (lambda_per_periods == 1)
				{
					if ((part_id == num_branches - 1) && (sub_step_id == num_sub_steps_per_part - 1))
					{
						lambda_lpn_now(ad, tmp, tr_id, tr_id - num_trajectories / 2);
					}
				}
				evo_chars_std(ad, tr_id, dump_id);
				evo_chars_lpn(ad, tr_id, dump_id);
			}

#pragma omp parallel for
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				if (dump_evo_sep == 1)
				{
					dump_adr_single(ad, tr_id, true);
				}
			}
		}
	}

	delete tmp;
}

void Floq2SpinsCoreBehaviour::ex_period_obs_deep_lpn_per_period(struct AllData* ad, int period_id, int num_periods) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);

	int num_trajectories = cp->num_trajectories;

	int num_branches = md->num_ham_qj;
	int num_sub_steps_per_part = int(cp->params.find("deep_num_steps")->second);
	int num_sub_steps = num_branches * int(cp->params.find("deep_num_steps")->second);

	int dump_point_id = 0;
	int global_point_id = 0;

	int step_id = 0;

	CoreBehavior* tmp = new Floq2SpinsCoreBehaviour;

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		for (int sub_step_id = 0; sub_step_id < num_sub_steps_per_part; sub_step_id++)
		{
			step_id = part_id * num_sub_steps_per_part + sub_step_id;

			global_point_id = period_id * num_sub_steps + dump_point_id;
			dump_point_id++;
			int dump_id = global_point_id + 1;

			ed->curr_time = double(dump_id) / double(num_sub_steps) * md->T;

			one_sub_period_deep(ad, 0, part_id, 0);
			calc_chars_std(ad, 0);
			calc_chars_lpn(ad, 0, 0);
			evo_chars_std(ad, 0, dump_id);
			evo_chars_lpn(ad, 0, dump_id);

#pragma omp parallel for
			for (int tr_id = 1; tr_id < num_trajectories; tr_id++)
			{
				int thread_id = omp_get_thread_num();
				one_sub_period_deep(ad, tr_id, part_id, thread_id);
				calc_chars_std(ad, tr_id);
				lambda_lpn_per_periods(ad, tmp, tr_id, 0, num_sub_steps, global_point_id, num_periods);
				evo_chars_std(ad, tr_id, dump_id);
				evo_chars_lpn(ad, tr_id, dump_id);
			}

#pragma omp parallel for
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				if (dump_evo_sep == 1)
				{
					dump_adr_single(ad, tr_id, true);
				}
			}
		}
	}

	delete tmp;
}

void Floq2SpinsCoreBehaviour::ex_period_obs_deep_cd(AllData* ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

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

			one_sub_period_deep(ad, tr_id, part_id, th_id);
			calc_chars_std(ad, tr_id);
		}
	}
}

void Floq2SpinsCoreBehaviour::ex_period_obs_deep_sigma(AllData* ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;

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

			one_sub_period_deep(ad, tr_id, part_id, th_id);
			calc_chars_std(ad, tr_id);
		}
	}
}

void Floq2SpinsCoreBehaviour::rk_period(AllData* ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	for (int step_id = 0; step_id < cp->rk_ns; step_id++)
	{
		ed->times_all[tr_id] = double(period_id) * md->T + double(step_id) * ed->rk_step;
		rk_step_floq2spins(ad, tr_id, th_id, ed->rk_step);
	}
}

void Floq2SpinsCoreBehaviour::rk_period_trp_deep(AllData* ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int num_branches = md->num_ham_qj;
	int num_sub_steps_per_part = int(cp->params.find("deep_num_steps")->second);

	int step_id = 0;

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		for (int sub_step_id = 0; sub_step_id < num_sub_steps_per_part; sub_step_id++)
		{
			step_id = part_id * num_sub_steps_per_part + sub_step_id;

			for (int in_step_id = 0; in_step_id < cp->rk_ns; in_step_id++)
			{
				ed->times_all[tr_id] = double(period_id) * md->T + double(step_id * cp->rk_ns + in_step_id) * ed->rk_step;
				rk_step_floq2spins(ad, tr_id, th_id, ed->rk_step);
			}
		}
	}
}

void Floq2SpinsCoreBehaviour::rk_period_obs_deep(AllData* ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

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

			for (int in_step_id = 0; in_step_id < cp->rk_ns; in_step_id++)
			{
				ed->times_all[tr_id] = double(period_id) * md->T + double(step_id * cp->rk_ns + in_step_id) * ed->rk_step;
				rk_step_floq2spins(ad, tr_id, th_id, ed->rk_step);
			}

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

void Floq2SpinsCoreBehaviour::rk_period_obs_deep_lpn(AllData* ad, int period_id) const
{
	stringstream msg;
	msg << "Error: need to add functionality" << endl;
	Error(msg.str());
}

void Floq2SpinsCoreBehaviour::rk_period_obs_deep_cd(AllData* ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

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

			for (int in_step_id = 0; in_step_id < cp->rk_ns; in_step_id++)
			{
				ed->times_all[tr_id] = double(period_id) * md->T + double(step_id * cp->rk_ns + in_step_id) * ed->rk_step;
				rk_step_floq2spins(ad, tr_id, th_id, ed->rk_step);
			}

			calc_chars_std(ad, tr_id);
		}
	}
}

void Floq2SpinsCoreBehaviour::rk_period_obs_deep_sigma(AllData* ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

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

			dump_point_id++;

			for (int in_step_id = 0; in_step_id < cp->rk_ns; in_step_id++)
			{
				ed->times_all[tr_id] = double(period_id) * md->T + double(step_id * cp->rk_ns + in_step_id) * ed->rk_step;
				rk_step_floq2spins(ad, tr_id, th_id, ed->rk_step);
			}

			calc_chars_std(ad, tr_id);
		}
	}
}

void Floq2SpinsCoreBehaviour::calc_chars_std_start(AllData* ad, int tr_id) const
{
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int sys_size = md->sys_size;
	MKL_Complex16* phi = &(ed->phi_all[tr_id * sys_size]);
	double* adr = &(ed->abs_diag_rho_all[tr_id * sys_size]);

	double norm = norm_square(phi, sys_size);

	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		adr[st_id] = mult_scalar_double(mult_scalar_complex(&phi[st_id], &phi[st_id], 1), 1.0 / norm).real;
	}

	MKL_Complex16 spec = get_spec_flq2sp(ad, tr_id);

	ed->norm[tr_id] = norm;
	ed->spec[tr_id] = spec;

	calc_random_obs(ad, tr_id);
}

void Floq2SpinsCoreBehaviour::calc_chars_std(AllData* ad, int tr_id) const
{
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int sys_size = md->sys_size;
	MKL_Complex16* phi = &(ed->phi_all[tr_id * sys_size]);
	double* adr = &(ed->abs_diag_rho_all[tr_id * sys_size]);

	double norm = norm_square(phi, sys_size);

	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		adr[st_id] = mult_scalar_double(mult_scalar_complex(&phi[st_id], &phi[st_id], 1), 1.0 / norm).real;
	}

	MKL_Complex16 spec = get_spec_flq2sp(ad, tr_id);

	ed->norm[tr_id] = norm;
	ed->spec[tr_id] = spec;

	calc_random_obs(ad, tr_id);
}

void Floq2SpinsCoreBehaviour::calc_chars_lpn_start(AllData* ad, int tr_id, int base_tr_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	double lambda = 0.0;
	double lambda_now = 0.0;
	MKL_Complex16 spec_lpn = ed->spec[tr_id];

	double delta_s = this->calc_delta_f(ad, tr_id, base_tr_id); // Important! Here we use calc_delta_f not calc_delta_s

	ed->lambda[tr_id] = lambda;
	ed->lambda_now[tr_id] = lambda_now;
	ed->delta_s[tr_id] = delta_s;

	ed->spec_lpn[tr_id] = spec_lpn;

	calc_random_obs_lpn_start(ad, tr_id);
}

void Floq2SpinsCoreBehaviour::calc_chars_lpn(AllData* ad, int tr_id, int base_tr_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int sys_size = md->sys_size;
	MKL_Complex16* phi = &(ed->phi_all[tr_id * sys_size]);
	double* adr = &(ed->abs_diag_rho_all[tr_id * sys_size]);
	double norm_2 = norm_square(phi, sys_size);
	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		adr[st_id] = mult_scalar_double(mult_scalar_complex(&phi[st_id], &phi[st_id], 1), 1.0 / norm_2).real;
	}

	MKL_Complex16 spec_lpn = get_spec_flq2sp(ad, tr_id);

	ed->spec_lpn[tr_id] = spec_lpn;

	calc_random_obs_lpn(ad, tr_id);
}

double Floq2SpinsCoreBehaviour::calc_T(AllData* ad) const
{
	MainData* md = ad->md;
	double T = md->T;
	return T;
}

void Floq2SpinsCoreBehaviour::evo_chars_std(AllData* ad, int tr_id, int dump_id) const
{
	ConfigParam* cp = ad->cp;
	ExpData* ed = ad->ed;

	int dump_num_total = ed->dump_num_total;

	int index = tr_id * dump_num_total + dump_id;

	ed->norm_evo[index] = ed->norm[tr_id];
	ed->spec_evo[index] = ed->spec[tr_id];

	int num_random_obs = int(cp->params.find("num_random_obs")->second);
	for (int obs_id = 0; obs_id < num_random_obs; obs_id++)
	{
		ed->random_obs_evo[tr_id][obs_id].push_back(ed->random_obs[tr_id][obs_id]);
	}
}

void Floq2SpinsCoreBehaviour::evo_chars_lpn(AllData* ad, int tr_id, int dump_id) const
{
	ConfigParam* cp = ad->cp;
	ExpData* ed = ad->ed;

	int dump_num_total = ed->dump_num_total;

	int index = tr_id * dump_num_total + dump_id;

	ed->lambda_evo[index] = ed->lambda_now[tr_id];
	ed->spec_lpn_evo[index] = ed->spec_lpn[tr_id];

	int num_random_obs = int(cp->params.find("num_random_obs")->second);
	for (int obs_id = 0; obs_id < num_random_obs; obs_id++)
	{
		ed->random_obs_lpn_evo[tr_id][obs_id].push_back(ed->random_obs_lpn[tr_id][obs_id]);
	}
}

double Floq2SpinsCoreBehaviour::calc_delta_s(AllData* ad, int tr_id, int base_tr_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int sys_size = md->sys_size;

	int lpn_type = int(cp->params.find("lpn_type")->second);
	int num_random_obs = int(ad->cp->params.find("num_random_obs")->second);

	double delta_s = 0.0;
	if (lpn_type == -1)
	{
		MKL_Complex16 base = ed->spec_lpn[base_tr_id];
		MKL_Complex16 var = ed->spec_lpn[tr_id];
		double tmp = pow((base.real - var.real), 2) + pow((base.imag - var.imag), 2);
		delta_s = sqrt(tmp);
	}
	else if (lpn_type >= 0 && lpn_type < num_random_obs)
	{
		complex<double> base = ed->random_obs_lpn[base_tr_id][lpn_type];
		complex<double> var = ed->random_obs_lpn[tr_id][lpn_type];
		double tmp = pow((base.real() - var.real()), 2) + pow((base.imag() - var.imag()), 2);
		delta_s = sqrt(tmp);
	}
	else
	{
		stringstream msg;
		msg << "Error: Wrong lpn_type: " << lpn_type << endl;
		Error(msg.str());
	}

	return delta_s;
}

double Floq2SpinsCoreBehaviour::calc_delta_f(AllData* ad, int tr_id, int base_tr_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int sys_size = md->sys_size;

	int lpn_type = int(cp->params.find("lpn_type")->second);
	int num_random_obs = int(ad->cp->params.find("num_random_obs")->second);

	double delta_f = 0.0;
	if (lpn_type == -1)
	{
		MKL_Complex16 base = ed->spec[base_tr_id];
		MKL_Complex16 var = ed->spec[tr_id];
		double tmp = pow((base.real - var.real), 2) + pow((base.imag - var.imag), 2);
		delta_f = sqrt(tmp);
	}
	else if (lpn_type >= 0 && lpn_type < num_random_obs)
	{
		complex<double> base = ed->random_obs[base_tr_id][lpn_type];
		complex<double> var = ed->random_obs[tr_id][lpn_type];
		double tmp = pow((base.real() - var.real()), 2) + pow((base.imag() - var.imag()), 2);
		delta_f = sqrt(tmp);
	}
	else
	{
		stringstream msg;
		msg << "Error: Wrong lpn_type: " << lpn_type << endl;
		Error(msg.str());
	}

	return delta_f;
}

void Floq2SpinsCoreBehaviour::calc_ci(AllData* ad, int tr_id) const
{
	stringstream msg;
	msg << "Error: need to add functionality" << endl;
	Error(msg.str());
}

void Floq2SpinsCoreBehaviour::dump_std(AllData* ad) const
{
	RunParam* rp = ad->rp;
	ConfigParam* cp = ad->cp;
	ExpData* ed = ad->ed;

	int dump_obs = int(cp->params.find("dump_obs")->second);

	if (dump_obs == 1)
	{
		int num_trajectories = cp->num_trajectories;

		double* norm = ed->norm;
		MKL_Complex16* spec = ed->spec;

		string fn;

		fn = rp->path + "norm" + cp->fn_suffix;
		save_double_data(fn, norm, num_trajectories, 16, false);

		fn = rp->path + "spec" + cp->fn_suffix;
		save_complex_data(fn, spec, num_trajectories, 16, false);
	}
}

void Floq2SpinsCoreBehaviour::dump_lpn(AllData* ad) const
{
	RunParam* rp = ad->rp;
	ConfigParam* cp = ad->cp;
	ExpData* ed = ad->ed;

	int num_trajectories = cp->num_trajectories;

	int dump_obs = int(cp->params.find("dump_obs")->second);

	string fn;

	double* lambda = ed->lambda_now;
	fn = rp->path + "lambda" + cp->fn_suffix;
	save_double_data(fn, lambda, num_trajectories, 16, false);

	if (dump_obs == 1)
	{
		int* num_renorms = ed->num_renorms;
		fn = rp->path + "num_renorms" + cp->fn_suffix;
		save_int_data(fn, num_renorms, num_trajectories, false);

		MKL_Complex16* spec_lpn = ed->spec_lpn;
		fn = rp->path + "spec_lpn" + cp->fn_suffix;
		save_complex_data(fn, spec_lpn, num_trajectories, 16, false);

		int save_lambdas = int(ad->cp->params.find("save_lambdas")->second);
		if (save_lambdas > 0)
		{
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				fn = rp->path + "delta_f_by_delta_s_" + to_string(tr_id) + cp->fn_suffix;
				save_double_vector(fn, ed->lambdas[tr_id], 16, false);

				fn = rp->path + "deltas_s_" + to_string(tr_id) + cp->fn_suffix;
				save_double_vector(fn, ed->deltas_s[tr_id], 16, false);

				fn = rp->path + "deltas_f_" + to_string(tr_id) + cp->fn_suffix;
				save_double_vector(fn, ed->deltas_f[tr_id], 16, false);
			}
		}
	}
}

void Floq2SpinsCoreBehaviour::dump_std_evo(AllData* ad) const
{
	RunParam* rp = ad->rp;
	ConfigParam* cp = ad->cp;
	ExpData* ed = ad->ed;

	int num_trajectories = cp->num_trajectories;
	int dump_num_total = ed->dump_num_total;

	int dump_obs = int(cp->params.find("dump_obs")->second);
	int jump = int(cp->params.find("jump")->second);

	if (dump_obs == 1)
	{
		int* dump_periods = ed->dump_periods;

		double* norm_evo = ed->norm_evo;
		MKL_Complex16* spec_evo = ed->spec_evo;

		string fn;

		fn = rp->path + "periods" + cp->fn_suffix;
		save_int_data(fn, dump_periods, dump_num_total, false);

		fn = rp->path + "norm_evo" + cp->fn_suffix;
		save_2d_inv_double_data(fn, norm_evo, dump_num_total, num_trajectories, 16, false);

		fn = rp->path + "spec_evo" + cp->fn_suffix;
		save_2d_inv_complex_data(fn, spec_evo, dump_num_total, num_trajectories, 16, false);

		int num_random_obs = int(cp->params.find("num_random_obs")->second);
		int random_obs_seed = int(cp->params.find("random_obs_seed")->second);
		for (int obs_id = 0; obs_id < num_random_obs; obs_id++)
		{
			fn = rp->path + "random_obs_evo_" + to_string(random_obs_seed + obs_id) + cp->fn_suffix;
			std::vector<std::vector<std::complex<double>>> random_obs_evo_to_dump;
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				random_obs_evo_to_dump.push_back(ed->random_obs_evo[tr_id][obs_id]);
			}
			save_2d_vector(fn, random_obs_evo_to_dump, 16);
		}

		if (jump > 0)
		{
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				fn = rp->path + "jump_times_" + to_string(tr_id) + cp->fn_suffix;
				save_double_vector(fn, ed->jump_times[tr_id], 16, false);

				fn = rp->path + "diss_types_" + to_string(tr_id) + cp->fn_suffix;
				save_int_vector(fn, ed->diss_types[tr_id], false);

				if (jump > 1)
				{
					fn = rp->path + "jump_norms_" + to_string(tr_id) + cp->fn_suffix;
					save_double_vector(fn, ed->jump_norms[tr_id], 16, false);

					fn = rp->path + "jump_etas_" + to_string(tr_id) + cp->fn_suffix;
					save_double_vector(fn, ed->jump_etas[tr_id], 16, false);
				}
			}
		}
	}
}

void Floq2SpinsCoreBehaviour::dump_lpn_evo(AllData* ad) const
{
	RunParam* rp = ad->rp;
	ConfigParam* cp = ad->cp;
	ExpData* ed = ad->ed;

	int num_trajectories = cp->num_trajectories;
	int dump_num_total = ed->dump_num_total;

	int dump_obs = int(cp->params.find("dump_obs")->second);

	if (dump_obs == 1)
	{
		double* lambda_evo = ed->lambda_evo;
		MKL_Complex16* spec_lpn_evo = ed->spec_lpn_evo;

		string fn;

		fn = rp->path + "lambda_evo" + cp->fn_suffix;
		save_2d_inv_double_data(fn, lambda_evo, dump_num_total, num_trajectories, 16, false);

		fn = rp->path + "spec_lpn_evo" + cp->fn_suffix;
		save_2d_inv_complex_data(fn, spec_lpn_evo, dump_num_total, num_trajectories, 16, false);

		int num_random_obs = int(cp->params.find("num_random_obs")->second);
		int random_obs_seed = int(cp->params.find("random_obs_seed")->second);
		for (int obs_id = 0; obs_id < num_random_obs; obs_id++)
		{
			fn = rp->path + "random_obs_lpn_evo_" + to_string(random_obs_seed + obs_id) + cp->fn_suffix;
			std::vector<std::vector<std::complex<double>>> random_obs_evo_to_dump;
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				random_obs_evo_to_dump.push_back(ed->random_obs_lpn_evo[tr_id][obs_id]);
			}
			save_2d_vector(fn, random_obs_evo_to_dump, 16);
		}
	}
}





void FloqManySpinsCoreBehaviour::init_splits(AllData* ad) const
{
	RunParam* rp = ad->rp;
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;

	int num_threads = rp->num_threads;

	md->structure = init_split_structure_floqnspins(ad);
	md->splits = new Split[num_threads];
	for (int th_id = 0; th_id < num_threads; th_id++)
	{
		copy_struct_not_member(md->structure, &(md->splits)[th_id]);
	}

	cout << "Split initialization complete" << endl;
}

void FloqManySpinsCoreBehaviour::free_splits(AllData* ad) const
{
	RunParam* rp = ad->rp;
	MainData* md = ad->md;

	int num_threads = rp->num_threads;

	for (int i = 0; i < num_threads; i++)
	{
		delete_split_struct_not_member(&(md->splits[i]));
	}

	delete(md->splits);
	delete_split_struct(md->structure);
}

void FloqManySpinsCoreBehaviour::init_splits_deep(AllData* ad) const
{
	RunParam* rp = ad->rp;
	MainData* md = ad->md;

	int num_branches = md->num_ham_qj;
	int num_threads = rp->num_threads;

	int num_total = num_threads * num_branches;

	md->structure = init_split_structure_floqnspins_deep(ad);
	md->splits = new Split[num_total];

	for (int b_id = 0; b_id < num_branches; b_id++)
	{
		for (int th_id = 0; th_id < num_threads; th_id++)
		{
			int index = b_id * num_threads + th_id;
			copy_struct_not_member(&(md->structure)[b_id], &(md->splits)[index]);
		}
	}
}

void FloqManySpinsCoreBehaviour::free_splits_deep(AllData* ad) const
{
	RunParam* rp = ad->rp;
	MainData* md = ad->md;

	int num_branches = md->num_ham_qj;
	int num_threads = rp->num_threads;

	for (int b_id = 0; b_id < num_branches; b_id++)
	{
		for (int th_id = 0; th_id < num_threads; th_id++)
		{
			int index = b_id * num_threads + th_id;
			delete_split_struct_not_member(&(md->splits[index]));
		}
	}

	delete(md->splits);
	delete_split_struct(md->structure);
}

void FloqManySpinsCoreBehaviour::ex_period(AllData* ad, int tr_id, int th_id, int period_id) const
{
	MainData* md = ad->md;

	int sys_size = md->sys_size;
	Split* head = &(md->splits)[th_id];

	for (int b_id = 0; b_id < head->counter; b_id++)
	{
		one_period_branch(ad, head, tr_id, &(head->next)[b_id]);
	}
}

void FloqManySpinsCoreBehaviour::ex_period_trp_deep(AllData* ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;

	int num_branches = md->num_ham_qj;
	int num_sub_steps = int(cp->params.find("deep_num_steps")->second);

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		for (int sub_step_id = 0; sub_step_id < num_sub_steps; sub_step_id++)
		{
			one_sub_period_deep(ad, tr_id, part_id, th_id);
		}
	}
}

void FloqManySpinsCoreBehaviour::ex_period_obs_deep(AllData* ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

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

			one_sub_period_deep(ad, tr_id, part_id, th_id);
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

void FloqManySpinsCoreBehaviour::ex_period_obs_deep_lpn(AllData* ad, int period_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);

	int num_trajectories = cp->num_trajectories;

	int num_branches = md->num_ham_qj;
	int num_sub_steps_per_part = int(cp->params.find("deep_num_steps")->second);
	int num_sub_steps = num_branches * int(cp->params.find("deep_num_steps")->second);

	int dump_point_id = 0;
	int global_point_id = 0;

	int step_id = 0;

	CoreBehavior* tmp = new FloqManySpinsCoreBehaviour;

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		for (int sub_step_id = 0; sub_step_id < num_sub_steps_per_part; sub_step_id++)
		{
			step_id = part_id * num_sub_steps_per_part + sub_step_id;

			global_point_id = period_id * num_sub_steps + dump_point_id;
			dump_point_id++;
			int dump_id = global_point_id + 1;

			ed->curr_time = double(dump_id) / double(num_sub_steps) * md->T;

			one_sub_period_deep(ad, 0, part_id, 0);
			calc_chars_std(ad, 0);
			calc_chars_lpn(ad, 0, 0);
			evo_chars_std(ad, 0, dump_id);
			evo_chars_lpn(ad, 0, dump_id);

#pragma omp parallel for
			for (int tr_id = 1; tr_id < num_trajectories; tr_id++)
			{
				int thread_id = omp_get_thread_num();
				one_sub_period_deep(ad, tr_id, part_id, thread_id);
				calc_chars_std(ad, tr_id);
				lambda_lpn(ad, tmp, tr_id, 0);
				evo_chars_std(ad, tr_id, dump_id);
				evo_chars_lpn(ad, tr_id, dump_id);
			}

#pragma omp parallel for
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				if (dump_evo_sep == 1)
				{
					dump_adr_single(ad, tr_id, true);
				}
			}
		}
	}

	delete tmp;
}

void FloqManySpinsCoreBehaviour::ex_period_obs_deep_mult_lpn(AllData* ad, int period_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);

	int num_trajectories = cp->num_trajectories;

	int num_branches = md->num_ham_qj;
	int num_sub_steps_per_part = int(cp->params.find("deep_num_steps")->second);
	int num_sub_steps = num_branches * int(cp->params.find("deep_num_steps")->second);

	int lambda_per_periods = int(cp->params.find("lambda_per_periods")->second);

	int dump_point_id = 0;
	int global_point_id = 0;

	int step_id = 0;

	CoreBehavior* tmp = new FloqManySpinsCoreBehaviour;

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		for (int sub_step_id = 0; sub_step_id < num_sub_steps_per_part; sub_step_id++)
		{
			step_id = part_id * num_sub_steps_per_part + sub_step_id;

			global_point_id = period_id * num_sub_steps + dump_point_id;
			dump_point_id++;
			int dump_id = global_point_id + 1;

			ed->curr_time = double(dump_id) / double(num_sub_steps) * md->T;

#pragma omp parallel for
			for (int tr_id = 0; tr_id < num_trajectories / 2; tr_id++)
			{
				int thread_id = omp_get_thread_num();
				one_sub_period_deep(ad, tr_id, part_id, thread_id);
				calc_chars_std(ad, tr_id);
				calc_chars_lpn(ad, tr_id, tr_id);
				evo_chars_std(ad, tr_id, dump_id);
				evo_chars_lpn(ad, tr_id, dump_id);
			}

#pragma omp parallel for
			for (int tr_id = num_trajectories / 2; tr_id < num_trajectories; tr_id++)
			{
				int thread_id = omp_get_thread_num();
				one_sub_period_deep(ad, tr_id, part_id, thread_id);
				calc_chars_std(ad, tr_id);
				if (lambda_per_periods == 0)
				{
					lambda_lpn(ad, tmp, tr_id, tr_id - num_trajectories / 2);
				}
				if (lambda_per_periods == 1)
				{
					if ((part_id == num_branches - 1) && (sub_step_id == num_sub_steps_per_part - 1))
					{
						lambda_lpn_now(ad, tmp, tr_id, tr_id - num_trajectories / 2);
					}
				}
				evo_chars_std(ad, tr_id, dump_id);
				evo_chars_lpn(ad, tr_id, dump_id);
			}

#pragma omp parallel for
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				if (dump_evo_sep == 1)
				{
					dump_adr_single(ad, tr_id, true);
				}
			}
		}
	}

	delete tmp;
}

void FloqManySpinsCoreBehaviour::ex_period_obs_deep_lpn_per_period(struct AllData* ad, int period_id, int num_periods) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);

	int num_trajectories = cp->num_trajectories;

	int num_branches = md->num_ham_qj;
	int num_sub_steps_per_part = int(cp->params.find("deep_num_steps")->second);
	int num_sub_steps = num_branches * int(cp->params.find("deep_num_steps")->second);

	int dump_point_id = 0;
	int global_point_id = 0;

	int step_id = 0;

	CoreBehavior* tmp = new FloqManySpinsCoreBehaviour;

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		for (int sub_step_id = 0; sub_step_id < num_sub_steps_per_part; sub_step_id++)
		{
			step_id = part_id * num_sub_steps_per_part + sub_step_id;

			global_point_id = period_id * num_sub_steps + dump_point_id;
			dump_point_id++;
			int dump_id = global_point_id + 1;

			ed->curr_time = double(dump_id) / double(num_sub_steps) * md->T;

			one_sub_period_deep(ad, 0, part_id, 0);
			calc_chars_std(ad, 0);
			calc_chars_lpn(ad, 0, 0);
			evo_chars_std(ad, 0, dump_id);
			evo_chars_lpn(ad, 0, dump_id);

#pragma omp parallel for
			for (int tr_id = 1; tr_id < num_trajectories; tr_id++)
			{
				int thread_id = omp_get_thread_num();
				one_sub_period_deep(ad, tr_id, part_id, thread_id);
				calc_chars_std(ad, tr_id);
				lambda_lpn_per_periods(ad, tmp, tr_id, 0, num_sub_steps, global_point_id, num_periods);
				evo_chars_std(ad, tr_id, dump_id);
				evo_chars_lpn(ad, tr_id, dump_id);
			}

#pragma omp parallel for
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				if (dump_evo_sep == 1)
				{
					dump_adr_single(ad, tr_id, true);
				}
			}
		}
	}

	delete tmp;
}

void FloqManySpinsCoreBehaviour::ex_period_obs_deep_cd(AllData* ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

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

			one_sub_period_deep(ad, tr_id, part_id, th_id);
			calc_chars_std(ad, tr_id);
		}
	}
}

void FloqManySpinsCoreBehaviour::ex_period_obs_deep_sigma(AllData* ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;

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

			one_sub_period_deep(ad, tr_id, part_id, th_id);
			calc_chars_std(ad, tr_id);
		}
	}
}

void FloqManySpinsCoreBehaviour::rk_period(AllData* ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	for (int step_id = 0; step_id < cp->rk_ns; step_id++)
	{
		ed->times_all[tr_id] = double(period_id) * md->T + double(step_id) * ed->rk_step;
		rk_step_floqnspins(ad, tr_id, th_id, ed->rk_step);
	}
}

void FloqManySpinsCoreBehaviour::rk_period_trp_deep(AllData* ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int num_branches = md->num_ham_qj;
	int num_sub_steps_per_part = int(cp->params.find("deep_num_steps")->second);

	int step_id = 0;

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		for (int sub_step_id = 0; sub_step_id < num_sub_steps_per_part; sub_step_id++)
		{
			step_id = part_id * num_sub_steps_per_part + sub_step_id;

			for (int in_step_id = 0; in_step_id < cp->rk_ns; in_step_id++)
			{
				ed->times_all[tr_id] = double(period_id) * md->T + double(step_id * cp->rk_ns + in_step_id) * ed->rk_step;
				rk_step_floqnspins(ad, tr_id, th_id, ed->rk_step);
			}
		}
	}
}

void FloqManySpinsCoreBehaviour::rk_period_obs_deep(AllData* ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

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

			for (int in_step_id = 0; in_step_id < cp->rk_ns; in_step_id++)
			{
				ed->times_all[tr_id] = double(period_id) * md->T + double(step_id * cp->rk_ns + in_step_id) * ed->rk_step;
				rk_step_floqnspins(ad, tr_id, th_id, ed->rk_step);
			}

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

void FloqManySpinsCoreBehaviour::rk_period_obs_deep_lpn(AllData* ad, int period_id) const
{
	stringstream msg;
	msg << "Error: need to add functionality" << endl;
	Error(msg.str());
}

void FloqManySpinsCoreBehaviour::rk_period_obs_deep_cd(AllData* ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

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

			for (int in_step_id = 0; in_step_id < cp->rk_ns; in_step_id++)
			{
				ed->times_all[tr_id] = double(period_id) * md->T + double(step_id * cp->rk_ns + in_step_id) * ed->rk_step;
				rk_step_floqnspins(ad, tr_id, th_id, ed->rk_step);
			}

			calc_chars_std(ad, tr_id);
		}
	}
}

void FloqManySpinsCoreBehaviour::rk_period_obs_deep_sigma(AllData* ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

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

			dump_point_id++;

			for (int in_step_id = 0; in_step_id < cp->rk_ns; in_step_id++)
			{
				ed->times_all[tr_id] = double(period_id) * md->T + double(step_id * cp->rk_ns + in_step_id) * ed->rk_step;
				rk_step_floqnspins(ad, tr_id, th_id, ed->rk_step);
			}

			calc_chars_std(ad, tr_id);
		}
	}
}

void FloqManySpinsCoreBehaviour::calc_chars_std_start(AllData* ad, int tr_id) const
{
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int sys_size = md->sys_size;
	MKL_Complex16* phi = &(ed->phi_all[tr_id * sys_size]);
	double* adr = &(ed->abs_diag_rho_all[tr_id * sys_size]);

	double norm = norm_square(phi, sys_size);

	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		adr[st_id] = mult_scalar_double(mult_scalar_complex(&phi[st_id], &phi[st_id], 1), 1.0 / norm).real;
	}

	MKL_Complex16 spec = get_spec_flqnsp(ad, tr_id);

	ed->norm[tr_id] = norm;
	ed->spec[tr_id] = spec;

	calc_random_obs(ad, tr_id);
}

void FloqManySpinsCoreBehaviour::calc_chars_std(AllData* ad, int tr_id) const
{
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int sys_size = md->sys_size;
	MKL_Complex16* phi = &(ed->phi_all[tr_id * sys_size]);
	double* adr = &(ed->abs_diag_rho_all[tr_id * sys_size]);

	double norm = norm_square(phi, sys_size);

	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		adr[st_id] = mult_scalar_double(mult_scalar_complex(&phi[st_id], &phi[st_id], 1), 1.0 / norm).real;
	}

	MKL_Complex16 spec = get_spec_flqnsp(ad, tr_id);

	ed->norm[tr_id] = norm;
	ed->spec[tr_id] = spec;

	calc_random_obs(ad, tr_id);
}

void FloqManySpinsCoreBehaviour::calc_chars_lpn_start(AllData* ad, int tr_id, int base_tr_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	double lambda = 0.0;
	double lambda_now = 0.0;
	MKL_Complex16 spec_lpn = ed->spec[tr_id];

	double delta_s = this->calc_delta_f(ad, tr_id, base_tr_id); // Important! Here we use calc_delta_f not calc_delta_s

	ed->lambda[tr_id] = lambda;
	ed->lambda_now[tr_id] = lambda_now;
	ed->delta_s[tr_id] = delta_s;

	ed->spec_lpn[tr_id] = spec_lpn;

	calc_random_obs_lpn_start(ad, tr_id);
}

void FloqManySpinsCoreBehaviour::calc_chars_lpn(AllData* ad, int tr_id, int base_tr_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int sys_size = md->sys_size;
	MKL_Complex16* phi = &(ed->phi_all[tr_id * sys_size]);
	double* adr = &(ed->abs_diag_rho_all[tr_id * sys_size]);
	double norm_2 = norm_square(phi, sys_size);
	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		adr[st_id] = mult_scalar_double(mult_scalar_complex(&phi[st_id], &phi[st_id], 1), 1.0 / norm_2).real;
	}

	MKL_Complex16 spec_lpn = get_spec_flqnsp(ad, tr_id);

	ed->spec_lpn[tr_id] = spec_lpn;

	calc_random_obs_lpn(ad, tr_id);
}

double FloqManySpinsCoreBehaviour::calc_T(AllData* ad) const
{
	MainData* md = ad->md;
	double T = md->T;
	return T;
}

void FloqManySpinsCoreBehaviour::evo_chars_std(AllData* ad, int tr_id, int dump_id) const
{
	ConfigParam* cp = ad->cp;
	ExpData* ed = ad->ed;

	int dump_num_total = ed->dump_num_total;

	int index = tr_id * dump_num_total + dump_id;

	ed->norm_evo[index] = ed->norm[tr_id];
	ed->spec_evo[index] = ed->spec[tr_id];

	int num_random_obs = int(cp->params.find("num_random_obs")->second);
	for (int obs_id = 0; obs_id < num_random_obs; obs_id++)
	{
		ed->random_obs_evo[tr_id][obs_id].push_back(ed->random_obs[tr_id][obs_id]);
	}
}

void FloqManySpinsCoreBehaviour::evo_chars_lpn(AllData* ad, int tr_id, int dump_id) const
{
	ConfigParam* cp = ad->cp;
	ExpData* ed = ad->ed;

	int dump_num_total = ed->dump_num_total;

	int index = tr_id * dump_num_total + dump_id;

	ed->lambda_evo[index] = ed->lambda_now[tr_id];
	ed->spec_lpn_evo[index] = ed->spec_lpn[tr_id];

	int num_random_obs = int(cp->params.find("num_random_obs")->second);
	for (int obs_id = 0; obs_id < num_random_obs; obs_id++)
	{
		ed->random_obs_lpn_evo[tr_id][obs_id].push_back(ed->random_obs_lpn[tr_id][obs_id]);
	}
}

double FloqManySpinsCoreBehaviour::calc_delta_s(AllData* ad, int tr_id, int base_tr_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int sys_size = md->sys_size;

	int lpn_type = int(cp->params.find("lpn_type")->second);
	int num_random_obs = int(ad->cp->params.find("num_random_obs")->second);

	double delta_s = 0.0;
	if (lpn_type == -1)
	{
		MKL_Complex16 base = ed->spec_lpn[base_tr_id];
		MKL_Complex16 var = ed->spec_lpn[tr_id];
		double tmp = pow((base.real - var.real), 2) + pow((base.imag - var.imag), 2);
		delta_s = sqrt(tmp);
	}
	else if (lpn_type >= 0 && lpn_type < num_random_obs)
	{
		complex<double> base = ed->random_obs_lpn[base_tr_id][lpn_type];
		complex<double> var = ed->random_obs_lpn[tr_id][lpn_type];
		double tmp = pow((base.real() - var.real()), 2) + pow((base.imag() - var.imag()), 2);
		delta_s = sqrt(tmp);
	}
	else
	{
		stringstream msg;
		msg << "Error: Wrong lpn_type: " << lpn_type << endl;
		Error(msg.str());
	}

	return delta_s;
}

double FloqManySpinsCoreBehaviour::calc_delta_f(AllData* ad, int tr_id, int base_tr_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int sys_size = md->sys_size;

	int lpn_type = int(cp->params.find("lpn_type")->second);
	int num_random_obs = int(ad->cp->params.find("num_random_obs")->second);

	double delta_f = 0.0;
	if (lpn_type == -1)
	{
		MKL_Complex16 base = ed->spec[base_tr_id];
		MKL_Complex16 var = ed->spec[tr_id];
		double tmp = pow((base.real - var.real), 2) + pow((base.imag - var.imag), 2);
		delta_f = sqrt(tmp);
	}
	else if (lpn_type >= 0 && lpn_type < num_random_obs)
	{
		complex<double> base = ed->random_obs[base_tr_id][lpn_type];
		complex<double> var = ed->random_obs[tr_id][lpn_type];
		double tmp = pow((base.real() - var.real()), 2) + pow((base.imag() - var.imag()), 2);
		delta_f = sqrt(tmp);
	}
	else
	{
		stringstream msg;
		msg << "Error: Wrong lpn_type: " << lpn_type << endl;
		Error(msg.str());
	}

	return delta_f;
}

void FloqManySpinsCoreBehaviour::calc_ci(AllData* ad, int tr_id) const
{
	stringstream msg;
	msg << "Error: need to add functionality" << endl;
	Error(msg.str());
}

void FloqManySpinsCoreBehaviour::dump_std(AllData* ad) const
{
	RunParam* rp = ad->rp;
	ConfigParam* cp = ad->cp;
	ExpData* ed = ad->ed;

	int dump_obs = int(cp->params.find("dump_obs")->second);

	if (dump_obs == 1)
	{
		int num_trajectories = cp->num_trajectories;

		double* norm = ed->norm;
		MKL_Complex16* spec = ed->spec;

		string fn;

		fn = rp->path + "norm" + cp->fn_suffix;
		save_double_data(fn, norm, num_trajectories, 16, false);

		fn = rp->path + "spec" + cp->fn_suffix;
		save_complex_data(fn, spec, num_trajectories, 16, false);
	}
}

void FloqManySpinsCoreBehaviour::dump_lpn(AllData* ad) const
{
	RunParam* rp = ad->rp;
	ConfigParam* cp = ad->cp;
	ExpData* ed = ad->ed;

	int num_trajectories = cp->num_trajectories;

	int dump_obs = int(cp->params.find("dump_obs")->second);

	string fn;

	double* lambda = ed->lambda_now;
	fn = rp->path + "lambda" + cp->fn_suffix;
	save_double_data(fn, lambda, num_trajectories, 16, false);

	if (dump_obs == 1)
	{
		int* num_renorms = ed->num_renorms;
		fn = rp->path + "num_renorms" + cp->fn_suffix;
		save_int_data(fn, num_renorms, num_trajectories, false);

		MKL_Complex16* spec_lpn = ed->spec_lpn;
		fn = rp->path + "spec_lpn" + cp->fn_suffix;
		save_complex_data(fn, spec_lpn, num_trajectories, 16, false);

		int save_lambdas = int(ad->cp->params.find("save_lambdas")->second);
		if (save_lambdas > 0)
		{
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				fn = rp->path + "delta_f_by_delta_s_" + to_string(tr_id) + cp->fn_suffix;
				save_double_vector(fn, ed->lambdas[tr_id], 16, false);

				fn = rp->path + "deltas_s_" + to_string(tr_id) + cp->fn_suffix;
				save_double_vector(fn, ed->deltas_s[tr_id], 16, false);

				fn = rp->path + "deltas_f_" + to_string(tr_id) + cp->fn_suffix;
				save_double_vector(fn, ed->deltas_f[tr_id], 16, false);
			}
		}
	}
}

void FloqManySpinsCoreBehaviour::dump_std_evo(AllData* ad) const
{
	RunParam* rp = ad->rp;
	ConfigParam* cp = ad->cp;
	ExpData* ed = ad->ed;

	int num_trajectories = cp->num_trajectories;
	int dump_num_total = ed->dump_num_total;

	int dump_obs = int(cp->params.find("dump_obs")->second);
	int jump = int(cp->params.find("jump")->second);

	if (dump_obs == 1)
	{
		int* dump_periods = ed->dump_periods;

		double* norm_evo = ed->norm_evo;
		MKL_Complex16* spec_evo = ed->spec_evo;

		string fn;

		fn = rp->path + "periods" + cp->fn_suffix;
		save_int_data(fn, dump_periods, dump_num_total, false);

		fn = rp->path + "norm_evo" + cp->fn_suffix;
		save_2d_inv_double_data(fn, norm_evo, dump_num_total, num_trajectories, 16, false);

		fn = rp->path + "spec_evo" + cp->fn_suffix;
		save_2d_inv_complex_data(fn, spec_evo, dump_num_total, num_trajectories, 16, false);

		int num_random_obs = int(cp->params.find("num_random_obs")->second);
		int random_obs_seed = int(cp->params.find("random_obs_seed")->second);
		for (int obs_id = 0; obs_id < num_random_obs; obs_id++)
		{
			fn = rp->path + "random_obs_evo_" + to_string(random_obs_seed + obs_id) + cp->fn_suffix;
			std::vector<std::vector<std::complex<double>>> random_obs_evo_to_dump;
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				random_obs_evo_to_dump.push_back(ed->random_obs_evo[tr_id][obs_id]);
			}
			save_2d_vector(fn, random_obs_evo_to_dump, 16);
		}

		if (jump > 0)
		{
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				fn = rp->path + "jump_times_" + to_string(tr_id) + cp->fn_suffix;
				save_double_vector(fn, ed->jump_times[tr_id], 16, false);

				fn = rp->path + "diss_types_" + to_string(tr_id) + cp->fn_suffix;
				save_int_vector(fn, ed->diss_types[tr_id], false);

				if (jump > 1)
				{
					fn = rp->path + "jump_norms_" + to_string(tr_id) + cp->fn_suffix;
					save_double_vector(fn, ed->jump_norms[tr_id], 16, false);

					fn = rp->path + "jump_etas_" + to_string(tr_id) + cp->fn_suffix;
					save_double_vector(fn, ed->jump_etas[tr_id], 16, false);
				}
			}
		}
	}
}

void FloqManySpinsCoreBehaviour::dump_lpn_evo(AllData* ad) const
{
	RunParam* rp = ad->rp;
	ConfigParam* cp = ad->cp;
	ExpData* ed = ad->ed;

	int num_trajectories = cp->num_trajectories;
	int dump_num_total = ed->dump_num_total;

	int dump_obs = int(cp->params.find("dump_obs")->second);

	if (dump_obs == 1)
	{
		double* lambda_evo = ed->lambda_evo;
		MKL_Complex16* spec_lpn_evo = ed->spec_lpn_evo;

		string fn;

		fn = rp->path + "lambda_evo" + cp->fn_suffix;
		save_2d_inv_double_data(fn, lambda_evo, dump_num_total, num_trajectories, 16, false);

		fn = rp->path + "spec_lpn_evo" + cp->fn_suffix;
		save_2d_inv_complex_data(fn, spec_lpn_evo, dump_num_total, num_trajectories, 16, false);

		int num_random_obs = int(cp->params.find("num_random_obs")->second);
		int random_obs_seed = int(cp->params.find("random_obs_seed")->second);
		for (int obs_id = 0; obs_id < num_random_obs; obs_id++)
		{
			fn = rp->path + "random_obs_lpn_evo_" + to_string(random_obs_seed + obs_id) + cp->fn_suffix;
			std::vector<std::vector<std::complex<double>>> random_obs_evo_to_dump;
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				random_obs_evo_to_dump.push_back(ed->random_obs_lpn_evo[tr_id][obs_id]);
			}
			save_2d_vector(fn, random_obs_evo_to_dump, 16);
		}
	}
}



void FloqSpinPhCoreBehaviour::init_splits(AllData* ad) const
{
	RunParam* rp = ad->rp;
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;

	int num_threads = rp->num_threads;

	md->structure = init_split_structure_flq_sp_ph(ad);
	md->splits = new Split[num_threads];
	for (int th_id = 0; th_id < num_threads; th_id++)
	{
		copy_struct_not_member(md->structure, &(md->splits)[th_id]);
	}

	cout << "Split initialization complete" << endl;
}

void FloqSpinPhCoreBehaviour::free_splits(AllData* ad) const
{
	RunParam* rp = ad->rp;
	MainData* md = ad->md;

	int num_threads = rp->num_threads;

	for (int i = 0; i < num_threads; i++)
	{
		delete_split_struct_not_member(&(md->splits[i]));
	}

	delete(md->splits);
	delete_split_struct(md->structure);
}

void FloqSpinPhCoreBehaviour::init_splits_deep(AllData* ad) const
{
	RunParam* rp = ad->rp;
	MainData* md = ad->md;

	int num_branches = md->num_ham_qj;
	int num_threads = rp->num_threads;

	int num_total = num_threads * num_branches;

	md->structure = init_split_structure_flq_sp_ph_deep(ad);
	md->splits = new Split[num_total];

	for (int b_id = 0; b_id < num_branches; b_id++)
	{
		for (int th_id = 0; th_id < num_threads; th_id++)
		{
			int index = b_id * num_threads + th_id;
			copy_struct_not_member(&(md->structure)[b_id], &(md->splits)[index]);
		}
	}
}

void FloqSpinPhCoreBehaviour::free_splits_deep(AllData* ad) const
{
	RunParam* rp = ad->rp;
	MainData* md = ad->md;

	int num_branches = md->num_ham_qj;
	int num_threads = rp->num_threads;

	for (int b_id = 0; b_id < num_branches; b_id++)
	{
		for (int th_id = 0; th_id < num_threads; th_id++)
		{
			int index = b_id * num_threads + th_id;
			delete_split_struct_not_member(&(md->splits[index]));
		}
	}

	delete(md->splits);
	delete_split_struct(md->structure);
}

void FloqSpinPhCoreBehaviour::ex_period(AllData* ad, int tr_id, int th_id, int period_id) const
{
	MainData* md = ad->md;

	int sys_size = md->sys_size;
	Split* head = &(md->splits)[th_id];

	for (int b_id = 0; b_id < head->counter; b_id++)
	{
		one_period_branch(ad, head, tr_id, &(head->next)[b_id]);
	}
}

void FloqSpinPhCoreBehaviour::ex_period_trp_deep(AllData* ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;

	int num_branches = md->num_ham_qj;
	int num_sub_steps = int(cp->params.find("deep_num_steps")->second);

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		for (int sub_step_id = 0; sub_step_id < num_sub_steps; sub_step_id++)
		{
			one_sub_period_deep(ad, tr_id, part_id, th_id);
		}
	}
}

void FloqSpinPhCoreBehaviour::ex_period_obs_deep(AllData* ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

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

			one_sub_period_deep(ad, tr_id, part_id, th_id);
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

void FloqSpinPhCoreBehaviour::ex_period_obs_deep_lpn(AllData* ad, int period_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);

	int num_trajectories = cp->num_trajectories;

	int num_branches = md->num_ham_qj;
	int num_sub_steps_per_part = int(cp->params.find("deep_num_steps")->second);
	int num_sub_steps = num_branches * int(cp->params.find("deep_num_steps")->second);

	int dump_point_id = 0;
	int global_point_id = 0;

	int step_id = 0;

	CoreBehavior* tmp = new FloqSpinPhCoreBehaviour;

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		for (int sub_step_id = 0; sub_step_id < num_sub_steps_per_part; sub_step_id++)
		{
			step_id = part_id * num_sub_steps_per_part + sub_step_id;

			global_point_id = period_id * num_sub_steps + dump_point_id;
			dump_point_id++;
			int dump_id = global_point_id + 1;

			ed->curr_time = double(dump_id) / double(num_sub_steps) * md->T;

			one_sub_period_deep(ad, 0, part_id, 0);
			calc_chars_std(ad, 0);
			calc_chars_lpn(ad, 0, 0);
			evo_chars_std(ad, 0, dump_id);
			evo_chars_lpn(ad, 0, dump_id);

#pragma omp parallel for
			for (int tr_id = 1; tr_id < num_trajectories; tr_id++)
			{
				int thread_id = omp_get_thread_num();
				one_sub_period_deep(ad, tr_id, part_id, thread_id);
				calc_chars_std(ad, tr_id);
				lambda_lpn(ad, tmp, tr_id, 0);
				evo_chars_std(ad, tr_id, dump_id);
				evo_chars_lpn(ad, tr_id, dump_id);
			}

#pragma omp parallel for
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				if (dump_evo_sep == 1)
				{
					dump_adr_single(ad, tr_id, true);
				}
			}
		}
	}

	delete tmp;
}

void FloqSpinPhCoreBehaviour::ex_period_obs_deep_mult_lpn(AllData* ad, int period_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);

	int num_trajectories = cp->num_trajectories;

	int num_branches = md->num_ham_qj;
	int num_sub_steps_per_part = int(cp->params.find("deep_num_steps")->second);
	int num_sub_steps = num_branches * int(cp->params.find("deep_num_steps")->second);

	int lambda_per_periods = int(cp->params.find("lambda_per_periods")->second);

	int dump_point_id = 0;
	int global_point_id = 0;

	int step_id = 0;

	CoreBehavior* tmp = new FloqSpinPhCoreBehaviour;

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		for (int sub_step_id = 0; sub_step_id < num_sub_steps_per_part; sub_step_id++)
		{
			step_id = part_id * num_sub_steps_per_part + sub_step_id;

			global_point_id = period_id * num_sub_steps + dump_point_id;
			dump_point_id++;
			int dump_id = global_point_id + 1;

			ed->curr_time = double(dump_id) / double(num_sub_steps) * md->T;

#pragma omp parallel for
			for (int tr_id = 0; tr_id < num_trajectories / 2; tr_id++)
			{
				int thread_id = omp_get_thread_num();
				one_sub_period_deep(ad, tr_id, part_id, thread_id);
				calc_chars_std(ad, tr_id);
				calc_chars_lpn(ad, tr_id, tr_id);
				evo_chars_std(ad, tr_id, dump_id);
				evo_chars_lpn(ad, tr_id, dump_id);
			}

#pragma omp parallel for
			for (int tr_id = num_trajectories / 2; tr_id < num_trajectories; tr_id++)
			{
				int thread_id = omp_get_thread_num();
				one_sub_period_deep(ad, tr_id, part_id, thread_id);
				calc_chars_std(ad, tr_id);
				if (lambda_per_periods == 0)
				{
					lambda_lpn(ad, tmp, tr_id, tr_id - num_trajectories / 2);
				}
				if (lambda_per_periods == 1)
				{
					if ((part_id == num_branches - 1) && (sub_step_id == num_sub_steps_per_part - 1))
					{
						lambda_lpn_now(ad, tmp, tr_id, tr_id - num_trajectories / 2);
					}
				}
				evo_chars_std(ad, tr_id, dump_id);
				evo_chars_lpn(ad, tr_id, dump_id);
			}

#pragma omp parallel for
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				if (dump_evo_sep == 1)
				{
					dump_adr_single(ad, tr_id, true);
				}
			}
		}
	}

	delete tmp;
}

void FloqSpinPhCoreBehaviour::ex_period_obs_deep_lpn_per_period(struct AllData* ad, int period_id, int num_periods) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);

	int num_trajectories = cp->num_trajectories;

	int num_branches = md->num_ham_qj;
	int num_sub_steps_per_part = int(cp->params.find("deep_num_steps")->second);
	int num_sub_steps = num_branches * int(cp->params.find("deep_num_steps")->second);

	int dump_point_id = 0;
	int global_point_id = 0;

	int step_id = 0;

	CoreBehavior* tmp = new FloqSpinPhCoreBehaviour;

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		for (int sub_step_id = 0; sub_step_id < num_sub_steps_per_part; sub_step_id++)
		{
			step_id = part_id * num_sub_steps_per_part + sub_step_id;

			global_point_id = period_id * num_sub_steps + dump_point_id;
			dump_point_id++;
			int dump_id = global_point_id + 1;

			ed->curr_time = double(dump_id) / double(num_sub_steps) * md->T;

			one_sub_period_deep(ad, 0, part_id, 0);
			calc_chars_std(ad, 0);
			calc_chars_lpn(ad, 0, 0);
			evo_chars_std(ad, 0, dump_id);
			evo_chars_lpn(ad, 0, dump_id);

#pragma omp parallel for
			for (int tr_id = 1; tr_id < num_trajectories; tr_id++)
			{
				int thread_id = omp_get_thread_num();
				one_sub_period_deep(ad, tr_id, part_id, thread_id);
				calc_chars_std(ad, tr_id);
				lambda_lpn_per_periods(ad, tmp, tr_id, 0, num_sub_steps, global_point_id, num_periods);
				evo_chars_std(ad, tr_id, dump_id);
				evo_chars_lpn(ad, tr_id, dump_id);
			}

#pragma omp parallel for
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				if (dump_evo_sep == 1)
				{
					dump_adr_single(ad, tr_id, true);
				}
			}
		}
	}

	delete tmp;
}

void FloqSpinPhCoreBehaviour::ex_period_obs_deep_cd(AllData* ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

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

			one_sub_period_deep(ad, tr_id, part_id, th_id);
			calc_chars_std(ad, tr_id);
		}
	}
}

void FloqSpinPhCoreBehaviour::ex_period_obs_deep_sigma(AllData* ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;

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

			one_sub_period_deep(ad, tr_id, part_id, th_id);
			calc_chars_std(ad, tr_id);
		}
	}
}

void FloqSpinPhCoreBehaviour::rk_period(AllData* ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	for (int step_id = 0; step_id < cp->rk_ns; step_id++)
	{
		ed->times_all[tr_id] = double(period_id) * md->T + double(step_id) * ed->rk_step;
		rk_step_flq_sp_ph(ad, tr_id, th_id, ed->rk_step);
	}
}

void FloqSpinPhCoreBehaviour::rk_period_trp_deep(AllData* ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int num_branches = md->num_ham_qj;
	int num_sub_steps_per_part = int(cp->params.find("deep_num_steps")->second);

	int step_id = 0;

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		for (int sub_step_id = 0; sub_step_id < num_sub_steps_per_part; sub_step_id++)
		{
			step_id = part_id * num_sub_steps_per_part + sub_step_id;

			for (int in_step_id = 0; in_step_id < cp->rk_ns; in_step_id++)
			{
				ed->times_all[tr_id] = double(period_id) * md->T + double(step_id * cp->rk_ns + in_step_id) * ed->rk_step;
				rk_step_flq_sp_ph(ad, tr_id, th_id, ed->rk_step);
			}
		}
	}
}

void FloqSpinPhCoreBehaviour::rk_period_obs_deep(AllData* ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

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

			for (int in_step_id = 0; in_step_id < cp->rk_ns; in_step_id++)
			{
				ed->times_all[tr_id] = double(period_id) * md->T + double(step_id * cp->rk_ns + in_step_id) * ed->rk_step;
				rk_step_flq_sp_ph(ad, tr_id, th_id, ed->rk_step);
			}

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

void FloqSpinPhCoreBehaviour::rk_period_obs_deep_lpn(AllData* ad, int period_id) const
{
	stringstream msg;
	msg << "Error: need to add functionality" << endl;
	Error(msg.str());
}

void FloqSpinPhCoreBehaviour::rk_period_obs_deep_cd(AllData* ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

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

			for (int in_step_id = 0; in_step_id < cp->rk_ns; in_step_id++)
			{
				ed->times_all[tr_id] = double(period_id) * md->T + double(step_id * cp->rk_ns + in_step_id) * ed->rk_step;
				rk_step_flq_sp_ph(ad, tr_id, th_id, ed->rk_step);
			}

			calc_chars_std(ad, tr_id);
		}
	}
}

void FloqSpinPhCoreBehaviour::rk_period_obs_deep_sigma(AllData* ad, int tr_id, int th_id, int period_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

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

			dump_point_id++;

			for (int in_step_id = 0; in_step_id < cp->rk_ns; in_step_id++)
			{
				ed->times_all[tr_id] = double(period_id) * md->T + double(step_id * cp->rk_ns + in_step_id) * ed->rk_step;
				rk_step_flq_sp_ph(ad, tr_id, th_id, ed->rk_step);
			}

			calc_chars_std(ad, tr_id);
		}
	}
}

void FloqSpinPhCoreBehaviour::calc_chars_std_start(AllData* ad, int tr_id) const
{
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int sys_size = md->sys_size;
	MKL_Complex16* phi = &(ed->phi_all[tr_id * sys_size]);
	double* adr = &(ed->abs_diag_rho_all[tr_id * sys_size]);

	double norm = norm_square(phi, sys_size);

	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		adr[st_id] = mult_scalar_double(mult_scalar_complex(&phi[st_id], &phi[st_id], 1), 1.0 / norm).real;
	}

	MKL_Complex16 spec = get_spec_flq_sp_ph(ad, tr_id);

	ed->norm[tr_id] = norm;
	ed->spec[tr_id] = spec;

	calc_random_obs(ad, tr_id);
}

void FloqSpinPhCoreBehaviour::calc_chars_std(AllData* ad, int tr_id) const
{
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int sys_size = md->sys_size;
	MKL_Complex16* phi = &(ed->phi_all[tr_id * sys_size]);
	double* adr = &(ed->abs_diag_rho_all[tr_id * sys_size]);

	double norm = norm_square(phi, sys_size);

	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		adr[st_id] = mult_scalar_double(mult_scalar_complex(&phi[st_id], &phi[st_id], 1), 1.0 / norm).real;
	}

	MKL_Complex16 spec = get_spec_flq_sp_ph(ad, tr_id);

	ed->norm[tr_id] = norm;
	ed->spec[tr_id] = spec;

	calc_random_obs(ad, tr_id);
}

void FloqSpinPhCoreBehaviour::calc_chars_lpn_start(AllData* ad, int tr_id, int base_tr_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	double lambda = 0.0;
	double lambda_now = 0.0;
	MKL_Complex16 spec_lpn = ed->spec[tr_id];

	double delta_s = this->calc_delta_f(ad, tr_id, base_tr_id); // Important! Here we use calc_delta_f not calc_delta_s

	ed->lambda[tr_id] = lambda;
	ed->lambda_now[tr_id] = lambda_now;
	ed->delta_s[tr_id] = delta_s;

	ed->spec_lpn[tr_id] = spec_lpn;

	calc_random_obs_lpn_start(ad, tr_id);
}

void FloqSpinPhCoreBehaviour::calc_chars_lpn(AllData* ad, int tr_id, int base_tr_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int sys_size = md->sys_size;
	MKL_Complex16* phi = &(ed->phi_all[tr_id * sys_size]);
	double* adr = &(ed->abs_diag_rho_all[tr_id * sys_size]);
	double norm_2 = norm_square(phi, sys_size);
	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		adr[st_id] = mult_scalar_double(mult_scalar_complex(&phi[st_id], &phi[st_id], 1), 1.0 / norm_2).real;
	}

	MKL_Complex16 spec_lpn = get_spec_flq_sp_ph(ad, tr_id);

	ed->spec_lpn[tr_id] = spec_lpn;

	calc_random_obs_lpn(ad, tr_id);
}

double FloqSpinPhCoreBehaviour::calc_T(AllData* ad) const
{
	MainData* md = ad->md;
	double T = md->T;
	return T;
}

void FloqSpinPhCoreBehaviour::evo_chars_std(AllData* ad, int tr_id, int dump_id) const
{
	ConfigParam* cp = ad->cp;
	ExpData* ed = ad->ed;

	int dump_num_total = ed->dump_num_total;

	int index = tr_id * dump_num_total + dump_id;

	ed->norm_evo[index] = ed->norm[tr_id];
	ed->spec_evo[index] = ed->spec[tr_id];

	int num_random_obs = int(cp->params.find("num_random_obs")->second);
	for (int obs_id = 0; obs_id < num_random_obs; obs_id++)
	{
		ed->random_obs_evo[tr_id][obs_id].push_back(ed->random_obs[tr_id][obs_id]);
	}
}

void FloqSpinPhCoreBehaviour::evo_chars_lpn(AllData* ad, int tr_id, int dump_id) const
{
	ConfigParam* cp = ad->cp;
	ExpData* ed = ad->ed;

	int dump_num_total = ed->dump_num_total;

	int index = tr_id * dump_num_total + dump_id;

	ed->lambda_evo[index] = ed->lambda_now[tr_id];
	ed->spec_lpn_evo[index] = ed->spec_lpn[tr_id];

	int num_random_obs = int(cp->params.find("num_random_obs")->second);
	for (int obs_id = 0; obs_id < num_random_obs; obs_id++)
	{
		ed->random_obs_lpn_evo[tr_id][obs_id].push_back(ed->random_obs_lpn[tr_id][obs_id]);
	}
}

double FloqSpinPhCoreBehaviour::calc_delta_s(AllData* ad, int tr_id, int base_tr_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int sys_size = md->sys_size;

	int lpn_type = int(cp->params.find("lpn_type")->second);
	int num_random_obs = int(ad->cp->params.find("num_random_obs")->second);

	double delta_s = 0.0;
	if (lpn_type == -1)
	{
		MKL_Complex16 base = ed->spec_lpn[base_tr_id];
		MKL_Complex16 var = ed->spec_lpn[tr_id];
		double tmp = pow((base.real - var.real), 2) + pow((base.imag - var.imag), 2);
		delta_s = sqrt(tmp);
	}
	else if (lpn_type >= 0 && lpn_type < num_random_obs)
	{
		complex<double> base = ed->random_obs_lpn[base_tr_id][lpn_type];
		complex<double> var = ed->random_obs_lpn[tr_id][lpn_type];
		double tmp = pow((base.real() - var.real()), 2) + pow((base.imag() - var.imag()), 2);
		delta_s = sqrt(tmp);
	}
	else
	{
		stringstream msg;
		msg << "Error: Wrong lpn_type: " << lpn_type << endl;
		Error(msg.str());
	}

	return delta_s;
}

double FloqSpinPhCoreBehaviour::calc_delta_f(AllData* ad, int tr_id, int base_tr_id) const
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int sys_size = md->sys_size;

	int lpn_type = int(cp->params.find("lpn_type")->second);
	int num_random_obs = int(ad->cp->params.find("num_random_obs")->second);

	double delta_f = 0.0;
	if (lpn_type == -1)
	{
		MKL_Complex16 base = ed->spec[base_tr_id];
		MKL_Complex16 var = ed->spec[tr_id];
		double tmp = pow((base.real - var.real), 2) + pow((base.imag - var.imag), 2);
		delta_f = sqrt(tmp);
	}
	else if (lpn_type >= 0 && lpn_type < num_random_obs)
	{
		complex<double> base = ed->random_obs[base_tr_id][lpn_type];
		complex<double> var = ed->random_obs[tr_id][lpn_type];
		double tmp = pow((base.real() - var.real()), 2) + pow((base.imag() - var.imag()), 2);
		delta_f = sqrt(tmp);
	}
	else
	{
		stringstream msg;
		msg << "Error: Wrong lpn_type: " << lpn_type << endl;
		Error(msg.str());
	}

	return delta_f;
}

void FloqSpinPhCoreBehaviour::calc_ci(AllData* ad, int tr_id) const
{
	stringstream msg;
	msg << "Error: need to add functionality" << endl;
	Error(msg.str());
}

void FloqSpinPhCoreBehaviour::dump_std(AllData* ad) const
{
	RunParam* rp = ad->rp;
	ConfigParam* cp = ad->cp;
	ExpData* ed = ad->ed;

	int dump_obs = int(cp->params.find("dump_obs")->second);

	if (dump_obs == 1)
	{
		int num_trajectories = cp->num_trajectories;

		double* norm = ed->norm;
		MKL_Complex16* spec = ed->spec;

		string fn;

		fn = rp->path + "norm" + cp->fn_suffix;
		save_double_data(fn, norm, num_trajectories, 16, false);

		fn = rp->path + "spec" + cp->fn_suffix;
		save_complex_data(fn, spec, num_trajectories, 16, false);
	}
}

void FloqSpinPhCoreBehaviour::dump_lpn(AllData* ad) const
{
	RunParam* rp = ad->rp;
	ConfigParam* cp = ad->cp;
	ExpData* ed = ad->ed;

	int num_trajectories = cp->num_trajectories;

	int dump_obs = int(cp->params.find("dump_obs")->second);

	string fn;

	double* lambda = ed->lambda_now;
	fn = rp->path + "lambda" + cp->fn_suffix;
	save_double_data(fn, lambda, num_trajectories, 16, false);

	if (dump_obs == 1)
	{
		int* num_renorms = ed->num_renorms;
		fn = rp->path + "num_renorms" + cp->fn_suffix;
		save_int_data(fn, num_renorms, num_trajectories, false);

		MKL_Complex16* spec_lpn = ed->spec_lpn;
		fn = rp->path + "spec_lpn" + cp->fn_suffix;
		save_complex_data(fn, spec_lpn, num_trajectories, 16, false);

		int save_lambdas = int(ad->cp->params.find("save_lambdas")->second);
		if (save_lambdas > 0)
		{
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				fn = rp->path + "delta_f_by_delta_s_" + to_string(tr_id) + cp->fn_suffix;
				save_double_vector(fn, ed->lambdas[tr_id], 16, false);

				fn = rp->path + "deltas_s_" + to_string(tr_id) + cp->fn_suffix;
				save_double_vector(fn, ed->deltas_s[tr_id], 16, false);

				fn = rp->path + "deltas_f_" + to_string(tr_id) + cp->fn_suffix;
				save_double_vector(fn, ed->deltas_f[tr_id], 16, false);
			}
		}
	}
}

void FloqSpinPhCoreBehaviour::dump_std_evo(AllData* ad) const
{
	RunParam* rp = ad->rp;
	ConfigParam* cp = ad->cp;
	ExpData* ed = ad->ed;

	int num_trajectories = cp->num_trajectories;
	int dump_num_total = ed->dump_num_total;

	int dump_obs = int(cp->params.find("dump_obs")->second);
	int jump = int(cp->params.find("jump")->second);

	if (dump_obs == 1)
	{
		int* dump_periods = ed->dump_periods;

		double* norm_evo = ed->norm_evo;
		MKL_Complex16* spec_evo = ed->spec_evo;

		string fn;

		fn = rp->path + "periods" + cp->fn_suffix;
		save_int_data(fn, dump_periods, dump_num_total, false);

		fn = rp->path + "norm_evo" + cp->fn_suffix;
		save_2d_inv_double_data(fn, norm_evo, dump_num_total, num_trajectories, 16, false);

		fn = rp->path + "spec_evo" + cp->fn_suffix;
		save_2d_inv_complex_data(fn, spec_evo, dump_num_total, num_trajectories, 16, false);

		int num_random_obs = int(cp->params.find("num_random_obs")->second);
		int random_obs_seed = int(cp->params.find("random_obs_seed")->second);
		for (int obs_id = 0; obs_id < num_random_obs; obs_id++)
		{
			fn = rp->path + "random_obs_evo_" + to_string(random_obs_seed + obs_id) + cp->fn_suffix;
			std::vector<std::vector<std::complex<double>>> random_obs_evo_to_dump;
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				random_obs_evo_to_dump.push_back(ed->random_obs_evo[tr_id][obs_id]);
			}
			save_2d_vector(fn, random_obs_evo_to_dump, 16);
		}

		if (jump > 0)
		{
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				fn = rp->path + "jump_times_" + to_string(tr_id) + cp->fn_suffix;
				save_double_vector(fn, ed->jump_times[tr_id], 16, false);

				fn = rp->path + "diss_types_" + to_string(tr_id) + cp->fn_suffix;
				save_int_vector(fn, ed->diss_types[tr_id], false);

				if (jump > 1)
				{
					fn = rp->path + "jump_norms_" + to_string(tr_id) + cp->fn_suffix;
					save_double_vector(fn, ed->jump_norms[tr_id], 16, false);

					fn = rp->path + "jump_etas_" + to_string(tr_id) + cp->fn_suffix;
					save_double_vector(fn, ed->jump_etas[tr_id], 16, false);
				}
			}
		}
	}
}

void FloqSpinPhCoreBehaviour::dump_lpn_evo(AllData* ad) const
{
	RunParam* rp = ad->rp;
	ConfigParam* cp = ad->cp;
	ExpData* ed = ad->ed;

	int num_trajectories = cp->num_trajectories;
	int dump_num_total = ed->dump_num_total;

	int dump_obs = int(cp->params.find("dump_obs")->second);

	if (dump_obs == 1)
	{
		double* lambda_evo = ed->lambda_evo;
		MKL_Complex16* spec_lpn_evo = ed->spec_lpn_evo;

		string fn;

		fn = rp->path + "lambda_evo" + cp->fn_suffix;
		save_2d_inv_double_data(fn, lambda_evo, dump_num_total, num_trajectories, 16, false);

		fn = rp->path + "spec_lpn_evo" + cp->fn_suffix;
		save_2d_inv_complex_data(fn, spec_lpn_evo, dump_num_total, num_trajectories, 16, false);

		int num_random_obs = int(cp->params.find("num_random_obs")->second);
		int random_obs_seed = int(cp->params.find("random_obs_seed")->second);
		for (int obs_id = 0; obs_id < num_random_obs; obs_id++)
		{
			fn = rp->path + "random_obs_lpn_evo_" + to_string(random_obs_seed + obs_id) + cp->fn_suffix;
			std::vector<std::vector<std::complex<double>>> random_obs_evo_to_dump;
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				random_obs_evo_to_dump.push_back(ed->random_obs_lpn_evo[tr_id][obs_id]);
			}
			save_2d_vector(fn, random_obs_evo_to_dump, 16);
		}
	}
}







Split * init_split_structure_dimer(AllData * ad)
{
	MainData * md = ad->md;

	Split * head = new Split[1];

	double T = md->T;
	int N = md->sys_size;
	int num_branches = md->num_ham_qj;

	head->prev = 0;
	head->type = false;
	head->dt = T;
	head->counter = num_branches;
	head->N = N;
	head->next = new Split[num_branches];

	for (int i = 0; i < num_branches; i++)
	{
		(head->next)[i].prev = head;
		init_split_branches(&((head->next)[i]), i, ad);
	}

	head->steps = md->num_diss;

	head->matrix = new MKL_Complex16[head->steps * N * N];
	head->g = new double[head->steps];

	for (int diss_id = 0; diss_id < head->steps; diss_id++)
	{
		head->g[diss_id] = 1.0;

		for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
		{
			for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
			{
				int index_xtd = diss_id * (md->sys_size * md->sys_size) + st_id_1 * md->sys_size + st_id_2;
				int index = st_id_1 * md->sys_size + st_id_2;

				head->matrix[index_xtd].real = md->dissipators[diss_id][index].real;
				head->matrix[index_xtd].imag = md->dissipators[diss_id][index].imag;
			}
		}
	}

	return head;
}

Split * init_split_structure_dimer_deep(AllData * ad)
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;

	int deep_num_steps = int(cp->params.find("deep_num_steps")->second);

	double T = 0.5 * md->T / double(deep_num_steps);
	int N = md->sys_size;
	int num_branches = md->num_ham_qj;

	Split * head = new Split[num_branches];

	for (int br_id = 0; br_id < num_branches; br_id++)
	{
		Split * branch = &head[br_id];
		branch->prev = 0;
		branch->type = false;
		branch->dt = T;
		branch->counter = 2;
		branch->N = N;
		branch->next = new Split[2];

		for (unsigned int i = 0; i < 2; i++)
		{
			(branch->next)[i].prev = branch;
			init_split_branches(&((branch->next)[i]), br_id, ad);
		}

		branch->steps = md->num_diss;

		branch->matrix = new MKL_Complex16[branch->steps * N * N];
		branch->g = new double[branch->steps];

		for (int diss_id = 0; diss_id < branch->steps; diss_id++)
		{
			branch->g[diss_id] = 1.0;

			for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
			{
				for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
				{
					int index_xtd = diss_id * (md->sys_size * md->sys_size) + st_id_1 * md->sys_size + st_id_2;
					int index = st_id_1 * md->sys_size + st_id_2;

					branch->matrix[index_xtd].real = md->dissipators[diss_id][index].real;
					branch->matrix[index_xtd].imag = md->dissipators[diss_id][index].imag;
				}
			}
		}
	}

	return head;
}

Split* init_split_structure_dimersync(AllData* ad)
{
	MainData* md = ad->md;

	Split* head = new Split[1];

	double T = md->T;
	int N = md->sys_size;
	int num_branches = md->num_ham_qj;

	head->prev = 0;
	head->type = false;
	head->dt = T;
	head->counter = num_branches;
	head->N = N;
	head->next = new Split[num_branches];

	for (int i = 0; i < num_branches; i++)
	{
		(head->next)[i].prev = head;
		init_split_branches(&((head->next)[i]), i, ad);
	}

	head->steps = md->num_diss;

	head->matrix = new MKL_Complex16[head->steps * N * N];
	head->g = new double[head->steps];

	for (int diss_id = 0; diss_id < head->steps; diss_id++)
	{
		head->g[diss_id] = 1.0;

		for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
		{
			for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
			{
				int index_xtd = diss_id * (md->sys_size * md->sys_size) + st_id_1 * md->sys_size + st_id_2;
				int index = st_id_1 * md->sys_size + st_id_2;

				head->matrix[index_xtd].real = md->dissipators[diss_id][index].real;
				head->matrix[index_xtd].imag = md->dissipators[diss_id][index].imag;
			}
		}
	}

	return head;
}

Split* init_split_structure_dimersync_deep(AllData* ad)
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;

	int deep_num_steps = int(cp->params.find("deep_num_steps")->second);

	double T = 0.5 * md->T / double(deep_num_steps);
	int N = md->sys_size;
	int num_branches = md->num_ham_qj;

	Split* head = new Split[num_branches];

	for (int br_id = 0; br_id < num_branches; br_id++)
	{
		Split* branch = &head[br_id];
		branch->prev = 0;
		branch->type = false;
		branch->dt = T;
		branch->counter = 2;
		branch->N = N;
		branch->next = new Split[2];

		for (unsigned int i = 0; i < 2; i++)
		{
			(branch->next)[i].prev = branch;
			init_split_branches(&((branch->next)[i]), br_id, ad);
		}

		branch->steps = md->num_diss;

		branch->matrix = new MKL_Complex16[branch->steps * N * N];
		branch->g = new double[branch->steps];

		for (int diss_id = 0; diss_id < branch->steps; diss_id++)
		{
			branch->g[diss_id] = 1.0;

			for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
			{
				for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
				{
					int index_xtd = diss_id * (md->sys_size * md->sys_size) + st_id_1 * md->sys_size + st_id_2;
					int index = st_id_1 * md->sys_size + st_id_2;

					branch->matrix[index_xtd].real = md->dissipators[diss_id][index].real;
					branch->matrix[index_xtd].imag = md->dissipators[diss_id][index].imag;
				}
			}
		}
	}

	return head;
}

Split * init_split_structure_jcs(AllData * ad)
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;

	double jcs_drv_part_1 = double(cp->params.find("jcs_drv_part_1")->second);
	double jcs_drv_part_2 = double(cp->params.find("jcs_drv_part_2")->second);

	double T_1 = jcs_drv_part_1 * md->T;
	double T_2 = jcs_drv_part_2 * md->T;

	int N = md->sys_size;
	int num_branches = md->num_ham_qj;

	Split * head = new Split[num_branches];

	for (int br_id = 0; br_id < num_branches; br_id++)
	{
		Split * branch = &head[br_id];
		branch->prev = 0;
		branch->type = false;

		if (br_id == 0)
		{
			branch->dt = T_1;
		}
		else if (br_id == 1)
		{
			branch->dt = T_2;
		}
		else
		{
			stringstream msg;
			msg << "Error: Wrong br_id" << endl;
			Error(msg.str());
		}

		branch->counter = 2;
		branch->N = N;
		branch->next = new Split[2];

		for (unsigned int i = 0; i < 2; i++)
		{
			(branch->next)[i].prev = branch;
			init_split_branches(&((branch->next)[i]), br_id, ad);
		}

		branch->steps = md->num_diss;

		branch->matrix = new MKL_Complex16[branch->steps * N * N];
		branch->g = new double[branch->steps];

		for (int diss_id = 0; diss_id < branch->steps; diss_id++)
		{
			branch->g[diss_id] = 1.0;

			for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
			{
				for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
				{
					int index_xtd = diss_id * (md->sys_size * md->sys_size) + st_id_1 * md->sys_size + st_id_2;
					int index = st_id_1 * md->sys_size + st_id_2;

					branch->matrix[index_xtd].real = md->dissipators[diss_id][index].real;
					branch->matrix[index_xtd].imag = md->dissipators[diss_id][index].imag;
				}
			}
		}
	}

	return head;
}

Split * init_split_structure_jcs_deep(AllData * ad)
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;

	double jcs_drv_part_1 = double(cp->params.find("jcs_drv_part_1")->second);
	double jcs_drv_part_2 = double(cp->params.find("jcs_drv_part_2")->second);

	int deep_num_steps = int(cp->params.find("deep_num_steps")->second);

	double T_1 = jcs_drv_part_1 * md->T / double(deep_num_steps);
	double T_2 = jcs_drv_part_2 * md->T / double(deep_num_steps);

	int N = md->sys_size;
	int num_branches = md->num_ham_qj;

	Split * head = new Split[num_branches];

	for (int br_id = 0; br_id < num_branches; br_id++)
	{
		Split * branch = &head[br_id];
		branch->prev = 0;
		branch->type = false;

		if (br_id == 0)
		{
			branch->dt = T_1;
		}
		else if (br_id == 1)
		{
			branch->dt = T_2;
		}
		else
		{
			stringstream msg;
			msg << "Error: Wrong br_id" << endl;
			Error(msg.str());
		}

		branch->counter = 2;
		branch->N = N;
		branch->next = new Split[2];

		for (unsigned int i = 0; i < 2; i++)
		{
			(branch->next)[i].prev = branch;
			init_split_branches(&((branch->next)[i]), br_id, ad);
		}

		branch->steps = md->num_diss;

		branch->matrix = new MKL_Complex16[branch->steps * N * N];
		branch->g = new double[branch->steps];

		for (int diss_id = 0; diss_id < branch->steps; diss_id++)
		{
			branch->g[diss_id] = 1.0;

			for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
			{
				for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
				{
					int index_xtd = diss_id * (md->sys_size * md->sys_size) + st_id_1 * md->sys_size + st_id_2;
					int index = st_id_1 * md->sys_size + st_id_2;

					branch->matrix[index_xtd].real = md->dissipators[diss_id][index].real;
					branch->matrix[index_xtd].imag = md->dissipators[diss_id][index].imag;
				}
			}
		}
	}

	return head;
}

Split * init_split_structure_ps(AllData * ad)
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;

	double ps_drv_part_1 = double(cp->params.find("ps_drv_part_1")->second);
	double ps_drv_part_2 = double(cp->params.find("ps_drv_part_2")->second);

	double T_1 = ps_drv_part_1 * md->T;
	double T_2 = ps_drv_part_2 * md->T;

	int sys_size = md->sys_size;
	int num_branches = md->num_ham_qj;

	Split * head = new Split[num_branches];

	for (int br_id = 0; br_id < num_branches; br_id++)
	{
		Split * branch = &head[br_id];
		branch->prev = 0;
		branch->type = false;

		if (br_id == 0)
		{
			branch->dt = T_1;
		}
		else if (br_id == 1)
		{
			branch->dt = T_2;
		}
		else
		{
			stringstream msg;
			msg << "Error: Wrong br_id" << endl;
			Error(msg.str());
		}

		branch->counter = 2;
		branch->N = sys_size;
		branch->next = new Split[2];

		for (unsigned int i = 0; i < 2; i++)
		{
			(branch->next)[i].prev = branch;
			init_split_branches(&((branch->next)[i]), br_id, ad);
		}

		branch->steps = md->num_diss;

		branch->matrix = new MKL_Complex16[branch->steps * sys_size * sys_size];
		branch->g = new double[branch->steps];

		for (int diss_id = 0; diss_id < branch->steps; diss_id++)
		{
			branch->g[diss_id] = 1.0;

			for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
			{
				for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
				{
					int index_xtd = diss_id * (md->sys_size * md->sys_size) + st_id_1 * md->sys_size + st_id_2;
					int index = st_id_1 * md->sys_size + st_id_2;

					branch->matrix[index_xtd].real = md->dissipators[diss_id][index].real;
					branch->matrix[index_xtd].imag = md->dissipators[diss_id][index].imag;
				}
			}
		}
	}

	return head;
}

Split * init_split_structure_ps_deep(AllData * ad)
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;

	double ps_drv_part_1 = double(cp->params.find("ps_drv_part_1")->second);
	double ps_drv_part_2 = double(cp->params.find("ps_drv_part_2")->second);

	int deep_num_steps = int(cp->params.find("deep_num_steps")->second);

	double T_1 = ps_drv_part_1 * md->T / double(deep_num_steps);
	double T_2 = ps_drv_part_2 * md->T / double(deep_num_steps);

	int N = md->sys_size;
	int num_branches = md->num_ham_qj;

	Split * head = new Split[num_branches];

	for (int br_id = 0; br_id < num_branches; br_id++)
	{
		Split * branch = &head[br_id];
		branch->prev = 0;
		branch->type = false;

		if (br_id == 0)
		{
			branch->dt = T_1;
		}
		else if (br_id == 1)
		{
			branch->dt = T_2;
		}
		else
		{
			stringstream msg;
			msg << "Error: Wrong br_id" << endl;
			Error(msg.str());
		}

		branch->counter = 2;
		branch->N = N;
		branch->next = new Split[2];

		for (unsigned int i = 0; i < 2; i++)
		{
			(branch->next)[i].prev = branch;
			init_split_branches(&((branch->next)[i]), br_id, ad);
		}

		branch->steps = md->num_diss;

		branch->matrix = new MKL_Complex16[branch->steps * N * N];
		branch->g = new double[branch->steps];

		for (int diss_id = 0; diss_id < branch->steps; diss_id++)
		{
			branch->g[diss_id] = 1.0;

			for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
			{
				for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
				{
					int index_xtd = diss_id * (md->sys_size * md->sys_size) + st_id_1 * md->sys_size + st_id_2;
					int index = st_id_1 * md->sys_size + st_id_2;

					branch->matrix[index_xtd].real = md->dissipators[diss_id][index].real;
					branch->matrix[index_xtd].imag = md->dissipators[diss_id][index].imag;
				}
			}
		}
	}

	return head;
}

Split * init_split_structure_mbl(AllData * ad)
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;

	double T = md->T;

	int sys_size = md->sys_size;
	int num_branches = md->num_ham_qj;

	Split * head = new Split[num_branches];

	for (int br_id = 0; br_id < num_branches; br_id++)
	{
		Split * branch = &head[br_id];
		branch->prev = 0;
		branch->type = false;
		branch->dt = T;
		branch->counter = 2;
		branch->N = sys_size;
		branch->next = new Split[2];

		for (unsigned int i = 0; i < 2; i++)
		{
			(branch->next)[i].prev = branch;
			init_split_branches(&((branch->next)[i]), br_id, ad);
		}

		branch->steps = md->num_diss;

		branch->matrix = new MKL_Complex16[branch->steps * sys_size * sys_size];
		branch->g = new double[branch->steps];

		for (int diss_id = 0; diss_id < branch->steps; diss_id++)
		{
			branch->g[diss_id] = 1.0;

			for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
			{
				for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
				{
					int index_xtd = diss_id * (md->sys_size * md->sys_size) + st_id_1 * md->sys_size + st_id_2;
					int index = st_id_1 * md->sys_size + st_id_2;

					branch->matrix[index_xtd].real = md->dissipators[diss_id][index].real;
					branch->matrix[index_xtd].imag = md->dissipators[diss_id][index].imag;
				}
			}
		}
	}

	return head;
}

Split * init_split_structure_mbl_deep(AllData * ad)
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;

	int deep_num_steps = int(cp->params.find("deep_num_steps")->second);

	double T = md->T / double(deep_num_steps);

	int N = md->sys_size;
	int num_branches = md->num_ham_qj;

	Split * head = new Split[num_branches];

	for (int br_id = 0; br_id < num_branches; br_id++)
	{
		Split * branch = &head[br_id];
		branch->prev = 0;
		branch->type = false;

		branch->dt = T;

		branch->counter = 2;
		branch->N = N;
		branch->next = new Split[2];

		for (unsigned int i = 0; i < 2; i++)
		{
			(branch->next)[i].prev = branch;
			init_split_branches(&((branch->next)[i]), br_id, ad);
		}

		branch->steps = md->num_diss;

		branch->matrix = new MKL_Complex16[branch->steps * N * N];
		branch->g = new double[branch->steps];

		for (int diss_id = 0; diss_id < branch->steps; diss_id++)
		{
			branch->g[diss_id] = 1.0;

			for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
			{
				for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
				{
					int index_xtd = diss_id * (md->sys_size * md->sys_size) + st_id_1 * md->sys_size + st_id_2;
					int index = st_id_1 * md->sys_size + st_id_2;

					branch->matrix[index_xtd].real = md->dissipators[diss_id][index].real;
					branch->matrix[index_xtd].imag = md->dissipators[diss_id][index].imag;
				}
			}
		}
	}

	return head;
}

Split* init_split_structure_lndham(AllData* ad)
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;

	double T = md->T;

	int sys_size = md->sys_size;
	int num_branches = md->num_ham_qj;

	Split* head = new Split[num_branches];

	for (int br_id = 0; br_id < num_branches; br_id++)
	{
		Split* branch = &head[br_id];
		branch->prev = 0;
		branch->type = false;
		branch->dt = T;
		branch->counter = 2;
		branch->N = sys_size;
		branch->next = new Split[2];

		for (unsigned int i = 0; i < 2; i++)
		{
			(branch->next)[i].prev = branch;
			init_split_branches(&((branch->next)[i]), br_id, ad);
		}

		branch->steps = md->num_diss;

		branch->matrix = new MKL_Complex16[branch->steps * sys_size * sys_size];
		branch->g = new double[branch->steps];

		for (int diss_id = 0; diss_id < branch->steps; diss_id++)
		{
			branch->g[diss_id] = 1.0;

			ds_mtx diss = ds_mtx(md->dissipators_eigen[diss_id]);

			for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
			{
				for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
				{
					int index_xtd = diss_id * (md->sys_size * md->sys_size) + st_id_1 * md->sys_size + st_id_2;
					int index = st_id_1 * md->sys_size + st_id_2;

					branch->matrix[index_xtd].real = diss(st_id_1, st_id_2).real();
					branch->matrix[index_xtd].imag = diss(st_id_1, st_id_2).imag();
				}
			}
		}
	}

	return head;
}

Split* init_split_structure_lndham_deep(AllData* ad)
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;

	int deep_num_steps = int(cp->params.find("deep_num_steps")->second);

	double T = md->T / double(deep_num_steps);

	int N = md->sys_size;
	int num_branches = md->num_ham_qj;

	Split* head = new Split[num_branches];

	for (int br_id = 0; br_id < num_branches; br_id++)
	{
		Split* branch = &head[br_id];
		branch->prev = 0;
		branch->type = false;

		branch->dt = T;

		branch->counter = 2;
		branch->N = N;
		branch->next = new Split[2];

		for (unsigned int i = 0; i < 2; i++)
		{
			(branch->next)[i].prev = branch;
			init_split_branches(&((branch->next)[i]), br_id, ad);
		}

		branch->steps = md->num_diss;

		branch->matrix = new MKL_Complex16[branch->steps * N * N];
		branch->g = new double[branch->steps];

		for (int diss_id = 0; diss_id < branch->steps; diss_id++)
		{
			branch->g[diss_id] = 1.0;

			ds_mtx diss = ds_mtx(md->dissipators_eigen[diss_id]);

			for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
			{
				for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
				{
					int index_xtd = diss_id * (md->sys_size * md->sys_size) + st_id_1 * md->sys_size + st_id_2;
					int index = st_id_1 * md->sys_size + st_id_2;

					branch->matrix[index_xtd].real = diss(st_id_1, st_id_2).real();
					branch->matrix[index_xtd].imag = diss(st_id_1, st_id_2).imag();
				}
			}
		}
	}

	return head;
}

Split* init_split_structure_integrable(AllData* ad)
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;

	double T = md->T;

	int sys_size = md->sys_size;
	int num_branches = md->num_ham_qj;

	Split* head = new Split[num_branches];

	for (int br_id = 0; br_id < num_branches; br_id++)
	{
		Split* branch = &head[br_id];
		branch->prev = 0;
		branch->type = false;
		branch->dt = T;
		branch->counter = 2;
		branch->N = sys_size;
		branch->next = new Split[2];

		for (unsigned int i = 0; i < 2; i++)
		{
			(branch->next)[i].prev = branch;
			init_split_branches(&((branch->next)[i]), br_id, ad);
		}

		branch->steps = md->num_diss;

		branch->matrix = new MKL_Complex16[branch->steps * sys_size * sys_size];
		branch->g = new double[branch->steps];

		for (int diss_id = 0; diss_id < branch->steps; diss_id++)
		{
			branch->g[diss_id] = 1.0;

			for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
			{
				for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
				{
					int index_xtd = diss_id * (md->sys_size * md->sys_size) + st_id_1 * md->sys_size + st_id_2;
					int index = st_id_1 * md->sys_size + st_id_2;

					branch->matrix[index_xtd].real = md->dissipators[diss_id][index].real;
					branch->matrix[index_xtd].imag = md->dissipators[diss_id][index].imag;
				}
			}
		}
	}

	return head;
}

Split* init_split_structure_integrable_deep(AllData* ad)
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;

	int deep_num_steps = int(cp->params.find("deep_num_steps")->second);

	double T = md->T / double(deep_num_steps);

	int N = md->sys_size;
	int num_branches = md->num_ham_qj;

	Split* head = new Split[num_branches];

	for (int br_id = 0; br_id < num_branches; br_id++)
	{
		Split* branch = &head[br_id];
		branch->prev = 0;
		branch->type = false;

		branch->dt = T;

		branch->counter = 2;
		branch->N = N;
		branch->next = new Split[2];

		for (unsigned int i = 0; i < 2; i++)
		{
			(branch->next)[i].prev = branch;
			init_split_branches(&((branch->next)[i]), br_id, ad);
		}

		branch->steps = md->num_diss;

		branch->matrix = new MKL_Complex16[branch->steps * N * N];
		branch->g = new double[branch->steps];

		for (int diss_id = 0; diss_id < branch->steps; diss_id++)
		{
			branch->g[diss_id] = 1.0;

			for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
			{
				for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
				{
					int index_xtd = diss_id * (md->sys_size * md->sys_size) + st_id_1 * md->sys_size + st_id_2;
					int index = st_id_1 * md->sys_size + st_id_2;

					branch->matrix[index_xtd].real = md->dissipators[diss_id][index].real;
					branch->matrix[index_xtd].imag = md->dissipators[diss_id][index].imag;
				}
			}
		}
	}

	return head;
}

Split* init_split_structure_floq2spins(AllData* ad)
{
	MainData* md = ad->md;

	Split* head = new Split[1];

	double T = md->T;
	int sys_size = md->sys_size;
	int num_branches = md->num_ham_qj;

	head->prev = 0;
	head->type = false;
	head->dt = T;
	head->counter = num_branches;
	head->N = sys_size;
	head->next = new Split[num_branches];

	for (int i = 0; i < num_branches; i++)
	{
		(head->next)[i].prev = head;
		init_split_branches(&((head->next)[i]), i, ad);
	}

	head->steps = md->num_diss;

	head->matrix = new MKL_Complex16[head->steps * sys_size * sys_size];
	head->g = new double[head->steps];

	for (int diss_id = 0; diss_id < head->steps; diss_id++)
	{
		head->g[diss_id] = 1.0;

		for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
		{
			for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
			{
				int index_xtd = diss_id * (md->sys_size * md->sys_size) + st_id_1 * md->sys_size + st_id_2;
				int index = st_id_1 * md->sys_size + st_id_2;

				head->matrix[index_xtd].real = md->dissipators[diss_id][index].real;
				head->matrix[index_xtd].imag = md->dissipators[diss_id][index].imag;
			}
		}
	}

	return head;
}

Split* init_split_structure_floq2spins_deep(AllData* ad)
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;

	int deep_num_steps = int(cp->params.find("deep_num_steps")->second);

	double T = 0.5 * md->T / double(deep_num_steps);
	int sys_size = md->sys_size;
	int num_branches = md->num_ham_qj;

	Split* head = new Split[num_branches];

	for (int br_id = 0; br_id < num_branches; br_id++)
	{
		Split* branch = &head[br_id];
		branch->prev = 0;
		branch->type = false;
		branch->dt = T;
		branch->counter = 2;
		branch->N = sys_size;
		branch->next = new Split[2];

		for (unsigned int i = 0; i < 2; i++)
		{
			(branch->next)[i].prev = branch;
			init_split_branches(&((branch->next)[i]), br_id, ad);
		}

		branch->steps = md->num_diss;

		branch->matrix = new MKL_Complex16[branch->steps * sys_size * sys_size];
		branch->g = new double[branch->steps];

		for (int diss_id = 0; diss_id < branch->steps; diss_id++)
		{
			branch->g[diss_id] = 1.0;

			for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
			{
				for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
				{
					int index_xtd = diss_id * (md->sys_size * md->sys_size) + st_id_1 * md->sys_size + st_id_2;
					int index = st_id_1 * md->sys_size + st_id_2;

					branch->matrix[index_xtd].real = md->dissipators[diss_id][index].real;
					branch->matrix[index_xtd].imag = md->dissipators[diss_id][index].imag;
				}
			}
		}
	}

	return head;
}

Split* init_split_structure_floqnspins(AllData* ad)
{
	MainData* md = ad->md;

	Split* head = new Split[1];

	double T = md->T;
	int sys_size = md->sys_size;
	int num_branches = md->num_ham_qj;

	head->prev = 0;
	head->type = false;
	head->dt = T;
	head->counter = num_branches;
	head->N = sys_size;
	head->next = new Split[num_branches];

	for (int i = 0; i < num_branches; i++)
	{
		(head->next)[i].prev = head;
		init_split_branches(&((head->next)[i]), i, ad);
	}

	head->steps = md->num_diss;

	head->matrix = new MKL_Complex16[head->steps * sys_size * sys_size];
	head->g = new double[head->steps];

	for (int diss_id = 0; diss_id < head->steps; diss_id++)
	{
		head->g[diss_id] = 1.0;

		for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
		{
			for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
			{
				int index_xtd = diss_id * (md->sys_size * md->sys_size) + st_id_1 * md->sys_size + st_id_2;
				int index = st_id_1 * md->sys_size + st_id_2;

				head->matrix[index_xtd].real = md->dissipators[diss_id][index].real;
				head->matrix[index_xtd].imag = md->dissipators[diss_id][index].imag;
			}
		}
	}

	return head;
}

Split* init_split_structure_floqnspins_deep(AllData* ad)
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;

	int deep_num_steps = int(cp->params.find("deep_num_steps")->second);

	double T = 0.5 * md->T / double(deep_num_steps);
	int sys_size = md->sys_size;
	int num_branches = md->num_ham_qj;

	Split* head = new Split[num_branches];

	for (int br_id = 0; br_id < num_branches; br_id++)
	{
		Split* branch = &head[br_id];
		branch->prev = 0;
		branch->type = false;
		branch->dt = T;
		branch->counter = 2;
		branch->N = sys_size;
		branch->next = new Split[2];

		for (unsigned int i = 0; i < 2; i++)
		{
			(branch->next)[i].prev = branch;
			init_split_branches(&((branch->next)[i]), br_id, ad);
		}

		branch->steps = md->num_diss;

		branch->matrix = new MKL_Complex16[branch->steps * sys_size * sys_size];
		branch->g = new double[branch->steps];

		for (int diss_id = 0; diss_id < branch->steps; diss_id++)
		{
			branch->g[diss_id] = 1.0;

			for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
			{
				for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
				{
					int index_xtd = diss_id * (md->sys_size * md->sys_size) + st_id_1 * md->sys_size + st_id_2;
					int index = st_id_1 * md->sys_size + st_id_2;

					branch->matrix[index_xtd].real = md->dissipators[diss_id][index].real;
					branch->matrix[index_xtd].imag = md->dissipators[diss_id][index].imag;
				}
			}
		}
	}

	return head;
}

Split* init_split_structure_flq_sp_ph(AllData* ad)
{
	MainData* md = ad->md;

	Split* head = new Split[1];

	double T = md->T;
	int sys_size = md->sys_size;
	int num_branches = md->num_ham_qj;

	head->prev = 0;
	head->type = false;
	head->dt = T;
	head->counter = num_branches;
	head->N = sys_size;
	head->next = new Split[num_branches];

	for (int i = 0; i < num_branches; i++)
	{
		(head->next)[i].prev = head;
		init_split_branches(&((head->next)[i]), i, ad);
	}

	head->steps = md->num_diss;

	head->matrix = new MKL_Complex16[head->steps * sys_size * sys_size];
	head->g = new double[head->steps];

	for (int diss_id = 0; diss_id < head->steps; diss_id++)
	{
		head->g[diss_id] = 1.0;

		for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
		{
			for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
			{
				int index_xtd = diss_id * (md->sys_size * md->sys_size) + st_id_1 * md->sys_size + st_id_2;
				int index = st_id_1 * md->sys_size + st_id_2;

				head->matrix[index_xtd].real = md->dissipators[diss_id][index].real;
				head->matrix[index_xtd].imag = md->dissipators[diss_id][index].imag;
			}
		}
	}

	return head;
}

Split* init_split_structure_flq_sp_ph_deep(AllData* ad)
{
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;

	int deep_num_steps = int(cp->params.find("deep_num_steps")->second);

	double T = 0.5 * md->T / double(deep_num_steps);
	int sys_size = md->sys_size;
	int num_branches = md->num_ham_qj;

	Split* head = new Split[num_branches];

	for (int br_id = 0; br_id < num_branches; br_id++)
	{
		Split* branch = &head[br_id];
		branch->prev = 0;
		branch->type = false;
		branch->dt = T;
		branch->counter = 2;
		branch->N = sys_size;
		branch->next = new Split[2];

		for (unsigned int i = 0; i < 2; i++)
		{
			(branch->next)[i].prev = branch;
			init_split_branches(&((branch->next)[i]), br_id, ad);
		}

		branch->steps = md->num_diss;

		branch->matrix = new MKL_Complex16[branch->steps * sys_size * sys_size];
		branch->g = new double[branch->steps];

		for (int diss_id = 0; diss_id < branch->steps; diss_id++)
		{
			branch->g[diss_id] = 1.0;

			for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
			{
				for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
				{
					int index_xtd = diss_id * (md->sys_size * md->sys_size) + st_id_1 * md->sys_size + st_id_2;
					int index = st_id_1 * md->sys_size + st_id_2;

					branch->matrix[index_xtd].real = md->dissipators[diss_id][index].real;
					branch->matrix[index_xtd].imag = md->dissipators[diss_id][index].imag;
				}
			}
		}
	}

	return head;
}



void rk_right_part_dimer(AllData * ad, int sub_step, int tr_id, int th_id)
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;

	double dimer_prm_E = double(cp->params.find("dimer_prm_E")->second);

	int dimer_drv_type = int(cp->params.find("dimer_drv_type")->second);
	double dimer_drv_ampl = double(cp->params.find("dimer_drv_ampl")->second);
	double dimer_drv_freq = double(cp->params.find("dimer_drv_freq")->second);
	double dimer_drv_phase = double(cp->params.find("dimer_drv_phase")->second);

	double time = ed->times_all[tr_id];

	MKL_Complex16 * arg = ed->args[th_id];
	MKL_Complex16 * non_drv_tmp = ed->non_drv_tmp[th_id];
	MKL_Complex16 * drv_tmp = ed->drv_tmp[th_id];

	MKL_Complex16 * k = ed->k1[th_id];
	if (sub_step == 1)
	{
		k = ed->k1[th_id];
	}
	else if (sub_step == 2)
	{
		k = ed->k2[th_id];
	}
	else if (sub_step == 3)
	{
		k = ed->k3[th_id];
	}
	else if (sub_step == 4)
	{
		k = ed->k4[th_id];
	}

	double E = 0.0;
	if (dimer_drv_type == 0)
	{
		double mod_time = fmod(time, md->T);
		double half_T = md->T * 0.5;
		if (mod_time < half_T)
		{
			E = dimer_prm_E + dimer_drv_ampl;
		}
		else
		{
			E = dimer_prm_E - dimer_drv_ampl;
		}
	}
	else if (dimer_drv_type == 1)
	{
		E = dimer_prm_E + dimer_drv_ampl * sin(dimer_drv_freq * time + dimer_drv_phase);
	}
	else
	{
		double mod_time = fmod(time, md->T);
		double half_T = md->T * 0.5;
		if (mod_time < half_T)
		{
			E = dimer_prm_E + dimer_drv_ampl;
		}
		else
		{
			E = dimer_prm_E - dimer_drv_ampl;
		}
	}

	MKL_Complex16 ZERO = { 0.0, 0.0 };
	MKL_Complex16 ONE = { 1.0, 0.0 };
	cblas_zgemv(CblasRowMajor, CblasNoTrans, sys_size, sys_size, &ONE, md->non_drv_part, sys_size, arg, 1, &ZERO, non_drv_tmp, 1);
	cblas_zgemv(CblasRowMajor, CblasNoTrans, sys_size, sys_size, &ONE, md->drv_part, sys_size, arg, 1, &ZERO, drv_tmp, 1);

	for (int st_id = 0; st_id < md->sys_size; st_id++)
	{
		k[st_id].real = +(ed->non_drv_tmp[th_id][st_id].imag + E * ed->drv_tmp[th_id][st_id].imag);
		k[st_id].imag = -(ed->non_drv_tmp[th_id][st_id].real + E * ed->drv_tmp[th_id][st_id].real);
	}
}

void rk_right_part_dimersync(AllData* ad, int sub_step, int tr_id, int th_id)
{
	RunParam* rp = ad->rp;
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int sys_size = md->sys_size;

	double prm_E = double(cp->params.find("dimersync_prm_E")->second);
	double prm_J = double(cp->params.find("dimersync_prm_J")->second);

	double drv_ampl_1 = double(cp->params.find("dimersync_drv_ampl_1")->second);
	double drv_freq_1 = double(cp->params.find("dimersync_drv_freq_1")->second);
	double drv_phase_1 = double(cp->params.find("dimersync_drv_phase_1")->second);
	double drv_ampl_2 = double(cp->params.find("dimersync_drv_ampl_2")->second);
	double drv_freq_2 = double(cp->params.find("dimersync_drv_freq_2")->second);
	double drv_phase_2 = double(cp->params.find("dimersync_drv_phase_2")->second);

	double time = ed->times_all[tr_id];

	MKL_Complex16* arg = ed->args[th_id];
	MKL_Complex16* non_drv_tmp = ed->non_drv_tmp[th_id];
	MKL_Complex16* drv_tmp = ed->drv_tmp[th_id];
	MKL_Complex16* drv_tmp_2 = ed->drv_tmp_2[th_id];

	MKL_Complex16* k = ed->k1[th_id];
	if (sub_step == 1)
	{
		k = ed->k1[th_id];
	}
	else if (sub_step == 2)
	{
		k = ed->k2[th_id];
	}
	else if (sub_step == 3)
	{
		k = ed->k3[th_id];
	}
	else if (sub_step == 4)
	{
		k = ed->k4[th_id];
	}

	double E = prm_E + drv_ampl_1 * sin(drv_freq_1 * time + drv_phase_1);
	double J = prm_J * drv_ampl_2 * sin(drv_freq_2 * time + drv_phase_2);

	MKL_Complex16 ZERO = { 0.0, 0.0 };
	MKL_Complex16 ONE = { 1.0, 0.0 };
	cblas_zgemv(CblasRowMajor, CblasNoTrans, sys_size, sys_size, &ONE, md->non_drv_part, sys_size, arg, 1, &ZERO, non_drv_tmp, 1);
	cblas_zgemv(CblasRowMajor, CblasNoTrans, sys_size, sys_size, &ONE, md->drv_part, sys_size, arg, 1, &ZERO, drv_tmp, 1);
	cblas_zgemv(CblasRowMajor, CblasNoTrans, sys_size, sys_size, &ONE, md->drv_part_2, sys_size, arg, 1, &ZERO, drv_tmp_2, 1);

	for (int st_id = 0; st_id < md->sys_size; st_id++)
	{
		k[st_id].real = +(ed->non_drv_tmp[th_id][st_id].imag + E * ed->drv_tmp[th_id][st_id].imag + J * ed->drv_tmp_2[th_id][st_id].imag);
		k[st_id].imag = -(ed->non_drv_tmp[th_id][st_id].real + E * ed->drv_tmp[th_id][st_id].real + J * ed->drv_tmp_2[th_id][st_id].real);
	}
}

void rk_right_part_jcs(AllData * ad, int sub_step, int tr_id, int th_id)
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;

	double jcs_drv_part_1 = double(cp->params.find("jcs_drv_part_1")->second);
	double jcs_drv_part_2 = double(cp->params.find("jcs_drv_part_2")->second);
	double jcs_drv_ampl = double(cp->params.find("jcs_drv_ampl")->second);

	double T_1 = jcs_drv_part_1 * md->T;
	double T_2 = jcs_drv_part_2 * md->T;
	double T = T_1 + T_2;

	double time = ed->times_all[tr_id];

	MKL_Complex16 * arg = ed->args[th_id];
	MKL_Complex16 * non_drv_tmp = ed->non_drv_tmp[th_id];
	MKL_Complex16 * drv_tmp = ed->drv_tmp[th_id];

	MKL_Complex16 * k = ed->k1[th_id];
	if (sub_step == 1)
	{
		k = ed->k1[th_id];
	}
	else if (sub_step == 2)
	{
		k = ed->k2[th_id];
	}
	else if (sub_step == 3)
	{
		k = ed->k3[th_id];
	}
	else if (sub_step == 4)
	{
		k = ed->k4[th_id];
	}

	double E = 0.0;
	double mod_time = fmod(time, T);
	if (mod_time < T_1)
	{
		E = jcs_drv_ampl;
	}
	else
	{
		E = 0.0;
	}

	MKL_Complex16 ZERO = { 0.0, 0.0 };
	MKL_Complex16 ONE = { 1.0, 0.0 };
	cblas_zgemv(CblasRowMajor, CblasNoTrans, sys_size, sys_size, &ONE, md->non_drv_part, sys_size, arg, 1, &ZERO, non_drv_tmp, 1);
	cblas_zgemv(CblasRowMajor, CblasNoTrans, sys_size, sys_size, &ONE, md->drv_part, sys_size, arg, 1, &ZERO, drv_tmp, 1);

	for (int st_id = 0; st_id < md->sys_size; st_id++)
	{
		k[st_id].real = +(ed->non_drv_tmp[th_id][st_id].imag + E * ed->drv_tmp[th_id][st_id].imag);
		k[st_id].imag = -(ed->non_drv_tmp[th_id][st_id].real + E * ed->drv_tmp[th_id][st_id].real);
	}
}

void rk_right_part_ps(AllData * ad, int sub_step, int tr_id, int th_id)
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;

	double ps_drv_part_1 = double(cp->params.find("ps_drv_part_1")->second);
	double ps_drv_part_2 = double(cp->params.find("ps_drv_part_2")->second);
	double ps_drv_ampl = double(cp->params.find("ps_drv_ampl")->second);

	double T_1 = ps_drv_part_1 * md->T;
	double T_2 = ps_drv_part_2 * md->T;
	double T = T_1 + T_2;

	double time = ed->times_all[tr_id];

	MKL_Complex16 * arg = ed->args[th_id];
	MKL_Complex16 * non_drv_tmp = ed->non_drv_tmp[th_id];
	MKL_Complex16 * drv_tmp = ed->drv_tmp[th_id];

	MKL_Complex16 * k = ed->k1[th_id];
	if (sub_step == 1)
	{
		k = ed->k1[th_id];
	}
	else if (sub_step == 2)
	{
		k = ed->k2[th_id];
	}
	else if (sub_step == 3)
	{
		k = ed->k3[th_id];
	}
	else if (sub_step == 4)
	{
		k = ed->k4[th_id];
	}

	double E = 0.0;
	double mod_time = fmod(time, T);
	if (mod_time < T_1)
	{
		E = ps_drv_ampl;
	}
	else
	{
		E = 0.0;
	}

	MKL_Complex16 ZERO = { 0.0, 0.0 };
	MKL_Complex16 ONE = { 1.0, 0.0 };
	cblas_zgemv(CblasRowMajor, CblasNoTrans, sys_size, sys_size, &ONE, md->non_drv_part, sys_size, arg, 1, &ZERO, non_drv_tmp, 1);
	cblas_zgemv(CblasRowMajor, CblasNoTrans, sys_size, sys_size, &ONE, md->drv_part, sys_size, arg, 1, &ZERO, drv_tmp, 1);

	for (int st_id = 0; st_id < md->sys_size; st_id++)
	{
		k[st_id].real = +(ed->non_drv_tmp[th_id][st_id].imag + E * ed->drv_tmp[th_id][st_id].imag);
		k[st_id].imag = -(ed->non_drv_tmp[th_id][st_id].real + E * ed->drv_tmp[th_id][st_id].real);
	}
}

void rk_right_part_mbl(AllData * ad, int sub_step, int tr_id, int th_id)
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;

	double T = md->T;

	double time = ed->times_all[tr_id];

	MKL_Complex16 * arg = ed->args[th_id];
	MKL_Complex16 * non_drv_tmp = ed->non_drv_tmp[th_id];

	MKL_Complex16 * k = ed->k1[th_id];
	if (sub_step == 1)
	{
		k = ed->k1[th_id];
	}
	else if (sub_step == 2)
	{
		k = ed->k2[th_id];
	}
	else if (sub_step == 3)
	{
		k = ed->k3[th_id];
	}
	else if (sub_step == 4)
	{
		k = ed->k4[th_id];
	}

	MKL_Complex16 ZERO = { 0.0, 0.0 };
	MKL_Complex16 ONE = { 1.0, 0.0 };
	cblas_zgemv(CblasRowMajor, CblasNoTrans, sys_size, sys_size, &ONE, md->non_drv_part, sys_size, arg, 1, &ZERO, non_drv_tmp, 1);

	for (int st_id = 0; st_id < md->sys_size; st_id++)
	{
		k[st_id].real = +(ed->non_drv_tmp[th_id][st_id].imag);
		k[st_id].imag = -(ed->non_drv_tmp[th_id][st_id].real);
	}
}

void rk_right_part_lndham(AllData* ad, int sub_step, int tr_id, int th_id)
{
	RunParam* rp = ad->rp;
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int sys_size = md->sys_size;

	double T = md->T;

	double time = ed->times_all[tr_id];

	MKL_Complex16* arg = ed->args[th_id];
	MKL_Complex16* non_drv_tmp = ed->non_drv_tmp[th_id];

	MKL_Complex16* k = ed->k1[th_id];
	if (sub_step == 1)
	{
		k = ed->k1[th_id];
	}
	else if (sub_step == 2)
	{
		k = ed->k2[th_id];
	}
	else if (sub_step == 3)
	{
		k = ed->k3[th_id];
	}
	else if (sub_step == 4)
	{
		k = ed->k4[th_id];
	}

	MKL_Complex16 ZERO = { 0.0, 0.0 };
	MKL_Complex16 ONE = { 1.0, 0.0 };
	cblas_zgemv(CblasRowMajor, CblasNoTrans, sys_size, sys_size, &ONE, md->non_drv_part, sys_size, arg, 1, &ZERO, non_drv_tmp, 1);

	for (int st_id = 0; st_id < md->sys_size; st_id++)
	{
		k[st_id].real = +(ed->non_drv_tmp[th_id][st_id].imag);
		k[st_id].imag = -(ed->non_drv_tmp[th_id][st_id].real);
	}
}

void rk_right_part_integrable(AllData* ad, int sub_step, int tr_id, int th_id)
{
	RunParam* rp = ad->rp;
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int sys_size = md->sys_size;

	double T = md->T;

	double time = ed->times_all[tr_id];

	MKL_Complex16* arg = ed->args[th_id];
	MKL_Complex16* non_drv_tmp = ed->non_drv_tmp[th_id];

	MKL_Complex16* k = ed->k1[th_id];
	if (sub_step == 1)
	{
		k = ed->k1[th_id];
	}
	else if (sub_step == 2)
	{
		k = ed->k2[th_id];
	}
	else if (sub_step == 3)
	{
		k = ed->k3[th_id];
	}
	else if (sub_step == 4)
	{
		k = ed->k4[th_id];
	}

	MKL_Complex16 ZERO = { 0.0, 0.0 };
	MKL_Complex16 ONE = { 1.0, 0.0 };
	cblas_zgemv(CblasRowMajor, CblasNoTrans, sys_size, sys_size, &ONE, md->non_drv_part, sys_size, arg, 1, &ZERO, non_drv_tmp, 1);

	for (int st_id = 0; st_id < md->sys_size; st_id++)
	{
		k[st_id].real = +(ed->non_drv_tmp[th_id][st_id].imag);
		k[st_id].imag = -(ed->non_drv_tmp[th_id][st_id].real);
	}
}

void rk_right_part_floq2spins(AllData* ad, int sub_step, int tr_id, int th_id)
{
	RunParam* rp = ad->rp;
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int sys_size = md->sys_size;

	double prm_ampl = double(cp->params.find("flq2sp_ampl")->second);
	double prm_freq = double(cp->params.find("flq2sp_freq")->second);

	double time = ed->times_all[tr_id];

	MKL_Complex16* arg = ed->args[th_id];
	MKL_Complex16* non_drv_tmp = ed->non_drv_tmp[th_id];
	MKL_Complex16* drv_tmp = ed->drv_tmp[th_id];

	MKL_Complex16* k = ed->k1[th_id];
	if (sub_step == 1)
	{
		k = ed->k1[th_id];
	}
	else if (sub_step == 2)
	{
		k = ed->k2[th_id];
	}
	else if (sub_step == 3)
	{
		k = ed->k3[th_id];
	}
	else if (sub_step == 4)
	{
		k = ed->k4[th_id];
	}

	double modulation = prm_ampl * sin(prm_freq * time);

	MKL_Complex16 ZERO = { 0.0, 0.0 };
	MKL_Complex16 ONE = { 1.0, 0.0 };
	cblas_zgemv(CblasRowMajor, CblasNoTrans, sys_size, sys_size, &ONE, md->non_drv_part, sys_size, arg, 1, &ZERO, non_drv_tmp, 1);
	cblas_zgemv(CblasRowMajor, CblasNoTrans, sys_size, sys_size, &ONE, md->drv_part, sys_size, arg, 1, &ZERO, drv_tmp, 1);

	for (int st_id = 0; st_id < md->sys_size; st_id++)
	{
		k[st_id].real = +(ed->non_drv_tmp[th_id][st_id].imag + modulation * ed->drv_tmp[th_id][st_id].imag);
		k[st_id].imag = -(ed->non_drv_tmp[th_id][st_id].real + modulation * ed->drv_tmp[th_id][st_id].real);
	}
}

void rk_right_part_floqnspins(AllData* ad, int sub_step, int tr_id, int th_id)
{
	RunParam* rp = ad->rp;
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int sys_size = md->sys_size;

	double prm_ampl = double(cp->params.find("flqnsp_ampl")->second);
	double prm_freq = double(cp->params.find("flqnsp_freq")->second);

	double time = ed->times_all[tr_id];

	MKL_Complex16* arg = ed->args[th_id];
	MKL_Complex16* non_drv_tmp = ed->non_drv_tmp[th_id];
	MKL_Complex16* drv_tmp = ed->drv_tmp[th_id];

	MKL_Complex16* k = ed->k1[th_id];
	if (sub_step == 1)
	{
		k = ed->k1[th_id];
	}
	else if (sub_step == 2)
	{
		k = ed->k2[th_id];
	}
	else if (sub_step == 3)
	{
		k = ed->k3[th_id];
	}
	else if (sub_step == 4)
	{
		k = ed->k4[th_id];
	}

	double modulation = prm_ampl * sin(prm_freq * time);

	MKL_Complex16 ZERO = { 0.0, 0.0 };
	MKL_Complex16 ONE = { 1.0, 0.0 };
	cblas_zgemv(CblasRowMajor, CblasNoTrans, sys_size, sys_size, &ONE, md->non_drv_part, sys_size, arg, 1, &ZERO, non_drv_tmp, 1);
	cblas_zgemv(CblasRowMajor, CblasNoTrans, sys_size, sys_size, &ONE, md->drv_part, sys_size, arg, 1, &ZERO, drv_tmp, 1);

	for (int st_id = 0; st_id < md->sys_size; st_id++)
	{
		k[st_id].real = +(ed->non_drv_tmp[th_id][st_id].imag + modulation * ed->drv_tmp[th_id][st_id].imag);
		k[st_id].imag = -(ed->non_drv_tmp[th_id][st_id].real + modulation * ed->drv_tmp[th_id][st_id].real);
	}
}

void rk_right_part_flq_sp_ph(AllData* ad, int sub_step, int tr_id, int th_id)
{
	RunParam* rp = ad->rp;
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int sys_size = md->sys_size;

	double prm_ampl = double(cp->params.find("flq_sp_ph_ampl")->second);
	double prm_freq = double(cp->params.find("flq_sp_ph_freq")->second);

	double time = ed->times_all[tr_id];

	MKL_Complex16* arg = ed->args[th_id];
	MKL_Complex16* non_drv_tmp = ed->non_drv_tmp[th_id];
	MKL_Complex16* drv_tmp = ed->drv_tmp[th_id];

	MKL_Complex16* k = ed->k1[th_id];
	if (sub_step == 1)
	{
		k = ed->k1[th_id];
	}
	else if (sub_step == 2)
	{
		k = ed->k2[th_id];
	}
	else if (sub_step == 3)
	{
		k = ed->k3[th_id];
	}
	else if (sub_step == 4)
	{
		k = ed->k4[th_id];
	}

	double modulation = prm_ampl * sin(prm_freq * time);

	MKL_Complex16 ZERO = { 0.0, 0.0 };
	MKL_Complex16 ONE = { 1.0, 0.0 };
	cblas_zgemv(CblasRowMajor, CblasNoTrans, sys_size, sys_size, &ONE, md->non_drv_part, sys_size, arg, 1, &ZERO, non_drv_tmp, 1);
	cblas_zgemv(CblasRowMajor, CblasNoTrans, sys_size, sys_size, &ONE, md->drv_part, sys_size, arg, 1, &ZERO, drv_tmp, 1);

	for (int st_id = 0; st_id < md->sys_size; st_id++)
	{
		k[st_id].real = +(ed->non_drv_tmp[th_id][st_id].imag + modulation * ed->drv_tmp[th_id][st_id].imag);
		k[st_id].imag = -(ed->non_drv_tmp[th_id][st_id].real + modulation * ed->drv_tmp[th_id][st_id].real);
	}
}



void rk_int_dimer(AllData * ad, int tr_id, int th_id, double step)
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;

	set_init_args(ad, tr_id, th_id);

	rk_right_part_dimer(ad, 1, tr_id, th_id);
	arg_upd(ad, 1, tr_id, th_id);

	ed->times_all[tr_id] += step * 0.5;

	rk_right_part_dimer(ad, 2, tr_id, th_id);
	arg_upd(ad, 2, tr_id, th_id);

	rk_right_part_dimer(ad, 3, tr_id, th_id);
	arg_upd(ad, 3, tr_id, th_id);

	ed->times_all[tr_id] += step * 0.5;

	rk_right_part_dimer(ad, 4, tr_id, th_id);

	rk_final(ad, tr_id, th_id, step);
}

void rk_int_dimersync(AllData* ad, int tr_id, int th_id, double step)
{
	RunParam* rp = ad->rp;
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int sys_size = md->sys_size;

	set_init_args(ad, tr_id, th_id);

	rk_right_part_dimersync(ad, 1, tr_id, th_id);
	arg_upd(ad, 1, tr_id, th_id);

	ed->times_all[tr_id] += step * 0.5;

	rk_right_part_dimersync(ad, 2, tr_id, th_id);
	arg_upd(ad, 2, tr_id, th_id);

	rk_right_part_dimersync(ad, 3, tr_id, th_id);
	arg_upd(ad, 3, tr_id, th_id);

	ed->times_all[tr_id] += step * 0.5;

	rk_right_part_dimersync(ad, 4, tr_id, th_id);

	rk_final(ad, tr_id, th_id, step);
}

void rk_int_jcs(AllData * ad, int tr_id, int th_id, double step)
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;

	set_init_args(ad, tr_id, th_id);

	rk_right_part_jcs(ad, 1, tr_id, th_id);
	arg_upd(ad, 1, tr_id, th_id);

	ed->times_all[tr_id] += step * 0.5;

	rk_right_part_jcs(ad, 2, tr_id, th_id);
	arg_upd(ad, 2, tr_id, th_id);

	rk_right_part_jcs(ad, 3, tr_id, th_id);
	arg_upd(ad, 3, tr_id, th_id);

	ed->times_all[tr_id] += step * 0.5;

	rk_right_part_jcs(ad, 4, tr_id, th_id);

	rk_final(ad, tr_id, th_id, step);
}

void rk_int_ps(AllData * ad, int tr_id, int th_id, double step)
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;

	set_init_args(ad, tr_id, th_id);

	rk_right_part_ps(ad, 1, tr_id, th_id);
	arg_upd(ad, 1, tr_id, th_id);

	ed->times_all[tr_id] += step * 0.5;

	rk_right_part_ps(ad, 2, tr_id, th_id);
	arg_upd(ad, 2, tr_id, th_id);

	rk_right_part_ps(ad, 3, tr_id, th_id);
	arg_upd(ad, 3, tr_id, th_id);

	ed->times_all[tr_id] += step * 0.5;

	rk_right_part_ps(ad, 4, tr_id, th_id);

	rk_final(ad, tr_id, th_id, step);
}

void rk_int_mbl(AllData * ad, int tr_id, int th_id, double step)
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;

	set_init_args(ad, tr_id, th_id);

	rk_right_part_mbl(ad, 1, tr_id, th_id);
	arg_upd(ad, 1, tr_id, th_id);

	ed->times_all[tr_id] += step * 0.5;

	rk_right_part_mbl(ad, 2, tr_id, th_id);
	arg_upd(ad, 2, tr_id, th_id);

	rk_right_part_mbl(ad, 3, tr_id, th_id);
	arg_upd(ad, 3, tr_id, th_id);

	ed->times_all[tr_id] += step * 0.5;

	rk_right_part_mbl(ad, 4, tr_id, th_id);

	rk_final(ad, tr_id, th_id, step);
}

void rk_int_lndham(AllData* ad, int tr_id, int th_id, double step)
{
	RunParam* rp = ad->rp;
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int sys_size = md->sys_size;

	set_init_args(ad, tr_id, th_id);

	rk_right_part_lndham(ad, 1, tr_id, th_id);
	arg_upd(ad, 1, tr_id, th_id);

	ed->times_all[tr_id] += step * 0.5;

	rk_right_part_lndham(ad, 2, tr_id, th_id);
	arg_upd(ad, 2, tr_id, th_id);

	rk_right_part_lndham(ad, 3, tr_id, th_id);
	arg_upd(ad, 3, tr_id, th_id);

	ed->times_all[tr_id] += step * 0.5;

	rk_right_part_lndham(ad, 4, tr_id, th_id);

	rk_final(ad, tr_id, th_id, step);
}

void rk_int_integrable(AllData* ad, int tr_id, int th_id, double step)
{
	RunParam* rp = ad->rp;
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int sys_size = md->sys_size;

	set_init_args(ad, tr_id, th_id);

	rk_right_part_integrable(ad, 1, tr_id, th_id);
	arg_upd(ad, 1, tr_id, th_id);

	ed->times_all[tr_id] += step * 0.5;

	rk_right_part_integrable(ad, 2, tr_id, th_id);
	arg_upd(ad, 2, tr_id, th_id);

	rk_right_part_integrable(ad, 3, tr_id, th_id);
	arg_upd(ad, 3, tr_id, th_id);

	ed->times_all[tr_id] += step * 0.5;

	rk_right_part_integrable(ad, 4, tr_id, th_id);

	rk_final(ad, tr_id, th_id, step);
}

void rk_int_floq2spins(AllData* ad, int tr_id, int th_id, double step)
{
	RunParam* rp = ad->rp;
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int sys_size = md->sys_size;

	set_init_args(ad, tr_id, th_id);

	rk_right_part_floq2spins(ad, 1, tr_id, th_id);
	arg_upd(ad, 1, tr_id, th_id);

	ed->times_all[tr_id] += step * 0.5;

	rk_right_part_floq2spins(ad, 2, tr_id, th_id);
	arg_upd(ad, 2, tr_id, th_id);

	rk_right_part_floq2spins(ad, 3, tr_id, th_id);
	arg_upd(ad, 3, tr_id, th_id);

	ed->times_all[tr_id] += step * 0.5;

	rk_right_part_floq2spins(ad, 4, tr_id, th_id);

	rk_final(ad, tr_id, th_id, step);
}

void rk_int_floqnspins(AllData* ad, int tr_id, int th_id, double step)
{
	RunParam* rp = ad->rp;
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int sys_size = md->sys_size;

	set_init_args(ad, tr_id, th_id);

	rk_right_part_floqnspins(ad, 1, tr_id, th_id);
	arg_upd(ad, 1, tr_id, th_id);

	ed->times_all[tr_id] += step * 0.5;

	rk_right_part_floqnspins(ad, 2, tr_id, th_id);
	arg_upd(ad, 2, tr_id, th_id);

	rk_right_part_floqnspins(ad, 3, tr_id, th_id);
	arg_upd(ad, 3, tr_id, th_id);

	ed->times_all[tr_id] += step * 0.5;

	rk_right_part_floqnspins(ad, 4, tr_id, th_id);

	rk_final(ad, tr_id, th_id, step);
}

void rk_int_flq_sp_ph(AllData* ad, int tr_id, int th_id, double step)
{
	RunParam* rp = ad->rp;
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int sys_size = md->sys_size;

	set_init_args(ad, tr_id, th_id);

	rk_right_part_flq_sp_ph(ad, 1, tr_id, th_id);
	arg_upd(ad, 1, tr_id, th_id);

	ed->times_all[tr_id] += step * 0.5;

	rk_right_part_flq_sp_ph(ad, 2, tr_id, th_id);
	arg_upd(ad, 2, tr_id, th_id);

	rk_right_part_flq_sp_ph(ad, 3, tr_id, th_id);
	arg_upd(ad, 3, tr_id, th_id);

	ed->times_all[tr_id] += step * 0.5;

	rk_right_part_flq_sp_ph(ad, 4, tr_id, th_id);

	rk_final(ad, tr_id, th_id, step);
}





void rk_step_dimer(AllData * ad, int tr_id, int th_id, double step)
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;

	int jump = int(cp->params.find("jump")->second);

	VSLStreamStatePtr * stream = &(ed->streams[tr_id]);
	MKL_Complex16 * phi = &(ed->phi_all[tr_id * sys_size]);
	double * eta = &(ed->etas_all[tr_id]);

	double prev_norm = 0.0;
	double curr_norm = 0.0;
	double norm_diff = 0.0;
	double begin_part = 0.0;
	double end_part = 0.0;
	double begin_step = 0.0;
	double end_step = 0.0;

	save_phi_prev(ad, tr_id, th_id);
	prev_norm = norm_square(phi, sys_size);
	rk_int_dimer(ad, tr_id, th_id, step);

	if (is_norm_crossed(phi, eta, sys_size))
	{
		curr_norm = norm_square(phi, sys_size);
		norm_diff = prev_norm - curr_norm;
		begin_part = (prev_norm - *(eta)) / norm_diff;
		end_part = 1.0 - begin_part;
		begin_step = step * begin_part;
		end_step = step * end_part;

		restore_from_prev(ad, tr_id, th_id, step);

		rk_int_dimer(ad, tr_id, th_id, begin_step);

		if (jump > 0 && ed->is_obs == 1)
		{
			double jump_time = ed->times_all[tr_id];
			double jump_norm = norm_square(phi, sys_size);
			double jump_eta = *eta;

			ed->jump_times[tr_id].push_back(jump_time);
			ed->jump_norms[tr_id].push_back(jump_norm);
			ed->jump_etas[tr_id].push_back(jump_eta);

			ed->jumps_counts[tr_id]++;
		}

		rk_recovery(ad, tr_id, th_id);
		vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, *stream, 1, eta, 0.0, 1.0);
		while (*eta == 0.0)
		{
			vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, *stream, 1, eta, 0.0, 1.0);
		}

		rk_step_dimer(ad, tr_id, th_id, end_step);
	}
}

void rk_step_dimersync(AllData* ad, int tr_id, int th_id, double step)
{
	RunParam* rp = ad->rp;
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int sys_size = md->sys_size;

	int jump = int(cp->params.find("jump")->second);

	VSLStreamStatePtr* stream = &(ed->streams[tr_id]);
	MKL_Complex16* phi = &(ed->phi_all[tr_id * sys_size]);
	double* eta = &(ed->etas_all[tr_id]);

	double prev_norm = 0.0;
	double curr_norm = 0.0;
	double norm_diff = 0.0;
	double begin_part = 0.0;
	double end_part = 0.0;
	double begin_step = 0.0;
	double end_step = 0.0;

	save_phi_prev(ad, tr_id, th_id);
	prev_norm = norm_square(phi, sys_size);
	rk_int_dimersync(ad, tr_id, th_id, step);

	if (is_norm_crossed(phi, eta, sys_size))
	{
		curr_norm = norm_square(phi, sys_size);
		norm_diff = prev_norm - curr_norm;
		begin_part = (prev_norm - *(eta)) / norm_diff;
		end_part = 1.0 - begin_part;
		begin_step = step * begin_part;
		end_step = step * end_part;

		restore_from_prev(ad, tr_id, th_id, step);

		rk_int_dimersync(ad, tr_id, th_id, begin_step);

		if (jump > 0 && ed->is_obs == 1)
		{
			double jump_time = ed->times_all[tr_id];
			double jump_norm = norm_square(phi, sys_size);
			double jump_eta = *eta;

			ed->jump_times[tr_id].push_back(jump_time);
			ed->jump_norms[tr_id].push_back(jump_norm);
			ed->jump_etas[tr_id].push_back(jump_eta);

			ed->jumps_counts[tr_id]++;
		}

		rk_recovery(ad, tr_id, th_id);
		vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, *stream, 1, eta, 0.0, 1.0);
		while (*eta == 0.0)
		{
			vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, *stream, 1, eta, 0.0, 1.0);
		}

		rk_step_dimersync(ad, tr_id, th_id, end_step);
	}
}

void rk_step_jcs(AllData * ad, int tr_id, int th_id, double step)
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;

	int jump = int(cp->params.find("jump")->second);

	VSLStreamStatePtr * stream = &(ed->streams[tr_id]);
	MKL_Complex16 * phi = &(ed->phi_all[tr_id * sys_size]);
	double * eta = &(ed->etas_all[tr_id]);

	double prev_norm = 0.0;
	double curr_norm = 0.0;
	double norm_diff = 0.0;
	double begin_part = 0.0;
	double end_part = 0.0;
	double begin_step = 0.0;
	double end_step = 0.0;

	save_phi_prev(ad, tr_id, th_id);
	prev_norm = norm_square(phi, sys_size);
	rk_int_jcs(ad, tr_id, th_id, step);

	if (is_norm_crossed(phi, eta, sys_size))
	{
		curr_norm = norm_square(phi, sys_size);
		norm_diff = prev_norm - curr_norm;
		begin_part = (prev_norm - *(eta)) / norm_diff;
		end_part = 1.0 - begin_part;
		begin_step = step * begin_part;
		end_step = step * end_part;

		restore_from_prev(ad, tr_id, th_id, step);

		rk_int_jcs(ad, tr_id, th_id, begin_step);

		if (jump > 0 && ed->is_obs == 1)
		{
			double jump_time = ed->times_all[tr_id];
			double jump_norm = norm_square(phi, sys_size);
			double jump_eta = *eta;

			ed->jump_times[tr_id].push_back(jump_time);
			ed->jump_norms[tr_id].push_back(jump_norm);
			ed->jump_etas[tr_id].push_back(jump_eta);
		
			ed->jumps_counts[tr_id]++;
		}

		rk_recovery(ad, tr_id, th_id);
		vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, *stream, 1, eta, 0.0, 1.0);
		while (*eta == 0.0)
		{
			vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, *stream, 1, eta, 0.0, 1.0);
		}

		rk_step_jcs(ad, tr_id, th_id, end_step);
	}
}

void rk_step_ps(AllData * ad, int tr_id, int th_id, double step)
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;

	int jump = int(cp->params.find("jump")->second);

	VSLStreamStatePtr * stream = &(ed->streams[tr_id]);
	MKL_Complex16 * phi = &(ed->phi_all[tr_id * sys_size]);
	double * eta = &(ed->etas_all[tr_id]);

	double prev_norm = 0.0;
	double curr_norm = 0.0;
	double norm_diff = 0.0;
	double begin_part = 0.0;
	double end_part = 0.0;
	double begin_step = 0.0;
	double end_step = 0.0;

	save_phi_prev(ad, tr_id, th_id);
	prev_norm = norm_square(phi, sys_size);
	rk_int_ps(ad, tr_id, th_id, step);

	if (is_norm_crossed(phi, eta, sys_size))
	{
		curr_norm = norm_square(phi, sys_size);
		norm_diff = prev_norm - curr_norm;
		begin_part = (prev_norm - *(eta)) / norm_diff;
		end_part = 1.0 - begin_part;
		begin_step = step * begin_part;
		end_step = step * end_part;

		restore_from_prev(ad, tr_id, th_id, step);

		rk_int_ps(ad, tr_id, th_id, begin_step);

		if (jump > 0 && ed->is_obs == 1)
		{
			double jump_time = ed->times_all[tr_id];
			double jump_norm = norm_square(phi, sys_size);
			double jump_eta = *eta;

			ed->jump_times[tr_id].push_back(jump_time);
			ed->jump_norms[tr_id].push_back(jump_norm);
			ed->jump_etas[tr_id].push_back(jump_eta);

			ed->jumps_counts[tr_id]++;
		}

		rk_recovery(ad, tr_id, th_id);
		vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, *stream, 1, eta, 0.0, 1.0);
		while (*eta == 0.0)
		{
			vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, *stream, 1, eta, 0.0, 1.0);
		}

		rk_step_ps(ad, tr_id, th_id, end_step);
	}
}

void rk_step_mbl(AllData * ad, int tr_id, int th_id, double step)
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;

	int jump = int(cp->params.find("jump")->second);

	VSLStreamStatePtr * stream = &(ed->streams[tr_id]);
	MKL_Complex16 * phi = &(ed->phi_all[tr_id * sys_size]);
	double * eta = &(ed->etas_all[tr_id]);

	double prev_norm = 0.0;
	double curr_norm = 0.0;
	double norm_diff = 0.0;
	double begin_part = 0.0;
	double end_part = 0.0;
	double begin_step = 0.0;
	double end_step = 0.0;

	save_phi_prev(ad, tr_id, th_id);
	prev_norm = norm_square(phi, sys_size);
	rk_int_mbl(ad, tr_id, th_id, step);

	if (is_norm_crossed(phi, eta, sys_size))
	{
		curr_norm = norm_square(phi, sys_size);
		norm_diff = prev_norm - curr_norm;
		begin_part = (prev_norm - *(eta)) / norm_diff;
		end_part = 1.0 - begin_part;
		begin_step = step * begin_part;
		end_step = step * end_part;

		restore_from_prev(ad, tr_id, th_id, step);

		rk_int_mbl(ad, tr_id, th_id, begin_step);

		if (jump > 0 && ed->is_obs == 1)
		{
			double jump_time = ed->times_all[tr_id];
			double jump_norm = norm_square(phi, sys_size);
			double jump_eta = *eta;

			ed->jump_times[tr_id].push_back(jump_time);
			ed->jump_norms[tr_id].push_back(jump_norm);
			ed->jump_etas[tr_id].push_back(jump_eta);

			ed->jumps_counts[tr_id]++;
		}

		rk_recovery(ad, tr_id, th_id);
		vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, *stream, 1, eta, 0.0, 1.0);
		while (*eta == 0.0)
		{
			vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, *stream, 1, eta, 0.0, 1.0);
		}

		rk_step_mbl(ad, tr_id, th_id, end_step);
	}
}

void rk_step_lndham(AllData* ad, int tr_id, int th_id, double step)
{
	RunParam* rp = ad->rp;
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int sys_size = md->sys_size;

	int jump = int(cp->params.find("jump")->second);

	VSLStreamStatePtr* stream = &(ed->streams[tr_id]);
	MKL_Complex16* phi = &(ed->phi_all[tr_id * sys_size]);
	double* eta = &(ed->etas_all[tr_id]);

	double prev_norm = 0.0;
	double curr_norm = 0.0;
	double norm_diff = 0.0;
	double begin_part = 0.0;
	double end_part = 0.0;
	double begin_step = 0.0;
	double end_step = 0.0;

	save_phi_prev(ad, tr_id, th_id);
	prev_norm = norm_square(phi, sys_size);
	rk_int_lndham(ad, tr_id, th_id, step);

	if (is_norm_crossed(phi, eta, sys_size))
	{
		curr_norm = norm_square(phi, sys_size);
		norm_diff = prev_norm - curr_norm;
		begin_part = (prev_norm - *(eta)) / norm_diff;
		end_part = 1.0 - begin_part;
		begin_step = step * begin_part;
		end_step = step * end_part;

		restore_from_prev(ad, tr_id, th_id, step);

		rk_int_lndham(ad, tr_id, th_id, begin_step);

		if (jump > 0 && ed->is_obs == 1)
		{
			double jump_time = ed->times_all[tr_id];
			double jump_norm = norm_square(phi, sys_size);
			double jump_eta = *eta;

			ed->jump_times[tr_id].push_back(jump_time);
			ed->jump_norms[tr_id].push_back(jump_norm);
			ed->jump_etas[tr_id].push_back(jump_eta);

			ed->jumps_counts[tr_id]++;
		}

		rk_recovery(ad, tr_id, th_id);
		vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, *stream, 1, eta, 0.0, 1.0);
		while (*eta == 0.0)
		{
			vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, *stream, 1, eta, 0.0, 1.0);
		}

		rk_step_lndham(ad, tr_id, th_id, end_step);
	}
}

void rk_step_integrable(AllData* ad, int tr_id, int th_id, double step)
{
	RunParam* rp = ad->rp;
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int sys_size = md->sys_size;

	int jump = int(cp->params.find("jump")->second);

	VSLStreamStatePtr* stream = &(ed->streams[tr_id]);
	MKL_Complex16* phi = &(ed->phi_all[tr_id * sys_size]);
	double* eta = &(ed->etas_all[tr_id]);

	double prev_norm = 0.0;
	double curr_norm = 0.0;
	double norm_diff = 0.0;
	double begin_part = 0.0;
	double end_part = 0.0;
	double begin_step = 0.0;
	double end_step = 0.0;

	save_phi_prev(ad, tr_id, th_id);
	prev_norm = norm_square(phi, sys_size);
	rk_int_integrable(ad, tr_id, th_id, step);

	if (is_norm_crossed(phi, eta, sys_size))
	{
		curr_norm = norm_square(phi, sys_size);
		norm_diff = prev_norm - curr_norm;
		begin_part = (prev_norm - *(eta)) / norm_diff;
		end_part = 1.0 - begin_part;
		begin_step = step * begin_part;
		end_step = step * end_part;

		restore_from_prev(ad, tr_id, th_id, step);

		rk_int_integrable(ad, tr_id, th_id, begin_step);

		if (jump > 0 && ed->is_obs == 1)
		{
			double jump_time = ed->times_all[tr_id];
			double jump_norm = norm_square(phi, sys_size);
			double jump_eta = *eta;

			ed->jump_times[tr_id].push_back(jump_time);
			ed->jump_norms[tr_id].push_back(jump_norm);
			ed->jump_etas[tr_id].push_back(jump_eta);

			ed->jumps_counts[tr_id]++;
		}

		rk_recovery(ad, tr_id, th_id);
		vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, *stream, 1, eta, 0.0, 1.0);
		while (*eta == 0.0)
		{
			vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, *stream, 1, eta, 0.0, 1.0);
		}

		rk_step_integrable(ad, tr_id, th_id, end_step);
	}
}

void rk_step_floq2spins(AllData* ad, int tr_id, int th_id, double step)
{
	RunParam* rp = ad->rp;
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int sys_size = md->sys_size;

	int jump = int(cp->params.find("jump")->second);

	VSLStreamStatePtr* stream = &(ed->streams[tr_id]);
	MKL_Complex16* phi = &(ed->phi_all[tr_id * sys_size]);
	double* eta = &(ed->etas_all[tr_id]);

	double prev_norm = 0.0;
	double curr_norm = 0.0;
	double norm_diff = 0.0;
	double begin_part = 0.0;
	double end_part = 0.0;
	double begin_step = 0.0;
	double end_step = 0.0;

	save_phi_prev(ad, tr_id, th_id);
	prev_norm = norm_square(phi, sys_size);
	rk_int_floq2spins(ad, tr_id, th_id, step);

	if (is_norm_crossed(phi, eta, sys_size))
	{
		curr_norm = norm_square(phi, sys_size);
		norm_diff = prev_norm - curr_norm;
		begin_part = (prev_norm - *(eta)) / norm_diff;
		end_part = 1.0 - begin_part;
		begin_step = step * begin_part;
		end_step = step * end_part;

		restore_from_prev(ad, tr_id, th_id, step);

		rk_int_floq2spins(ad, tr_id, th_id, begin_step);

		if (jump > 0 && ed->is_obs == 1)
		{
			double jump_time = ed->times_all[tr_id];
			double jump_norm = norm_square(phi, sys_size);
			double jump_eta = *eta;

			ed->jump_times[tr_id].push_back(jump_time);
			ed->jump_norms[tr_id].push_back(jump_norm);
			ed->jump_etas[tr_id].push_back(jump_eta);

			ed->jumps_counts[tr_id]++;
		}

		rk_recovery(ad, tr_id, th_id);
		vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, *stream, 1, eta, 0.0, 1.0);
		while (*eta == 0.0)
		{
			vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, *stream, 1, eta, 0.0, 1.0);
		}

		rk_step_floq2spins(ad, tr_id, th_id, end_step);
	}
}

void rk_step_floqnspins(AllData* ad, int tr_id, int th_id, double step)
{
	RunParam* rp = ad->rp;
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int sys_size = md->sys_size;

	int jump = int(cp->params.find("jump")->second);

	VSLStreamStatePtr* stream = &(ed->streams[tr_id]);
	MKL_Complex16* phi = &(ed->phi_all[tr_id * sys_size]);
	double* eta = &(ed->etas_all[tr_id]);

	double prev_norm = 0.0;
	double curr_norm = 0.0;
	double norm_diff = 0.0;
	double begin_part = 0.0;
	double end_part = 0.0;
	double begin_step = 0.0;
	double end_step = 0.0;

	save_phi_prev(ad, tr_id, th_id);
	prev_norm = norm_square(phi, sys_size);
	rk_int_floqnspins(ad, tr_id, th_id, step);

	if (is_norm_crossed(phi, eta, sys_size))
	{
		curr_norm = norm_square(phi, sys_size);
		norm_diff = prev_norm - curr_norm;
		begin_part = (prev_norm - *(eta)) / norm_diff;
		end_part = 1.0 - begin_part;
		begin_step = step * begin_part;
		end_step = step * end_part;

		restore_from_prev(ad, tr_id, th_id, step);

		rk_int_floqnspins(ad, tr_id, th_id, begin_step);

		if (jump > 0 && ed->is_obs == 1)
		{
			double jump_time = ed->times_all[tr_id];
			double jump_norm = norm_square(phi, sys_size);
			double jump_eta = *eta;

			ed->jump_times[tr_id].push_back(jump_time);
			ed->jump_norms[tr_id].push_back(jump_norm);
			ed->jump_etas[tr_id].push_back(jump_eta);

			ed->jumps_counts[tr_id]++;
		}

		rk_recovery(ad, tr_id, th_id);
		vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, *stream, 1, eta, 0.0, 1.0);
		while (*eta == 0.0)
		{
			vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, *stream, 1, eta, 0.0, 1.0);
		}

		rk_step_floqnspins(ad, tr_id, th_id, end_step);
	}
}

void rk_step_flq_sp_ph(AllData* ad, int tr_id, int th_id, double step)
{
	RunParam* rp = ad->rp;
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;
	ExpData* ed = ad->ed;

	int sys_size = md->sys_size;

	int jump = int(cp->params.find("jump")->second);

	VSLStreamStatePtr* stream = &(ed->streams[tr_id]);
	MKL_Complex16* phi = &(ed->phi_all[tr_id * sys_size]);
	double* eta = &(ed->etas_all[tr_id]);

	double prev_norm = 0.0;
	double curr_norm = 0.0;
	double norm_diff = 0.0;
	double begin_part = 0.0;
	double end_part = 0.0;
	double begin_step = 0.0;
	double end_step = 0.0;

	save_phi_prev(ad, tr_id, th_id);
	prev_norm = norm_square(phi, sys_size);
	rk_int_flq_sp_ph(ad, tr_id, th_id, step);

	if (is_norm_crossed(phi, eta, sys_size))
	{
		curr_norm = norm_square(phi, sys_size);
		norm_diff = prev_norm - curr_norm;
		begin_part = (prev_norm - *(eta)) / norm_diff;
		end_part = 1.0 - begin_part;
		begin_step = step * begin_part;
		end_step = step * end_part;

		restore_from_prev(ad, tr_id, th_id, step);

		rk_int_flq_sp_ph(ad, tr_id, th_id, begin_step);

		if (jump > 0 && ed->is_obs == 1)
		{
			double jump_time = ed->times_all[tr_id];
			double jump_norm = norm_square(phi, sys_size);
			double jump_eta = *eta;

			ed->jump_times[tr_id].push_back(jump_time);
			ed->jump_norms[tr_id].push_back(jump_norm);
			ed->jump_etas[tr_id].push_back(jump_eta);

			ed->jumps_counts[tr_id]++;
		}

		rk_recovery(ad, tr_id, th_id);
		vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, *stream, 1, eta, 0.0, 1.0);
		while (*eta == 0.0)
		{
			vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, *stream, 1, eta, 0.0, 1.0);
		}

		rk_step_flq_sp_ph(ad, tr_id, th_id, end_step);
	}
}


void calc_ci_double(AllData * ad, int tr_id)
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	double eps = double(cp->params.find("cd_eps")->second);

	double * curr_diff = new double[ed->cd_dim];
	double curr_norm = 0.0;
	double integral = 0.0;

	for (int cd_p_id_1 = 0; cd_p_id_1 < ed->cd_num_points; cd_p_id_1++)
	{
		for (int cd_p_id_2 = 0; cd_p_id_2 < ed->cd_num_points; cd_p_id_2++)
		{
			if (cd_p_id_1 != cd_p_id_2)
			{
				for (int cd_st_id = 0; cd_st_id < ed->cd_dim; cd_st_id++)
				{
					curr_diff[cd_st_id] = ed->cd_rec_data[tr_id][cd_p_id_1][cd_st_id] - ed->cd_rec_data[tr_id][cd_p_id_2][cd_st_id];
				}

				curr_norm = get_norm_cd(curr_diff, ed->cd_dim) / double(md->sys_size);

				if (curr_norm < eps)
				{
					integral += 1.0;
				}
			}
		}
	}

	integral /= (double(ed->cd_num_points) * double(ed->cd_num_points - 1));

	ed->cd_i[tr_id] = integral;

	delete[] curr_diff;
}

void dump_phi(AllData * ad)
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	ExpData * ed = ad->ed;
	MainData * md = ad->md;

	int dump_phi = int(cp->params.find("dump_phi")->second);

	if (dump_phi == 1)
	{
		string fn;

		int sys_size = md->sys_size;
		int num_trajectories = cp->num_trajectories;

		for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
		{
			fn = rp->path + "phi_" + to_string(tr_id) + cp->fn_suffix;
			save_complex_data(fn, &(ed->phi_all[tr_id * sys_size]), sys_size, 16, false);
		}
	}
}

void dump_phi_evo(AllData * ad, bool append)
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	ExpData * ed = ad->ed;
	MainData * md = ad->md;

	int dump_phi_evo = int(cp->params.find("dump_phi_evo")->second);

	if (dump_phi_evo == 1)
	{
		string fn;

		int sys_size = md->sys_size;
		int num_trajectories = cp->num_trajectories;

		for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
		{
			fn = rp->path + "phi_evo_" + to_string(tr_id) + cp->fn_suffix;
			save_complex_data(fn, &(ed->phi_all[tr_id * sys_size]), sys_size, 16, append);
		}
	}
}


void calc_random_obs(AllData * ad, int tr_id)
{
	int num_random_obs = int(ad->cp->params.find("num_random_obs")->second);
	if (num_random_obs > 0)
	{
		vector<complex<double>> random_obs = get_random_obs(ad, tr_id);
		ad->ed->random_obs[tr_id] = random_obs;
	}
}

void calc_random_obs_lpn_start(AllData * ad, int tr_id)
{
	int num_random_obs = int(ad->cp->params.find("num_random_obs")->second);
	if (num_random_obs > 0)
	{
		ad->ed->random_obs_lpn[tr_id] = ad->ed->random_obs[tr_id];
	}
}

void calc_random_obs_lpn(AllData * ad, int tr_id)
{
	int num_random_obs = int(ad->cp->params.find("num_random_obs")->second);
	if (num_random_obs > 0)
	{
		vector<complex<double>> random_obs = get_random_obs(ad, tr_id);
		ad->ed->random_obs_lpn[tr_id] = random_obs;
	}
}