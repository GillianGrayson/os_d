#include "dump.h"

void dump_adr_single(AllData * ad, int tr_id, bool append)
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);
	int dump_adr_sep = int(cp->params.find("dump_adr_sep")->second);

	if(dump_evo_sep == 1)
	{
		if (dump_adr_sep == 1)
		{
			int sys_size = md->sys_size;
			double * adr = &(ed->abs_diag_rho_all[tr_id * sys_size]);

			string fn = rp->path + "adr_" + to_string(tr_id) + cp->fn_suffix;
			save_double_data(fn, adr, sys_size, 16, append);
		}
	}
}

void dump_adr_avg(AllData * ad, bool append)
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int dump_adr_avg = int(cp->params.find("dump_adr_avg")->second);

	if (dump_adr_avg == 1)
	{
		int sys_size = md->sys_size;
		int num_trajectories = cp->num_trajectories;

		double * adr_avg = new double[sys_size];

		for (int st_id = 0; st_id < sys_size; st_id++)
		{
			adr_avg[st_id] = 0.0;
		}

		for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
		{
			double * adr = &(ed->abs_diag_rho_all[tr_id * sys_size]);
			for (int st_id = 0; st_id < sys_size; st_id++)
			{
				adr_avg[st_id] += adr[st_id] / double(num_trajectories);
			}
		}

		string fn = rp->path + "adr_avg" + cp->fn_suffix;
		save_double_data(fn, adr_avg, sys_size, 16, append);

		delete[] adr_avg;
	}
}

void dump_adr_avg_mult(AllData * ad, bool append, int begin_traj_id, int end_traj_id)
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int dump_adr_avg = int(cp->params.find("dump_adr_avg")->second);

	if (dump_adr_avg == 1)
	{
		int sys_size = md->sys_size;
		int num_trajectories = end_traj_id - begin_traj_id;

		double * adr_avg = new double[sys_size];

		for (int st_id = 0; st_id < sys_size; st_id++)
		{
			adr_avg[st_id] = 0.0;
		}

		for (int tr_id = begin_traj_id; tr_id < end_traj_id; tr_id++)
		{
			double * adr = &(ed->abs_diag_rho_all[tr_id * sys_size]);
			for (int st_id = 0; st_id < sys_size; st_id++)
			{
				adr_avg[st_id] += adr[st_id] / double(num_trajectories);
			}
		}

		string fn = rp->path + "adr_avg" + cp->fn_suffix;
		save_double_data(fn, adr_avg, sys_size, 16, append);

		delete[] adr_avg;
	}
}

void dump_cd(AllData * ad)
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	ExpData * ed = ad->ed;

	int num_trajectories = cp->num_trajectories;

	double * ci = ed->cd_i;

	string fn;

	fn = rp->path + "ci" + cp->fn_suffix;
	save_double_data(fn, ci, num_trajectories, 16, false);
}
