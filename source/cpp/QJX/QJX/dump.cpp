#include "dump.h"

void dump_adr_single(AllData * ad, int tr_id, bool append)
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int is_evo_dump_sep = int(cp->params.find("is_evo_dump_sep")->second);
	int is_adr_dump_sep = int(cp->params.find("is_adr_dump_sep")->second);

	if(is_evo_dump_sep == 1)
	{
		if (is_adr_dump_sep == 1)
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

	int is_adr_dump_avg = int(cp->params.find("is_adr_dump_avg")->second);

	if (is_adr_dump_avg == 1)
	{
		int sys_size = md->sys_size;
		int num_trajectories = cp->qj_num_trajectories;

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

void update_evo_std(AllData * ad, int dump_id)
{
	ConfigParam * cp = ad->cp;
	ExpData * ed = ad->ed;

	int num_trajectories = cp->qj_num_trajectories;
	int num_dumps_total = ed->num_dumps_total;

	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		ed->mean_evo[tr_id * num_dumps_total + dump_id] = ed->mean[tr_id];
		ed->dispersion_evo[tr_id * num_dumps_total + dump_id] = ed->dispersion[tr_id];
		ed->m2_evo[tr_id * num_dumps_total + dump_id] = ed->m2[tr_id];
	}
}

void update_evo_lpn(AllData * ad, int dump_id)
{
	ConfigParam * cp = ad->cp;
	ExpData * ed = ad->ed;

	int num_trajectories = cp->qj_num_trajectories;
	int num_dumps_total = ed->num_dumps_total;

	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		ed->energy_evo[tr_id * num_dumps_total + dump_id] = ed->energy[tr_id];
		ed->lambda_evo[tr_id * num_dumps_total + dump_id] = ed->lambda_now[tr_id];
		ed->mean_lpn_evo[tr_id * num_dumps_total + dump_id] = ed->mean_lpn[tr_id];
		ed->energy_lpn_evo[tr_id * num_dumps_total + dump_id] = ed->energy_lpn[tr_id];
	}
}

void dump_std(AllData * ad)
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	ExpData * ed = ad->ed;

	int is_obs_dump = int(cp->params.find("is_obs_dump")->second);

	if (is_obs_dump == 1)
	{
		int num_trajectories = cp->qj_num_trajectories;

		double * mean_start = ed->mean_start;
		double * mean = ed->mean;
		double * dispersion = ed->dispersion;
		double * m2 = ed->m2;

		string fn;

		fn = rp->path + "mean_start" + cp->fn_suffix;
		save_double_data(fn, mean_start, num_trajectories, 16, false);

		fn = rp->path + "mean" + cp->fn_suffix;
		save_double_data(fn, mean, num_trajectories, 16, false);

		fn = rp->path + "dispersion" + cp->fn_suffix;
		save_double_data(fn, dispersion, num_trajectories, 16, false);

		fn = rp->path + "m2" + cp->fn_suffix;
		save_double_data(fn, m2, num_trajectories, 16, false);
	}
}

void dump_lpn(AllData * ad)
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	ExpData * ed = ad->ed;

	int num_trajectories = cp->qj_num_trajectories;

	int is_obs_dump = int(cp->params.find("is_obs_dump")->second);

	if (is_obs_dump == 1)
	{

		double * energy = ed->energy;
		double * lambda = ed->lambda_now;
		double * mean_lpn = ed->mean_lpn;
		double * energy_lpn = ed->energy_lpn;

		string fn;

		fn = rp->path + "energy" + cp->fn_suffix;
		save_double_data(fn, energy, num_trajectories, 16, false);

		fn = rp->path + "lambda" + cp->fn_suffix;
		save_double_data(fn, lambda, num_trajectories, 16, false);

		fn = rp->path + "mean_lpn" + cp->fn_suffix;
		save_double_data(fn, mean_lpn, num_trajectories, 16, false);

		fn = rp->path + "energy_lpn" + cp->fn_suffix;
		save_double_data(fn, energy_lpn, num_trajectories, 16, false);
	}
}

void dump_cd(AllData * ad)
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	ExpData * ed = ad->ed;

	int num_trajectories = cp->qj_num_trajectories;

	double * ci = ed->cd_i;

	string fn;

	fn = rp->path + "ci" + cp->fn_suffix;
	save_double_data(fn, ci, num_trajectories, 16, false);
}

void dump_evo_std(AllData * ad)
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	ExpData * ed = ad->ed;

	int num_trajectories = cp->qj_num_trajectories;
	int num_dumps_total = ed->num_dumps_total;

	int is_obs_dump = int(cp->params.find("is_obs_dump")->second);

	if (is_obs_dump == 1)
	{

		int * dump_periods = ed->dump_periods;

		double * mean_evo = ed->mean_evo;
		double * dispersion_evo = ed->dispersion_evo;
		double * m2_evo = ed->m2_evo;

		string fn;

		fn = rp->path + "periods" + cp->fn_suffix;
		save_int_data(fn, dump_periods, num_dumps_total, false);

		fn = rp->path + "mean_evo" + cp->fn_suffix;
		save_2d_inv_double_data(fn, mean_evo, num_dumps_total, num_trajectories, 16, false);

		fn = rp->path + "dispersion_evo" + cp->fn_suffix;
		save_2d_inv_double_data(fn, dispersion_evo, num_dumps_total, num_trajectories, 16, false);

		fn = rp->path + "m2_evo" + cp->fn_suffix;
		save_2d_inv_double_data(fn, m2_evo, num_dumps_total, num_trajectories, 16, false);
	}
}

void dump_evo_lpn(AllData * ad)
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	ExpData * ed = ad->ed;

	int num_trajectories = cp->qj_num_trajectories;
	int num_dumps_total = ed->num_dumps_total;

	int is_obs_dump = int(cp->params.find("is_obs_dump")->second);

	if (is_obs_dump == 1)
	{

		double * energy_evo = ed->energy_evo;
		double * lambda_evo = ed->lambda_evo;
		double * mean_lpn_evo = ed->mean_lpn_evo;
		double * energy_lpn_evo = ed->energy_lpn_evo;

		string fn;

		fn = rp->path + "energy_evo" + cp->fn_suffix;
		save_2d_inv_double_data(fn, energy_evo, num_dumps_total, num_trajectories, 16, false);

		fn = rp->path + "lambda_evo" + cp->fn_suffix;
		save_2d_inv_double_data(fn, lambda_evo, num_dumps_total, num_trajectories, 16, false);

		fn = rp->path + "mean_lpn_evo" + cp->fn_suffix;
		save_2d_inv_double_data(fn, mean_lpn_evo, num_dumps_total, num_trajectories, 16, false);

		fn = rp->path + "energy_lpn_evo" + cp->fn_suffix;
		save_2d_inv_double_data(fn, energy_lpn_evo, num_dumps_total, num_trajectories, 16, false);
	}
}