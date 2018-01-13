#include "dump.h"

void dump_adr_single(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int tr_id)
{
	int is_adr_dump_sep = int(cp->params.find("is_adr_dump_sep")->second);

	if (is_adr_dump_sep == 1)
	{
		int sys_size = md->sys_size;
		double * adr = &(qjd->abs_diag_rho_all[tr_id * sys_size]);

		string fn = rp->path + "_adr_" + to_string(tr_id) + cp->fn_suffix;
		save_double_data(fn, adr, sys_size, 16, false);
	}
}

void dump_adr_avg(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd)
{
	int is_adr_dump_avg = int(cp->params.find("is_adr_dump_avg")->second);

	if (is_adr_dump_avg == 1)
	{
		int sys_size = md->sys_size;
		int num_trajectories = cp->qj_num_trajectories;

		double * adr_avg = new double[sys_size];

		for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
		{
			double * adr = &(qjd->abs_diag_rho_all[tr_id * sys_size]);
			for (int st_id = 0; st_id < sys_size; st_id++)
			{
				adr_avg[st_id] += adr[st_id] / double(num_trajectories);
			}
		}

		string fn = rp->path + "_adr_avg" + cp->fn_suffix;
		save_double_data(fn, adr_avg, sys_size, 16, false);

		delete[] adr_avg;
	}
}

void update_evo(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int dump_id)
{
	int num_trajectories = cp->qj_num_trajectories;
	int num_dumps_total = qjd->num_dumps_total;

	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		qjd->mean_evo[tr_id * num_dumps_total + dump_id] = qjd->mean[tr_id];
		qjd->dispersion_evo[tr_id * num_dumps_total + dump_id] = qjd->dispersion[tr_id];
		qjd->m2_evo[tr_id * num_dumps_total + dump_id] = qjd->m2[tr_id];
	}
}

void update_lpn_evo(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int dump_id)
{
	int num_trajectories = cp->qj_num_trajectories;
	int num_dumps_total = qjd->num_dumps_total;

	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		qjd->energy_evo[tr_id * num_dumps_total + dump_id] = qjd->energy[tr_id];
		qjd->lambda_evo[tr_id * num_dumps_total + dump_id] = qjd->lambda_now[tr_id];
		qjd->mean_lpn_evo[tr_id * num_dumps_total + dump_id] = qjd->mean_lpn[tr_id];
		qjd->energy_lpn_evo[tr_id * num_dumps_total + dump_id] = qjd->energy_lpn[tr_id];
	}
}

void dump_evo(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd)