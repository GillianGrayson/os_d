#include "qj_destructor.h"
#include "split_proc.h"

void LpnFreeBehaviour::free_data(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd) const
{
	free_splits(rp, cp, md, qjd);
	free_streams(rp, cp, md, qjd);
	free_streams_var(rp, cp, md, qjd);
	free_basic_data(rp, cp, md, qjd);
	free_dump_priods(rp, cp, md, qjd);
	free_obs_std(rp, cp, md, qjd);
	free_obs_lpn(rp, cp, md, qjd);
}

void StdFreeBehaviour::free_data(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd) const
{
	free_splits(rp, cp, md, qjd);
	free_streams(rp, cp, md, qjd);
	free_basic_data(rp, cp, md, qjd);
	free_dump_priods(rp, cp, md, qjd);
	free_obs_std(rp, cp, md, qjd);
}

void CorrDimFreeBehaviour::free_data(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd) const
{
	free_splits_deep(rp, cp, md, qjd);
	free_streams(rp, cp, md, qjd);
	free_basic_data(rp, cp, md, qjd);
	free_dump_priods(rp, cp, md, qjd);
	free_obs_std(rp, cp, md, qjd);
	free_obs_cd(rp, cp, md, qjd);
}

void free_splits_deep(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd)
{
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

void free_splits(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd)
{
	int num_threads = rp->num_threads;

	for (int i = 0; i < num_threads; i++)
	{
		delete_split_struct_not_member(&(md->splits[i]));
	}

	delete(md->splits);
	delete_split_struct(md->structure);
}

void free_streams(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd)
{
	VSLStreamStatePtr * streams = qjd->streams;
	delete[] streams;
}

void free_streams_var(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd)
{
	VSLStreamStatePtr * streams_var = qjd->streams_var;
	delete[] streams_var;
}

void free_basic_data(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd)
{
	delete[] qjd->phi_all;
	delete[] qjd->phi_all_aux;
	delete[] qjd->abs_diag_rho_all;
	delete[] qjd->times_all;
	delete[] qjd->etas_all;
}

void free_dump_priods(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd)
{
	delete[] qjd->dump_periods;
}

void free_obs_std(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd)
{
	delete[] qjd->mean_start;

	delete[] qjd->mean;
	delete[] qjd->dispersion;
	delete[] qjd->m2;

	delete[] qjd->mean_evo;
	delete[] qjd->dispersion_evo;
	delete[] qjd->m2_evo;
}

void free_obs_lpn(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd)
{
	delete[] qjd->delta_s;

	delete[] qjd->energy;
	delete[] qjd->lambda;
	delete[] qjd->lambda_now;
	delete[] qjd->mean_lpn;
	delete[] qjd->energy_lpn;

	delete[] qjd->energy_evo;
	delete[] qjd->lambda_evo;
	delete[] qjd->mean_lpn_evo;
	delete[] qjd->energy_lpn_evo;
}

void free_obs_cd(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd)
{
	int num_trajectories = cp->qj_num_trajectories;

	delete[] qjd->cd_i;

	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		for (int p_id = 0; p_id < qjd->cd_num_points; p_id++)
		{
			delete[] qjd->cd_rec_data[tr_id][p_id];
		}
		delete[] qjd->cd_rec_data[tr_id];
	}
	delete[] qjd->cd_rec_data;
}