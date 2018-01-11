#include "qj_destructor.h"
#include "split_proc.h"

void LyapunovMCFreeBehaviour::free_data(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd) const
{
	free_splits(rp, cp, md, qjd);
	free_streams(rp, cp, md, qjd);
	free_streams_var(rp, cp, md, qjd);
	free_basic_data(rp, cp, md, qjd);
	free_dump_priods(rp, cp, md, qjd);
	free_obs_std(rp, cp, md, qjd);
	free_obs_lpn(rp, cp, md, qjd);
}

void free_splits(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd)
{
	int num_threads = rp->num_threads;

	Split * head = md->structure;
	Split * heads = md->splits;

	for (int i = 0; i < num_threads; i++)
	{
		delete_split_struct_not_member(&heads[i]);
	}
	delete(heads);
	delete_split_struct(head);
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
	delete[] qjd->mean_lpn;
	delete[] qjd->energy_lpn;
	delete[] qjd->lambda_lpn;

	delete[] qjd->energy_evo;
	delete[] qjd->lambda_evo;
	delete[] qjd->mean_lpn_evo;
	delete[] qjd->energy_lpn_evo;
	delete[] qjd->lambda_lpn_evo;
}
