#include "qj_initiator.h"
#include "split_proc.h"

void LyapunovMCBehaviour::init_data(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd) const
{
	init_splits(rp, cp, md, qjd);
	init_streams(rp, cp, md, qjd);
	init_streams_var(rp, cp, md, qjd);
}

void init_splits(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd)
{
	int num_threads = rp->num_threads;

	Split * head = md->structure;

	head = init_split_structure(rp, cp, md);

	Split * heads = new Split[num_threads];
	for (int th_id = 0; th_id < num_threads; th_id++)
	{
		copy_struct_not_member(head, &heads[th_id]);
	}
}

void init_streams(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd)
{
	int num_trajectories = cp->qj_num_trajectories;
	int seed = cp->qj_seed;
	int mns = cp->qj_mns;

	VSLStreamStatePtr * streams = qjd->streams;

	streams = new VSLStreamStatePtr[num_trajectories];
	vslNewStream(&streams[0], VSL_BRNG_MCG59, 777);
	vslLeapfrogStream(streams[0], seed, mns);
}

void init_streams_var(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd)
{
	int num_trajectories = cp->qj_num_trajectories;
	int seed = cp->qj_seed;
	int mns = cp->qj_mns;

	VSLStreamStatePtr * streams_var = qjd->streams_var;

	VSLStreamStatePtr * streams_var = new VSLStreamStatePtr[num_trajectories];
	vslNewStream(&streams_var[0], VSL_BRNG_MCG31, 777);
	for (int tr_id = 1; tr_id < num_trajectories; tr_id++)
	{
		vslCopyStream(&streams_var[tr_id], streams_var[0]);
	}
	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		vslLeapfrogStream(streams_var[tr_id], tr_id, num_trajectories);
	}
}

void init_basic_data(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd)
{
	int num_trajectories = cp->qj_num_trajectories;
	int N = cp->params.find("N")->second;

	qjd->phi_all = new MKL_Complex16[num_trajectories * N];
	qjd->abs_diag_rho_all = new double[num_trajectories * N];
	qjd->eta_all = new double[num_trajectories];

	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		for (int st_id = 0; st_id < N; st_id++)
		{
			qjd->phi_all[tr_id * N + st_id].real = 0.0;
			qjd->phi_all[tr_id * N + st_id].imag = 0.0;

			qjd->abs_diag_rho_all[tr_id * N + st_id] = 0.0;
			qjd->abs_diag_rho_all[tr_id * N + st_id] = 0.0;
		}

		qjd->eta_all[tr_id] = 0.0;
	}
}

void init_obs_std(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd)
{
	int num_trajectories = cp->qj_num_trajectories;

	qjd->mean_start = new double[num_trajectories];
	qjd->mean = new double[num_trajectories];
	qjd->dispersion = new double[num_trajectories];
	qjd->m2 = new double[num_trajectories];

	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		qjd->mean_start[tr_id] = 0.0;
		qjd->mean[tr_id] = 0.0;
		qjd->dispersion[tr_id] = 0.0;
		qjd->m2[tr_id] = 0.0;
	}
}

void init_obs_lpn(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd)
{
	int num_trajectories = cp->qj_num_trajectories;

	qjd->energy = new double[num_trajectories];
	qjd->lambda = new double[num_trajectories];
	qjd->delta_s = new double[num_trajectories];

	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		qjd->energy[tr_id] = 0.0;
		qjd->lambda[tr_id] = 0.0;
		qjd->delta_s[tr_id] = 0.0;
	}
}

void init_obs_lpn_evo(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd)
{
	int num_trajectories = cp->qj_num_trajectories;

	qjd->energy = new double[num_trajectories];
	qjd->lambda = new double[num_trajectories];
	qjd->delta_s = new double[num_trajectories];

	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		qjd->energy[tr_id] = 0.0;
		qjd->lambda[tr_id] = 0.0;
		qjd->delta_s[tr_id] = 0.0;
	}
}


