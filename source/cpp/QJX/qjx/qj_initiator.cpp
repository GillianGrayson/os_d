#include "qj_initiator.h"
#include "split_proc.h"

void LyapunovMCInitBehaviour::init_data(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd) const
{
	init_splits(rp, cp, md, qjd);
	init_streams(rp, cp, md, qjd);
	init_streams_var(rp, cp, md, qjd);
	init_basic_data(rp, cp, md, qjd);
	init_dump_priods(rp, cp, md, qjd);
	init_obs_std(rp, cp, md, qjd);
	init_obs_lpn(rp, cp, md, qjd);
	init_obs_lpn_evo(rp, cp, md, qjd);

	init_start_state(rp, cp, md, qjd, 0);
}

void init_splits(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd)
{
	int num_threads = rp->num_threads;

	Split * structure = md->structure;
	Split * splits = md->splits;

	structure = init_split_structure(rp, cp, md);

	splits = new Split[num_threads];
	for (int th_id = 0; th_id < num_threads; th_id++)
	{
		copy_struct_not_member(structure, &splits[th_id]);
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
	int sys_size = md->sys_size;

	qjd->phi_all = new MKL_Complex16[num_trajectories * sys_size];
	qjd->phi_all_aux = new MKL_Complex16[num_trajectories * sys_size];
	qjd->abs_diag_rho_all = new double[num_trajectories * sys_size];
	qjd->times_all = new double[num_trajectories];
	qjd->etas_all = new double[num_trajectories];

	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		for (int st_id = 0; st_id < sys_size; st_id++)
		{
			qjd->phi_all[tr_id * sys_size + st_id].real = 0.0;
			qjd->phi_all[tr_id * sys_size + st_id].imag = 0.0;

			qjd->phi_all_aux[tr_id * sys_size + st_id].real = 0.0;
			qjd->phi_all_aux[tr_id * sys_size + st_id].imag = 0.0;

			qjd->abs_diag_rho_all[tr_id * sys_size + st_id] = 0.0;
			qjd->abs_diag_rho_all[tr_id * sys_size + st_id] = 0.0;
		}

		qjd->times_all[tr_id] = 0.0;
		qjd->etas_all[tr_id] = 0.0;
	}
}

void init_dump_priods(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd)
{
	qjd->dump_type = int(cp->params.find("dump_type")->second);
	int num_dumps = int(cp->params.find("num_dumps")->second);
	int num_obs_periods = cp->qj_num_obs_periods;

	if (qjd->dump_type == 0)
	{
		qjd->num_dumps_total = num_dumps + 1;
		qjd->dump_periods = new int[qjd->num_dumps_total];

		qjd->dump_periods[0] = 0;

		int dump_shift = num_obs_periods / num_dumps;
		for (int dump_id = 0; dump_id < num_dumps; dump_id++)
		{
			qjd->dump_periods[dump_id + 1] = (dump_id + 1) * dump_shift;
		}
	}
	else if (qjd->dump_type == 1)
	{
		qjd->num_dumps_total = num_dumps + 2;
		qjd->dump_periods = new int[qjd->num_dumps_total];

		qjd->dump_periods = 0;

		double begin_decade = log10(1.0);
		double end_decade = log10(double(num_obs_periods));
		double num_decades = end_decade - begin_decade;
		double num_decades_dump = double(num_dumps) / num_decades;
		for (int dump_id = 0; dump_id < num_dumps + 1; dump_id++)
		{
			int curr_val = int(pow(10.0, begin_decade) * pow(10.0, (1.0 / num_decades_dump) * double(dump_id)));

			if (curr_val > qjd->dump_periods[dump_id])
			{
				qjd->dump_periods[dump_id + 1] = curr_val;
			}
			else
			{
				qjd->dump_periods[dump_id + 1] = qjd->dump_periods[dump_id] + 1;
			}
		}
	}
}

void init_obs_std(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd)
{
	int num_trajectories = cp->qj_num_trajectories;
	int num_dumps_total = qjd->num_dumps_total;

	qjd->mean_start		= new double[num_trajectories];

	qjd->mean			= new double[num_trajectories];
	qjd->dispersion		= new double[num_trajectories];
	qjd->m2				= new double[num_trajectories];

	qjd->mean_evo		= new double[num_trajectories * num_dumps_total];
	qjd->dispersion_evo = new double[num_trajectories * num_dumps_total];
	qjd->m2_evo			= new double[num_trajectories * num_dumps_total];

	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		qjd->mean_start[tr_id]	= 0.0;

		qjd->mean[tr_id]		= 0.0;
		qjd->dispersion[tr_id]	= 0.0;
		qjd->m2[tr_id]			= 0.0;

		for (int dump_id = 0; dump_id < num_dumps_total; dump_id++)
		{
			qjd->mean_evo[tr_id * num_dumps_total + dump_id]		= 0.0;
			qjd->dispersion_evo[tr_id * num_dumps_total + dump_id]	= 0.0;
			qjd->m2_evo[tr_id * num_dumps_total + dump_id]			= 0.0;
		}
	}
}

void init_obs_lpn(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd)
{
	int num_trajectories = cp->qj_num_trajectories;
	int num_dumps_total = qjd->num_dumps_total;

	qjd->delta_s = new double[num_trajectories];

	qjd->energy = new double[num_trajectories];
	qjd->lambda = new double[num_trajectories];
	qjd->mean_lpn = new double[num_trajectories];
	qjd->energy_lpn = new double[num_trajectories];
	qjd->lambda_lpn = new double[num_trajectories];

	qjd->energy_evo = new double[num_trajectories * num_dumps_total];
	qjd->lambda_evo = new double[num_trajectories * num_dumps_total];
	qjd->mean_lpn_evo = new double[num_trajectories * num_dumps_total];
	qjd->energy_lpn_evo = new double[num_trajectories * num_dumps_total];
	qjd->lambda_lpn_evo = new double[num_trajectories * num_dumps_total];

	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		qjd->delta_s[tr_id] = 0.0;

		qjd->energy[tr_id] = 0.0;
		qjd->lambda[tr_id] = 0.0;
		qjd->mean_lpn[tr_id] = 0.0;
		qjd->energy_lpn[tr_id] = 0.0;
		qjd->lambda_lpn[tr_id] = 0.0;

		for (int dump_id = 0; dump_id < num_dumps_total; dump_id++)
		{
			qjd->energy_evo[tr_id * num_dumps_total + dump_id] = 0.0;
			qjd->lambda_evo[tr_id * num_dumps_total + dump_id] = 0.0;
			qjd->mean_lpn_evo[tr_id * num_dumps_total + dump_id] = 0.0;
			qjd->energy_lpn_evo[tr_id * num_dumps_total + dump_id] = 0.0;
			qjd->lambda_lpn_evo[tr_id * num_dumps_total + dump_id] = 0.0;
		}	
	}
}

void init_start_state(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int tr_id)
{
	int sys_size = md->sys_size;
	int start_type = int(cp->params.find("start_type")->second);
	int start_state = int(cp->params.find("start_state")->second);

	MKL_Complex16 * phi = &(qjd->phi_all[tr_id * sys_size]);

	if (start_type == 0)
	{
		phi[start_state].real = 1.0;
	}
	else
	{
		for (int i = 0; i < sys_size; i++)
		{
			phi[i].real = sqrt(1.0 / double(sys_size));
		}
	}
}

