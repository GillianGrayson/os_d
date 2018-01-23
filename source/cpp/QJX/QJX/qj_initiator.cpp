#include "qj_initiator.h"
#include "split_proc.h"

void LpnInitBehaviour::init_data(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd) const
{
	init_splits(rp, cp, md, qjd);
	init_streams(rp, cp, md, qjd);
	leap_frog_single_stream(rp, cp, md, qjd, 0);
	init_streams_var(rp, cp, md, qjd);
	init_basic_data(rp, cp, md, qjd);
	init_dump_periods(rp, cp, md, qjd);
	init_obs_std(rp, cp, md, qjd);
	init_obs_lpn(rp, cp, md, qjd);

	init_start_state(rp, cp, md, qjd, 0);
}

void StdInitBehaviour::init_data(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd) const
{
	int num_trajectories = cp->qj_num_trajectories;

	init_splits(rp, cp, md, qjd);
	init_streams(rp, cp, md, qjd);
	copy_streams(rp, cp, md, qjd);
	leap_frog_all_streams(rp, cp, md, qjd);
	init_basic_data(rp, cp, md, qjd);
	init_dump_periods(rp, cp, md, qjd);
	init_obs_std(rp, cp, md, qjd);

	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		init_start_state(rp, cp, md, qjd, tr_id);
	}
}

void CorrDimInitBehaviour::init_data(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd) const
{
	int num_trajectories = cp->qj_num_trajectories;

	int cd_dump_deep = int(cp->params.find("cd_dump_deep")->second);

	init_splits_cd(rp, cp, md, qjd);
	init_streams(rp, cp, md, qjd);
	copy_streams(rp, cp, md, qjd);
	leap_frog_all_streams(rp, cp, md, qjd);
	init_basic_data(rp, cp, md, qjd);
	
	if(cd_dump_deep == 1)
	{
		init_dump_periods_cd_deep(rp, cp, md, qjd);
	}
	else
	{
		init_dump_periods(rp, cp, md, qjd);
	}
	
	init_obs_std(rp, cp, md, qjd);
	init_obs_cd(rp, cp, md, qjd);

	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		init_start_state(rp, cp, md, qjd, tr_id);
	}
}

void init_splits_cd(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd)
{
	int num_branches = md->num_ham_qj;
	int num_threads = rp->num_threads;
	
	int num_total = num_threads * num_branches;

	md->structure = init_split_structure_cd(rp, cp, md);
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

void init_splits(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd)
{
	int num_threads = rp->num_threads;

	md->structure = init_split_structure(rp, cp, md);
	md->splits = new Split[num_threads];
	for (int th_id = 0; th_id < num_threads; th_id++)
	{
		copy_struct_not_member(md->structure, &(md->splits)[th_id]);
	}
}

void init_streams(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd)
{
	int num_trajectories = cp->qj_num_trajectories;
	int seed = cp->qj_seed;
	int mns = cp->qj_mns;

	qjd->streams = new VSLStreamStatePtr[num_trajectories];
	vslNewStream(&(qjd->streams)[0], VSL_BRNG_MCG59, 777);
}

void leap_frog_single_stream(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int tr_id)
{
	int seed = cp->qj_seed;
	int mns = cp->qj_mns;
	vslLeapfrogStream((qjd->streams)[tr_id], seed + tr_id, mns);
}

void leap_frog_all_streams(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd)
{
	int num_trajectories = cp->qj_num_trajectories;

	int seed = cp->qj_seed;
	int mns = cp->qj_mns;

	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		vslLeapfrogStream((qjd->streams)[tr_id], seed + tr_id, mns);
	}
}

void copy_streams(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd)
{
	int num_trajectories = cp->qj_num_trajectories;
	int seed = cp->qj_seed;
	int mns = cp->qj_mns;

	for (int tr_id = 1; tr_id < num_trajectories; tr_id++)
	{
		vslCopyStream(&(qjd->streams)[tr_id], (qjd->streams)[0]);
	}
}

void init_streams_var(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd)
{
	int num_trajectories = cp->qj_num_trajectories;

	qjd->streams_var = new VSLStreamStatePtr[num_trajectories];
	vslNewStream(&(qjd->streams_var)[0], VSL_BRNG_MCG31, 777);
	for (int tr_id = 1; tr_id < num_trajectories; tr_id++)
	{
		vslCopyStream(&(qjd->streams_var)[tr_id], (qjd->streams_var)[0]);
	}
	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		vslLeapfrogStream((qjd->streams_var)[tr_id], tr_id, num_trajectories);
	}
}

void init_basic_data(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd)
{
	int sys_size = md->sys_size;
	int num_trajectories = cp->qj_num_trajectories;
	
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

void init_dump_periods_cd_deep(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd)
{
	qjd->period_id = 0;

	qjd->dump_type = int(cp->params.find("dump_type")->second);
	int num_dumps = int(cp->params.find("num_dumps")->second);
	int num_obs_periods = cp->qj_num_obs_periods;

	int num_branches = md->num_ham_qj;
	int num_sub_steps = num_branches * int(cp->params.find("cd_num_sub_steps")->second);
	int num_dumps_total = num_sub_steps * cp->qj_num_obs_periods + 1;

	qjd->num_dumps_total = num_dumps_total;
	qjd->dump_periods = new int[qjd->num_dumps_total];

	qjd->dump_periods[0] = 0;

	int dump_shift = md->T / double(num_sub_steps);
	for (int period_id = 0; period_id < num_obs_periods; period_id++)
	{
		for (int step_id = 0; step_id < num_sub_steps; step_id++)
		{
			int index = period_id * num_sub_steps + step_id + 1;
			qjd->dump_periods[index] = period_id;
		}
	}
}

void init_dump_periods(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd)
{
	qjd->period_id = 0;

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
	int sys_size = md->sys_size;
	int num_trajectories = cp->qj_num_trajectories;
	int num_dumps_total = qjd->num_dumps_total;

	double prm_E = double(cp->params.find("prm_E")->second);
	double * hamiltonian = md->hamiltonian;
	double * hamiltonian_drv = md->hamiltonian_drv;
	qjd->max_energy = 0.0;
	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		int index = st_id * sys_size + st_id;
		double ham_val = (hamiltonian[index] + prm_E * hamiltonian_drv[index]);
		if (abs(ham_val) > qjd->max_energy)
		{
			qjd->max_energy = abs(ham_val);
		}
	}

	qjd->delta_s = new double[num_trajectories];

	qjd->energy = new double[num_trajectories];
	qjd->lambda = new double[num_trajectories];
	qjd->lambda_now = new double[num_trajectories];
	qjd->mean_lpn = new double[num_trajectories];
	qjd->energy_lpn = new double[num_trajectories];

	qjd->energy_evo = new double[num_trajectories * num_dumps_total];
	qjd->lambda_evo = new double[num_trajectories * num_dumps_total];
	qjd->mean_lpn_evo = new double[num_trajectories * num_dumps_total];
	qjd->energy_lpn_evo = new double[num_trajectories * num_dumps_total];

	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		qjd->delta_s[tr_id] = 0.0;

		qjd->energy[tr_id] = 0.0;
		qjd->lambda[tr_id] = 0.0;
		qjd->lambda_now[tr_id] = 0.0;
		qjd->mean_lpn[tr_id] = 0.0;
		qjd->energy_lpn[tr_id] = 0.0;

		for (int dump_id = 0; dump_id < num_dumps_total; dump_id++)
		{
			qjd->energy_evo[tr_id * num_dumps_total + dump_id] = 0.0;
			qjd->lambda_evo[tr_id * num_dumps_total + dump_id] = 0.0;
			qjd->mean_lpn_evo[tr_id * num_dumps_total + dump_id] = 0.0;
			qjd->energy_lpn_evo[tr_id * num_dumps_total + dump_id] = 0.0;
		}	
	}
}

void init_obs_cd(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd)
{
	int num_branches = md->num_ham_qj;
	int num_sub_steps = num_branches * int(cp->params.find("cd_num_sub_steps")->second);
	int cd_dim = int(cp->params.find("cd_dim")->second);

	int sys_size = md->sys_size;
	int num_trajectories = cp->qj_num_trajectories;
	int num_periods = cp->qj_num_obs_periods;

	qjd->cd_shift_size = 1;
	qjd->cd_dim = cd_dim;
	qjd->cd_num_points = num_periods * num_sub_steps - cd_dim + 1;

	qjd->cd_i = new double[num_trajectories];
	qjd->cd_rec_data = new double**[num_trajectories];
	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		qjd->cd_i[tr_id] = 0.0;
		qjd->cd_rec_data[tr_id] = new double*[qjd->cd_num_points];

		for (int p_id = 0; p_id < qjd->cd_num_points; p_id++)
		{
			qjd->cd_rec_data[tr_id][p_id] = new double[qjd->cd_dim];

			for (int st_id = 0; st_id < qjd->cd_dim; st_id++)
			{
				qjd->cd_rec_data[tr_id][p_id][st_id] = 0.0;
			}
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

