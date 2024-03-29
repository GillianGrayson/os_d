#include "newdel_exp.h"
#include "qj_proc.h"

void LpnNewDelBehaviour::init_data(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const
{
	pb->init_prop_data(ad, cb);
	init_streams(ad);
	leap_frog_single_stream(ad, 0);
	init_streams_var(ad);
	init_basic_data(ad);
	init_dump_periods(ad);
	init_obs_std(ad);
	init_obs_lpn(ad);

	init_start_state(ad, 0);
}

void LpnNewDelBehaviour::free_data(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const
{
	pb->free_prop_data(ad, cb);
	free_streams(ad);
	free_streams_var(ad);
	free_basic_data(ad);
	free_dump_priods(ad);
	free_obs_std(ad);
	free_obs_lpn(ad);
}

void LpnMultNewDelBehaviour::init_data(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const
{
	int num_trajectories = ad->cp->num_trajectories;

	pb->init_prop_data(ad, cb);
	init_streams(ad);
	copy_half_streams(ad);
	for (int tr_id = 0; tr_id < num_trajectories / 2; tr_id++)
	{
		leap_frog_single_stream(ad, tr_id);
	}
	init_streams_var(ad);
	init_basic_data(ad);
	init_dump_periods(ad);
	init_obs_std(ad);
	init_obs_lpn(ad);

	for (int tr_id = 0; tr_id < num_trajectories / 2; tr_id++)
	{
		init_start_state(ad, tr_id);
	}
}

void LpnMultNewDelBehaviour::free_data(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const
{
	pb->free_prop_data(ad, cb);
	free_streams(ad);
	free_streams_var(ad);
	free_basic_data(ad);
	free_dump_priods(ad);
	free_obs_std(ad);
	free_obs_lpn(ad);
}

void LpnMultDeepNewDelBehaviour::init_data(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const
{
	int num_trajectories = ad->cp->num_trajectories;

	pb->init_prop_data_deep(ad, cb);
	init_streams(ad);
	copy_half_streams(ad);
	for (int tr_id = 0; tr_id < num_trajectories / 2; tr_id++)
	{
		leap_frog_single_stream(ad, tr_id);
	}
	init_streams_var(ad);
	init_basic_data(ad);
	init_dump_periods_deep(ad);
	init_obs_std(ad);
	init_obs_lpn(ad);

	for (int tr_id = 0; tr_id < num_trajectories / 2; tr_id++)
	{
		init_start_state(ad, tr_id);
	}
}

void LpnMultDeepNewDelBehaviour::free_data(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const
{
	pb->free_prop_data_deep(ad, cb);
	free_streams(ad);
	free_streams_var(ad);
	free_basic_data(ad);
	free_dump_priods(ad);
	free_obs_std(ad);
	free_obs_lpn(ad);
}


void StdNewDelBehaviour::init_data(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const
{
	ConfigParam * cp = ad->cp;

	int num_trajectories = cp->num_trajectories;

	pb->init_prop_data(ad, cb);
	cout << "init_prop_data" << endl;
	init_streams(ad);
	copy_streams(ad);
	leap_frog_all_streams(ad);
	init_basic_data(ad);
	init_dump_periods(ad);
	init_obs_std(ad);

	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		init_start_state(ad, tr_id);
	}
}

void StdNewDelBehaviour::free_data(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const
{
	pb->free_prop_data(ad, cb);
	free_streams(ad);
	free_basic_data(ad);
	free_dump_priods(ad);
	free_obs_std(ad);
}

void CorrDimNewDelBehaviour::init_data(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const
{
	ConfigParam * cp = ad->cp;

	int num_trajectories = cp->num_trajectories;

	pb->init_prop_data_deep(ad, cb);
	init_streams(ad);
	copy_streams(ad);
	leap_frog_all_streams(ad);
	init_basic_data(ad);
	
	init_dump_periods(ad);
	
	init_obs_std(ad);
	init_obs_cd(ad);

	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		init_start_state(ad, tr_id);
	}
}

void CorrDimNewDelBehaviour::free_data(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const
{
	pb->free_prop_data_deep(ad, cb);
	free_streams(ad);
	free_basic_data(ad);
	free_dump_priods(ad);
	free_obs_std(ad);
	free_obs_cd(ad);
}

void SigmaNewDelBehaviour::init_data(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const
{
	ConfigParam * cp = ad->cp;

	int num_trajectories = cp->num_trajectories;

	pb->init_prop_data_deep(ad, cb);
	init_streams(ad);
	copy_streams(ad);
	leap_frog_all_streams(ad);

	init_basic_data(ad);

	init_dump_periods(ad);

	init_obs_std(ad);

	init_start_state(ad, 0);
}

void SigmaNewDelBehaviour::free_data(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const
{
	pb->free_prop_data_deep(ad, cb);
	free_streams(ad);
	free_basic_data(ad);
	free_dump_priods(ad);
	free_obs_std(ad);
}

void StdDeepNewDelBehaviour::init_data(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const
{
	ConfigParam * cp = ad->cp;

	int num_trajectories = cp->num_trajectories;

	pb->init_prop_data_deep(ad, cb);
	init_streams(ad);
	copy_streams(ad);
	leap_frog_all_streams(ad);
	init_basic_data(ad);
	init_dump_periods_deep(ad);
	init_obs_std(ad);

	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		init_start_state(ad, tr_id);
	}
}

void StdDeepNewDelBehaviour::free_data(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const
{
	pb->free_prop_data_deep(ad, cb);
	free_streams(ad);
	free_basic_data(ad);
	free_dump_priods(ad);
	free_obs_std(ad);
}

void LpnDeepNewDelBehaviour::init_data(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const
{
	pb->init_prop_data_deep(ad, cb);
	init_streams(ad);
	leap_frog_single_stream(ad, 0);
	init_streams_var(ad);
	init_basic_data(ad);
	init_dump_periods_deep(ad);
	init_obs_std(ad);
	init_obs_lpn(ad);

	init_start_state(ad, 0);
}

void LpnDeepNewDelBehaviour::free_data(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const
{
	pb->free_prop_data_deep(ad, cb);
	free_streams(ad);
	free_streams_var(ad);
	free_basic_data(ad);
	free_dump_priods(ad);
	free_obs_std(ad);
	free_obs_lpn(ad);
}

void init_streams(AllData * ad)
{
	ConfigParam * cp = ad->cp;
	ExpData * ed = ad->ed;

	int num_trajectories = cp->num_trajectories;
	int seed = cp->seed;
	int mns = cp->mns;

	ed->streams = new VSLStreamStatePtr[num_trajectories];
	vslNewStream(&(ed->streams)[0], VSL_BRNG_MCG59, 777);
}

void leap_frog_single_stream(AllData * ad, int tr_id)
{
	ConfigParam * cp = ad->cp;
	ExpData * ed = ad->ed;

	int seed = cp->seed;
	int mns = cp->mns;
	vslLeapfrogStream((ed->streams)[tr_id], seed + tr_id, mns);
}

void leap_frog_all_streams(AllData * ad)
{
	ConfigParam * cp = ad->cp;
	ExpData * ed = ad->ed;

	int num_trajectories = cp->num_trajectories;

	int seed = cp->seed;
	int mns = cp->mns;

	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		vslLeapfrogStream((ed->streams)[tr_id], seed + tr_id, mns);
	}
}

void copy_streams(AllData * ad)
{
	ConfigParam * cp = ad->cp;
	ExpData * ed = ad->ed;

	int num_trajectories = cp->num_trajectories;
	int seed = cp->seed;
	int mns = cp->mns;

	for (int tr_id = 1; tr_id < num_trajectories; tr_id++)
	{
		vslCopyStream(&(ed->streams)[tr_id], (ed->streams)[0]);
	}
}

void copy_half_streams(AllData * ad)
{
	ConfigParam * cp = ad->cp;
	ExpData * ed = ad->ed;

	int num_trajectories = cp->num_trajectories;
	int seed = cp->seed;
	int mns = cp->mns;

	for (int tr_id = 1; tr_id < num_trajectories / 2; tr_id++)
	{
		vslCopyStream(&(ed->streams)[tr_id], (ed->streams)[0]);
	}
}

void init_streams_var(AllData * ad)
{
	ConfigParam * cp = ad->cp;
	ExpData * ed = ad->ed;

	int num_trajectories = cp->num_trajectories;

	ed->streams_var = new VSLStreamStatePtr[num_trajectories];
	vslNewStream(&(ed->streams_var)[0], VSL_BRNG_MCG31, 777);
	for (int tr_id = 1; tr_id < num_trajectories; tr_id++)
	{
		vslCopyStream(&(ed->streams_var)[tr_id], (ed->streams_var)[0]);
	}
	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		vslLeapfrogStream((ed->streams_var)[tr_id], tr_id, num_trajectories);
	}
}

void init_basic_data(AllData * ad)
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;
	int num_trajectories = cp->num_trajectories;
	
	ed->phi_all = new MKL_Complex16[num_trajectories * sys_size];
	ed->phi_all_aux = new MKL_Complex16[num_trajectories * sys_size];
	ed->abs_diag_rho_all = new double[num_trajectories * sys_size];
	ed->times_all = new double[num_trajectories];
	ed->etas_all = new double[num_trajectories];
	ed->jumps_counts = new int[num_trajectories];

	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		for (int st_id = 0; st_id < sys_size; st_id++)
		{
			ed->phi_all[tr_id * sys_size + st_id].real = 0.0;
			ed->phi_all[tr_id * sys_size + st_id].imag = 0.0;

			ed->phi_all_aux[tr_id * sys_size + st_id].real = 0.0;
			ed->phi_all_aux[tr_id * sys_size + st_id].imag = 0.0;

			ed->abs_diag_rho_all[tr_id * sys_size + st_id] = 0.0;
			ed->abs_diag_rho_all[tr_id * sys_size + st_id] = 0.0;
		}

		ed->times_all[tr_id] = 0.0;
		ed->etas_all[tr_id] = 0.0;
		ed->jumps_counts[tr_id] = 0.0;
	}

	int jump = int(cp->params.find("jump")->second);
	if (jump > 0)
	{
		ed->jump_times = new vector<double>[num_trajectories];
		ed->diss_types = new vector<int>[num_trajectories];
		ed->jump_norms = new vector<double>[num_trajectories];
		ed->jump_etas = new vector<double>[num_trajectories];
	}
}

void init_dump_periods_deep(AllData * ad)
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	ed->curr_time = 0.0;

	ed->dump_type = int(cp->params.find("dump_type")->second);
	int dump_num = int(cp->params.find("dump_num")->second);
	int num_obs_periods = cp->num_obs_periods;

	int num_branches = md->num_ham_qj;
	int num_sub_steps = num_branches * int(cp->params.find("deep_num_steps")->second);
	int dump_num_total = num_sub_steps * cp->num_obs_periods + 1;

	ed->dump_num_total = dump_num_total;
	ed->dump_periods = new int[ed->dump_num_total];

	ed->dump_periods[0] = 0;

	for (int period_id = 0; period_id < num_obs_periods; period_id++)
	{
		for (int step_id = 0; step_id < num_sub_steps; step_id++)
		{
			int index = period_id * num_sub_steps + step_id + 1;
			ed->dump_periods[index] = period_id;
		}
	}
}

void init_dump_periods(AllData * ad)
{
	ConfigParam * cp = ad->cp;
	ExpData * ed = ad->ed;

	ed->curr_time = 0.0;

	ed->dump_type = int(cp->params.find("dump_type")->second);
	int dump_num = int(cp->params.find("dump_num")->second);
	int num_obs_periods = cp->num_obs_periods;

	if (ed->dump_type == 0)
	{
		ed->dump_num_total = dump_num + 1;
		ed->dump_periods = new int[ed->dump_num_total];

		ed->dump_periods[0] = 0;

		int dump_shift = num_obs_periods / dump_num;
		for (int dump_id = 0; dump_id < dump_num; dump_id++)
		{
			ed->dump_periods[dump_id + 1] = (dump_id + 1) * dump_shift;
		}
	}
	else if (ed->dump_type == 1)
	{
		ed->dump_num_total = dump_num + 2;
		ed->dump_periods = new int[ed->dump_num_total];

		double begin_decade = log10(1.0);
		double end_decade = log10(double(num_obs_periods));
		double num_decades = end_decade - begin_decade;
		double num_decades_dump = double(dump_num) / num_decades;
		for (int dump_id = 0; dump_id < dump_num + 1; dump_id++)
		{
			int curr_val = int(pow(10.0, begin_decade) * pow(10.0, (1.0 / num_decades_dump) * double(dump_id)));

			if (curr_val > ed->dump_periods[dump_id])
			{
				ed->dump_periods[dump_id + 1] = curr_val;
			}
			else
			{
				ed->dump_periods[dump_id + 1] = ed->dump_periods[dump_id] + 1;
			}
		}
	}
}

void init_obs_std(AllData * ad)
{
	ConfigParam * cp = ad->cp;
	ExpData * ed = ad->ed;

	int num_trajectories = cp->num_trajectories;
	int dump_num_total = ed->dump_num_total;

	int num_random_obs = int(cp->params.find("num_random_obs")->second);

	ed->mean_start = new double[num_trajectories];

	ed->norm			= new double[num_trajectories];
	ed->mean			= new double[num_trajectories];
	ed->dispersion		= new double[num_trajectories];
	ed->m2				= new double[num_trajectories];
	ed->energy			= new double[num_trajectories];
	ed->spec			= new MKL_Complex16[num_trajectories];
	ed->spec_2 = new MKL_Complex16[num_trajectories];
	ed->spec_3 = new MKL_Complex16[num_trajectories];
	
	ed->norm_evo		= new double[num_trajectories * dump_num_total];
	ed->mean_evo		= new double[num_trajectories * dump_num_total];
	ed->dispersion_evo	= new double[num_trajectories * dump_num_total];
	ed->m2_evo			= new double[num_trajectories * dump_num_total];
	ed->energy_evo		= new double[num_trajectories * dump_num_total];
	ed->spec_evo		= new MKL_Complex16[num_trajectories * dump_num_total];
	ed->spec_2_evo = new MKL_Complex16[num_trajectories * dump_num_total];
	ed->spec_3_evo = new MKL_Complex16[num_trajectories * dump_num_total];
	
	if (num_random_obs > 0)
	{
		ed->random_obs = std::vector<std::vector<std::complex<double>>>(num_trajectories, std::vector<std::complex<double>>(num_random_obs));
		ed->random_obs_evo = std::vector<std::vector<std::vector<std::complex<double>>>>(num_trajectories, std::vector<std::vector<std::complex<double>>>(num_random_obs));
	}

	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		ed->mean_start[tr_id] = 0.0;

		ed->norm[tr_id]			= 0.0;
		ed->mean[tr_id]			= 0.0;
		ed->dispersion[tr_id]	= 0.0;
		ed->m2[tr_id]			= 0.0;
		ed->energy[tr_id]		= 0.0;
		ed->spec[tr_id].real	= 0.0;
		ed->spec[tr_id].imag	= 0.0;
		ed->spec_2[tr_id].real = 0.0;
		ed->spec_2[tr_id].imag = 0.0;
		ed->spec_3[tr_id].real = 0.0;
		ed->spec_3[tr_id].imag = 0.0;

		for (int dump_id = 0; dump_id < dump_num_total; dump_id++)
		{
			ed->norm_evo[tr_id * dump_num_total + dump_id]			= 0.0;
			ed->mean_evo[tr_id * dump_num_total + dump_id]			= 0.0;
			ed->dispersion_evo[tr_id * dump_num_total + dump_id]	= 0.0;
			ed->m2_evo[tr_id * dump_num_total + dump_id]			= 0.0;
			ed->energy_evo[tr_id * dump_num_total + dump_id]		= 0.0;
			ed->spec_evo[tr_id * dump_num_total + dump_id].real		= 0.0;
			ed->spec_evo[tr_id * dump_num_total + dump_id].imag		= 0.0;
			ed->spec_2_evo[tr_id * dump_num_total + dump_id].real = 0.0;
			ed->spec_2_evo[tr_id * dump_num_total + dump_id].imag = 0.0;
			ed->spec_3_evo[tr_id * dump_num_total + dump_id].real = 0.0;
			ed->spec_3_evo[tr_id * dump_num_total + dump_id].imag = 0.0;
		}
	}

	ed->is_obs = 0;
}

void init_obs_lpn(AllData * ad)
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int num_random_obs = int(cp->params.find("num_random_obs")->second);

	int sys_size = md->sys_size;
	int num_trajectories = cp->num_trajectories;
	int dump_num_total = ed->dump_num_total;

	MKL_Complex16* hamiltonian = md->hamiltonian;
	ed->max_energy = 0.0;
	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		int index = st_id * sys_size + st_id;
		double ham_val = hamiltonian[index].real;
		if (abs(ham_val) > ed->max_energy)
		{
			ed->max_energy = abs(ham_val);
		}
	}

	ed->delta_s = new double[num_trajectories];
	ed->lambda_now = new double[num_trajectories];
	ed->lambda = new double[num_trajectories];

	ed->num_renorms = new int[num_trajectories];

	int save_lambdas = int(ad->cp->params.find("save_lambdas")->second);
	if (save_lambdas > 0)
	{
		ed->lambdas = new vector<double>[num_trajectories];
		ed->deltas_s = new vector<double>[num_trajectories];
		ed->deltas_f = new vector<double>[num_trajectories];
	}

	ed->mean_lpn = new double[num_trajectories];
	ed->energy_lpn = new double[num_trajectories];
	ed->spec_lpn = new MKL_Complex16[num_trajectories];
	ed->spec_2_lpn = new MKL_Complex16[num_trajectories];
	ed->spec_3_lpn = new MKL_Complex16[num_trajectories];

	ed->lambda_evo = new double[num_trajectories * dump_num_total];

	ed->mean_lpn_evo = new double[num_trajectories * dump_num_total];
	ed->energy_lpn_evo = new double[num_trajectories * dump_num_total];
	ed->spec_lpn_evo = new MKL_Complex16[num_trajectories * dump_num_total];
	ed->spec_2_lpn_evo = new MKL_Complex16[num_trajectories * dump_num_total];
	ed->spec_3_lpn_evo = new MKL_Complex16[num_trajectories * dump_num_total];

	if (num_random_obs > 0)
	{
		ed->random_obs_lpn = std::vector<std::vector<std::complex<double>>>(num_trajectories, std::vector<std::complex<double>>(num_random_obs));
		ed->random_obs_lpn_evo = std::vector<std::vector<std::vector<std::complex<double>>>>(num_trajectories, std::vector<std::vector<std::complex<double>>>(num_random_obs));
	}

	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		ed->delta_s[tr_id] = 0.0;
		ed->lambda_now[tr_id] = 0.0;
		ed->lambda[tr_id] = 0.0;

		ed->num_renorms[tr_id] = 0;

		ed->mean_lpn[tr_id] = 0.0;
		ed->energy_lpn[tr_id] = 0.0;
		ed->spec_lpn[tr_id].real = 0.0;
		ed->spec_lpn[tr_id].imag = 0.0;
		ed->spec_2_lpn[tr_id].real = 0.0;
		ed->spec_2_lpn[tr_id].imag = 0.0;
		ed->spec_3_lpn[tr_id].real = 0.0;
		ed->spec_3_lpn[tr_id].imag = 0.0;

		for (int dump_id = 0; dump_id < dump_num_total; dump_id++)
		{
			ed->lambda_evo[tr_id * dump_num_total + dump_id] = 0.0;

			ed->mean_lpn_evo[tr_id * dump_num_total + dump_id] = 0.0;
			ed->energy_lpn_evo[tr_id * dump_num_total + dump_id] = 0.0;
			ed->spec_lpn_evo[tr_id * dump_num_total + dump_id].real = 0.0;
			ed->spec_lpn_evo[tr_id * dump_num_total + dump_id].imag = 0.0;
			ed->spec_2_lpn_evo[tr_id * dump_num_total + dump_id].real = 0.0;
			ed->spec_2_lpn_evo[tr_id * dump_num_total + dump_id].imag = 0.0;
			ed->spec_3_lpn_evo[tr_id * dump_num_total + dump_id].real = 0.0;
			ed->spec_3_lpn_evo[tr_id * dump_num_total + dump_id].imag = 0.0;
		}	
	}
}

void init_obs_cd(AllData * ad)
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int num_branches = md->num_ham_qj;
	int num_sub_steps = num_branches * int(cp->params.find("deep_num_steps")->second);
	int cd_dim = int(cp->params.find("cd_dim")->second);

	int sys_size = md->sys_size;
	int num_trajectories = cp->num_trajectories;
	int num_periods = cp->num_obs_periods;

	ed->cd_shift_size = 1;
	ed->cd_dim = cd_dim;
	ed->cd_num_points = num_periods * num_sub_steps - cd_dim + 1;

	ed->cd_i = new double[num_trajectories];
	ed->cd_rec_data = new double**[num_trajectories];
	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		ed->cd_i[tr_id] = 0.0;
		ed->cd_rec_data[tr_id] = new double*[ed->cd_num_points];

		for (int p_id = 0; p_id < ed->cd_num_points; p_id++)
		{
			ed->cd_rec_data[tr_id][p_id] = new double[ed->cd_dim];

			for (int st_id = 0; st_id < ed->cd_dim; st_id++)
			{
				ed->cd_rec_data[tr_id][p_id][st_id] = 0.0;
			}
		}
	}
}

void init_start_state(AllData * ad, int tr_id)
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;
	int start_type = int(cp->params.find("start_type")->second);
	int start_state = int(cp->params.find("start_state")->second);

	MKL_Complex16 * phi = &(ed->phi_all[tr_id * sys_size]);

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

void free_streams(AllData * ad)
{
	ExpData * ed = ad->ed;

	VSLStreamStatePtr * streams = ed->streams;
	delete[] streams;
}

void free_streams_var(AllData * ad)
{
	ExpData * ed = ad->ed;

	VSLStreamStatePtr * streams_var = ed->streams_var;
	delete[] streams_var;
}

void free_basic_data(AllData * ad)
{
	ExpData * ed = ad->ed;
	ConfigParam * cp = ad->cp;

	int num_trajectories = cp->num_trajectories;

	delete[] ed->phi_all;
	delete[] ed->phi_all_aux;
	delete[] ed->abs_diag_rho_all;
	delete[] ed->times_all;
	delete[] ed->etas_all;
	delete[] ed->jumps_counts;

	int jump = int(cp->params.find("jump")->second);
	if (jump > 0)
	{
		for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
		{
			ed->jump_times[tr_id].clear();
			ed->diss_types[tr_id].clear();
			ed->jump_norms[tr_id].clear();
			ed->jump_etas[tr_id].clear();
		}
		delete[] ed->jump_times;
		delete[] ed->diss_types;
		delete[] ed->jump_norms;
		delete[] ed->jump_etas;
	}
}

void free_dump_priods(AllData * ad)
{
	ExpData * ed = ad->ed;

	delete[] ed->dump_periods;
}

void free_obs_std(AllData * ad)
{
	ExpData * ed = ad->ed;

	delete[] ed->norm;
	delete[] ed->mean_start;
	delete[] ed->mean;
	delete[] ed->dispersion;
	delete[] ed->m2;
	delete[] ed->energy;
	delete[] ed->spec;
	delete[] ed->spec_2;
	delete[] ed->spec_3;

	delete[] ed->norm_evo;
	delete[] ed->mean_evo;
	delete[] ed->dispersion_evo;
	delete[] ed->m2_evo;
	delete[] ed->energy_evo;
	delete[] ed->spec_evo;
	delete[] ed->spec_2_evo;
	delete[] ed->spec_3_evo;
}

void free_obs_lpn(AllData * ad)
{
	ExpData * ed = ad->ed;

	delete[] ed->delta_s;

	delete[] ed->lambda;
	delete[] ed->lambda_now;
	delete[] ed->mean_lpn;
	delete[] ed->energy_lpn;
	delete[] ed->spec_lpn;
	delete[] ed->spec_2_lpn;
	delete[] ed->spec_3_lpn;

	delete[] ed->num_renorms;
	int num_trajectories = ad->cp->num_trajectories;
	int save_lambdas = int(ad->cp->params.find("save_lambdas")->second);
	if (save_lambdas > 0)
	{
		for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
		{
			ed->lambdas[tr_id].clear();
			ed->deltas_s[tr_id].clear();
			ed->deltas_f[tr_id].clear();
		}
		delete[] ed->lambdas;
		delete[] ed->deltas_s;
		delete[] ed->deltas_f;
	}

	delete[] ed->lambda_evo;
	delete[] ed->mean_lpn_evo;
	delete[] ed->energy_lpn_evo;
	delete[] ed->spec_lpn_evo;
	delete[] ed->spec_2_lpn_evo;
	delete[] ed->spec_3_lpn_evo;
}

void free_obs_cd(AllData * ad)
{
	ConfigParam * cp = ad->cp;
	ExpData * ed = ad->ed;

	int num_trajectories = cp->num_trajectories;

	delete[] ed->cd_i;

	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		for (int p_id = 0; p_id < ed->cd_num_points; p_id++)
		{
			delete[] ed->cd_rec_data[tr_id][p_id];
		}
		delete[] ed->cd_rec_data[tr_id];
	}
	delete[] ed->cd_rec_data;
}
