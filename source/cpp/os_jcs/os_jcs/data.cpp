#include "data.h"
#include "f_basis.h"

void init_main_data(RunParam &rp, ConfigParam &cp, MainData &md)
{
	md.size = pow(2, cp.N);

	md.T = cp.t_0 + cp.t_1;

	md.step_t_0 = cp.t_0 / double(cp.num_steps_t_0);
	md.step_t_1 = cp.t_1 / double(cp.num_steps_t_1);

	md.disorder_J = new double[cp.N - 1];
	md.disorder_h_z = new double[cp.N];
	md.disorder_h_x = new double[cp.N];

	cout << "Generating disorders ... " << endl;

	VSLStreamStatePtr stream;
	vslNewStream(&stream, VSL_BRNG_MCG31, 77778888);
	vslLeapfrogStream(stream, cp.seed, cp.max_num_seeds);

	double J_begin = 0.5 * cp.J;
	double J_end = 1.5 * cp.J;

	vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, cp.N - 1, md.disorder_J, J_begin, J_end);
	vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, cp.N, md.disorder_h_z, 0.0, cp.h_z);
	vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, cp.N, md.disorder_h_x, 0.0, 1.0);

	for (int spin_id = 0; spin_id < cp.N; spin_id++)
	{
		md.disorder_h_x[spin_id] *= cp.h_x;
	}

	string disorder_fn = rp.path + "J_disorder" + file_name_suffix(cp, 4);
	save_double_data(disorder_fn, md.disorder_J, cp.N - 1, 16, false);

	disorder_fn = rp.path + "h_z_disorder" + file_name_suffix(cp, 4);
	save_double_data(disorder_fn, md.disorder_h_z, cp.N, 16, false);

	disorder_fn = rp.path + "h_x_disorder" + file_name_suffix(cp, 4);
	save_double_data(disorder_fn, md.disorder_h_x, cp.N, 16, false);
}

void delete_main_data(MainData &md)
{
	delete[] md.disorder_J;
	delete[] md.disorder_h_z;
	delete[] md.disorder_h_x;
}

void f_basis_prop_std(RunParam &rp, ConfigParam &cp, MainData &md)
{
	string fn;

	init_main_data(rp, cp, md);

	Model * model;
	model = createModel(md.size - 1, cp);

	f_basis_init(model, rp, cp, md);

	PropData pd;
	init_prop_data_std(rp, cp, md, pd);

	calcODE_std(model, rp, cp, md, pd);

	dump_prop_data_std(rp, cp, md, pd);

	free_prop_data_std(rp, cp, md, pd);

	if (rp.issmtx == 1)
	{
		fn = rp.path + "rho" + file_name_suffix(cp, 4);
		cout << "Saving rho to file:" << endl << fn << endl << endl;
		save_sparse_complex_mtx(fn, model->Rho, 16, false);
	}

	freeModel(model);
	delete_main_data(md);
}

void f_basis_prop_rate(RunParam &rp, ConfigParam &cp, MainData &md)
{
	string fn;

	init_main_data(rp, cp, md);

	Model * model;
	model = createModel(md.size - 1, cp);

	f_basis_init(model, rp, cp, md);

	PropData pd;
	init_prop_data_rate(rp, cp, md, pd);

	calcODE_rate(model, rp, cp, md, pd);

	dump_prop_data_rate(rp, cp, md, pd);

	free_prop_data_rate(rp, cp, md, pd);

	if (rp.issmtx == 1)
	{
		fn = rp.path + "rho" + file_name_suffix(cp, 4);
		cout << "Saving rho to file:" << endl << fn << endl << endl;
		save_sparse_complex_mtx(fn, model->Rho, 16, false);
	}

	freeModel(model);
	delete_main_data(md);
}

void f_basis_prop_deep(RunParam &rp, ConfigParam &cp, MainData &md)
{
	string fn;

	init_main_data(rp, cp, md);

	Model * model;
	model = createModel(md.size - 1, cp);

	f_basis_init(model, rp, cp, md);

	PropData pd;
	init_prop_data_deep(rp, cp, md, pd);

	calcODE_deep(model, rp, cp, md, pd);

	dump_prop_data_deep(rp, cp, md, pd);

	free_prop_data_deep(rp, cp, md, pd);

	if (rp.issmtx == 1)
	{
		fn = rp.path + "rho" + file_name_suffix(cp, 4);
		cout << "Saving rho to file:" << endl << fn << endl << endl;
		save_sparse_complex_mtx(fn, model->Rho, 16, false);
	}

	freeModel(model);
	delete_main_data(md);
}

void init_prop_data_std(RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd)
{
	if (cp.int_dt == 0)
	{
		pd.total_num_dumps = cp.int_dn + 1;
		pd.dump_periods = new int[cp.int_dn + 1];

		pd.dump_periods[0] = 0;

		int dump_shift = cp.num_periods / cp.int_dn;
		for (int dump_id = 0; dump_id < cp.int_dn; dump_id++)
		{
			pd.dump_periods[dump_id + 1] =  (dump_id + 1) * dump_shift;
		}
	}
	else if (cp.int_dt == 1)
	{
		pd.total_num_dumps = cp.int_dn + 2;
		pd.dump_periods = new int[cp.int_dn + 2];

		pd.dump_periods[0] = 0;

		double begin_decade = log10(1.0);
		double end_decade = log10(double(cp.num_periods));
		double num_decades = end_decade - begin_decade;
		double num_decades_dump = double(cp.int_dn) / num_decades;
		for (int dump_id = 0; dump_id < cp.int_dn + 1; dump_id++)
		{
			int curr_val = int(pow(10.0, begin_decade) * pow(10.0, (1.0 / num_decades_dump) * double(dump_id)));

			if (curr_val > pd.dump_periods[dump_id])
			{
				pd.dump_periods[dump_id + 1] = curr_val;
			}
			else
			{
				pd.dump_periods[dump_id + 1] = pd.dump_periods[dump_id] + 1;
			}
		}
	}
	else
	{
		stringstream msg;
		msg << "Wrong int_dt value: " << cp.int_dt << endl;
		Error(msg.str());
	}

	pd.sign_sigma_z_start = new int[cp.N];

	pd.magnetization = new double*[cp.N];
	pd.sigma_x = new double*[cp.N];
	pd.sigma_y = new double*[cp.N];
	pd.sigma_z = new double*[cp.N];
	
	for (int spin_id = 0; spin_id < cp.N; spin_id++)
	{
		pd.sign_sigma_z_start[spin_id] = 0;

		pd.magnetization[spin_id] = new double[pd.total_num_dumps];
		pd.sigma_x[spin_id] = new double[pd.total_num_dumps];
		pd.sigma_y[spin_id] = new double[pd.total_num_dumps];
		pd.sigma_z[spin_id] = new double[pd.total_num_dumps];

		for (int dump_id = 0; dump_id < pd.total_num_dumps; dump_id++)
		{
			pd.magnetization[spin_id][dump_id] = 0.0;
			pd.sigma_x[spin_id][dump_id] = 0.0;
			pd.sigma_y[spin_id][dump_id] = 0.0;
			pd.sigma_z[spin_id][dump_id] = 0.0;
		}
	}

	string fn = rp.path + "periods" + file_name_suffix(cp, 4);
	cout << "Saving times to file:" << endl << fn << endl << endl;
	save_int_data(fn, pd.dump_periods, pd.total_num_dumps, false);

	int N_mat = md.size * md.size - 1;

	pd.k1 = new double[N_mat];
	pd.k2 = new double[N_mat];
	pd.k3 = new double[N_mat];
	pd.k4 = new double[N_mat];
	pd.val = new double[N_mat];
	pd.tmp = new double[N_mat];
}

void init_prop_data_rate(RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd)
{
	init_prop_data_std(rp, cp, md, pd);

	pd.rate_result = new int[cp.N];
	pd.rate_flag = new int[cp.N];
	pd.rate_ampl_start = new double[cp.N];
	pd.rate_ampl_curr = new double[cp.N];
	pd.rate_max_curr = new double[cp.N];
	pd.rate_min_curr = new double[cp.N];

	for (int spin_id = 0; spin_id < cp.N; spin_id++)
	{
		pd.rate_result[spin_id] = cp.num_periods;
		pd.rate_flag[spin_id] = 0;
		pd.rate_ampl_start[spin_id] = 0.0;
		pd.rate_ampl_curr[spin_id] = 0.0;
		pd.rate_max_curr[spin_id] = -1.0e16;
		pd.rate_min_curr[spin_id] = 1.0e16;
	}
}

void init_prop_data_deep(RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd)
{
	pd.total_num_dumps = 1 + (cp.int_dn + cp.int_dn) * cp.num_periods;
	pd.dump_times = new double[pd.total_num_dumps];

	int dupm_step_t_0 = cp.num_steps_t_0 / cp.int_dn;
	int dupm_step_t_1 = cp.num_steps_t_1 / cp.int_dn;

	double curr_time = 0.0;
	int curr_dump_id = 0;

	pd.dump_times[curr_dump_id] = curr_time;
	curr_dump_id++;

	for (int per_id = 0; per_id < cp.num_periods; per_id++)
	{
		for (int t0_id = 0; t0_id < cp.num_steps_t_0; t0_id++)
		{
			curr_time += md.step_t_0;

			if (t0_id % dupm_step_t_0 == 0)
			{
				pd.dump_times[curr_dump_id] = curr_time;
				curr_dump_id++;
			}
		}

		curr_time = double(per_id) * md.T + cp.t_0;

		for (int t1_id = 0; t1_id < cp.num_steps_t_1; t1_id++)
		{
			curr_time += md.step_t_1;

			if (t1_id % dupm_step_t_1 == 0)
			{
				pd.dump_times[curr_dump_id] = curr_time;
				curr_dump_id++;
			}
		}

		curr_time = double(per_id + 1) * md.T;
	}

	pd.sigma_x = new double*[cp.N];
	pd.sigma_y = new double*[cp.N];
	pd.sigma_z = new double*[cp.N];

	for (int spin_id = 0; spin_id < cp.N; spin_id++)
	{
		pd.sigma_x[spin_id] = new double[pd.total_num_dumps];
		pd.sigma_y[spin_id] = new double[pd.total_num_dumps];
		pd.sigma_z[spin_id] = new double[pd.total_num_dumps];

		for (int dump_id = 0; dump_id < pd.total_num_dumps; dump_id++)
		{
			pd.sigma_x[spin_id][dump_id] = 0.0;
			pd.sigma_y[spin_id][dump_id] = 0.0;
			pd.sigma_z[spin_id][dump_id] = 0.0;
		}
	}

	string fn = rp.path + "periods" + file_name_suffix(cp, 4);
	cout << "Saving times to file:" << endl << fn << endl << endl;
	save_double_data(fn, pd.dump_times, pd.total_num_dumps, 16, false);

	int N_mat = md.size * md.size - 1;

	pd.k1 = new double[N_mat];
	pd.k2 = new double[N_mat];
	pd.k3 = new double[N_mat];
	pd.k4 = new double[N_mat];
	pd.val = new double[N_mat];
	pd.tmp = new double[N_mat];
}

void free_prop_data_std(RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd)
{
	delete[] pd.dump_periods;

	delete[] pd.sign_sigma_z_start;

	for (int spin_id = 0; spin_id < cp.N; spin_id++)
	{
		delete[] pd.magnetization[spin_id];
		delete[] pd.sigma_x[spin_id];
		delete[] pd.sigma_y[spin_id];
		delete[] pd.sigma_z[spin_id];
	}
	delete[] pd.magnetization;
	delete[] pd.sigma_x;
	delete[] pd.sigma_y;
	delete[] pd.sigma_z;

	delete[] pd.k1;
	delete[] pd.k2;
	delete[] pd.k3;
	delete[] pd.k4;
	delete[] pd.val;
	delete[] pd.tmp;
}

void free_prop_data_rate(RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd)
{
	free_prop_data_std(rp, cp, md, pd);

	delete[] pd.rate_result;
	delete[] pd.rate_flag;
	delete[] pd.rate_ampl_start;
	delete[] pd.rate_ampl_curr;
	delete[] pd.rate_max_curr;
	delete[] pd.rate_min_curr;
}

void free_prop_data_deep(RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd)
{
	delete[] pd.dump_times;

	for (int spin_id = 0; spin_id < cp.N; spin_id++)
	{
		delete[] pd.sigma_x[spin_id];
		delete[] pd.sigma_y[spin_id];
		delete[] pd.sigma_z[spin_id];
	}

	delete[] pd.sigma_x;
	delete[] pd.sigma_y;
	delete[] pd.sigma_z;

	delete[] pd.k1;
	delete[] pd.k2;
	delete[] pd.k3;
	delete[] pd.k4;
	delete[] pd.val;
	delete[] pd.tmp;
}

void dump_prop_data_std(RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd)
{
	cout << "Saving specs to files" << endl << endl;

	string fn = rp.path + "magnetization" + file_name_suffix(cp, 4);
	save_2d_double_data(fn, pd.magnetization, cp.N, pd.total_num_dumps, 16, false);

	fn = rp.path + "sigma_x" + file_name_suffix(cp, 4);
	save_2d_double_data(fn, pd.sigma_x, cp.N, pd.total_num_dumps, 16, false);

	fn = rp.path + "sigma_y" + file_name_suffix(cp, 4);
	save_2d_double_data(fn, pd.sigma_y, cp.N, pd.total_num_dumps, 16, false);

	fn = rp.path + "sigma_z" + file_name_suffix(cp, 4);
	save_2d_double_data(fn, pd.sigma_z, cp.N, pd.total_num_dumps, 16, false);
}

void dump_prop_data_rate(RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd)
{
	dump_prop_data_std(rp, cp, md, pd);
	cout << "Saving specs to files" << endl << endl;

	string fn = rp.path + "rate_periods" + file_name_suffix(cp, 4);
	save_int_data(fn, pd.rate_result, cp.N, false);
}

void dump_prop_data_deep(RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd)
{
	cout << "Saving specs to files" << endl << endl;

	string fn = rp.path;

	fn = rp.path + "sigma_x" + file_name_suffix(cp, 4);
	save_2d_double_data(fn, pd.sigma_x, cp.N, pd.total_num_dumps, 16, false);

	fn = rp.path + "sigma_y" + file_name_suffix(cp, 4);
	save_2d_double_data(fn, pd.sigma_y, cp.N, pd.total_num_dumps, 16, false);

	fn = rp.path + "sigma_z" + file_name_suffix(cp, 4);
	save_2d_double_data(fn, pd.sigma_z, cp.N, pd.total_num_dumps, 16, false);
}
