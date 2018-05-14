#include "data.h"
#include "f_basis.h"

void init_main_data(RunParam &rp, ConfigParam &cp, MainData &md)
{
	md.size = cp.N;

	md.T = cp.t_0 + cp.t_1;

	md.step_t_0 = cp.t_0 / double(cp.num_steps_t_0);
	md.step_t_1 = cp.t_1 / double(cp.num_steps_t_1);

	cp.g = 0.25 / cp.prm_alpha;
}

void delete_main_data(MainData &md)
{
	cout << "Hello" << endl;
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

	calcODE_trans(model, rp, cp, md, pd);
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

void f_basis_prop_deep(RunParam &rp, ConfigParam &cp, MainData &md)
{
	string fn;

	init_main_data(rp, cp, md);

	Model * model;
	model = createModel(md.size - 1, cp);

	f_basis_init(model, rp, cp, md);

	PropData pd;
	init_prop_data_deep(rp, cp, md, pd);

	calcODE_trans(model, rp, cp, md, pd);
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

		int dump_shift = cp.num_periods_obser / cp.int_dn;
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
		double end_decade = log10(double(cp.num_periods_obser));
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

void init_prop_data_deep(RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd)
{
	pd.total_num_dumps = 1 + (cp.int_dn + cp.int_dn) * cp.num_periods_obser;
	pd.deep_dump_times = new double[pd.total_num_dumps];

	pd.deep_evals = new double[pd.total_num_dumps * md.size];
	pd.deep_avg_rho = new MKL_Complex16[md.size * md.size];

	int dupm_step_t_0 = cp.num_steps_t_0 / cp.int_dn;
	int dupm_step_t_1 = cp.num_steps_t_1 / cp.int_dn;

	double curr_time = 0.0;
	int curr_dump_id = 0;

	pd.deep_dump_times[curr_dump_id] = curr_time;
	curr_dump_id++;

	for (int per_id = 0; per_id < cp.num_periods_obser; per_id++)
	{
		for (int t0_id = 0; t0_id < cp.num_steps_t_0; t0_id++)
		{
			curr_time += md.step_t_0;

			if (t0_id % dupm_step_t_0 == 0)
			{
				pd.deep_dump_times[curr_dump_id] = curr_time;
				curr_dump_id++;
			}
		}

		curr_time = double(per_id) * md.T + cp.t_0;

		for (int t1_id = 0; t1_id < cp.num_steps_t_1; t1_id++)
		{
			curr_time += md.step_t_1;

			if (t1_id % dupm_step_t_1 == 0)
			{
				pd.deep_dump_times[curr_dump_id] = curr_time;
				curr_dump_id++;
			}
		}

		curr_time = double(per_id + 1) * md.T;
	}

	for (int dump_id = 0; dump_id < pd.total_num_dumps; dump_id++)
	{
		for (int st_id = 0; st_id < md.size; st_id++)
		{
			pd.deep_evals[dump_id * md.size + st_id] = 0.0;
		}
	}

	for (int st_id_1 = 0; st_id_1 < md.size; st_id_1++)
	{
		for (int st_id_2 = 0; st_id_2 < md.size; st_id_2++)
		{
			pd.deep_avg_rho[st_id_1 * md.size + st_id_2].real = 0.0;
			pd.deep_avg_rho[st_id_1 * md.size + st_id_2].imag = 0.0;
		}
	}

	string fn = rp.path + "periods" + file_name_suffix(cp, 4);
	cout << "Saving times to file:" << endl << fn << endl << endl;
	save_double_data(fn, pd.deep_dump_times, pd.total_num_dumps, 16, false);

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

	delete[] pd.k1;
	delete[] pd.k2;
	delete[] pd.k3;
	delete[] pd.k4;
	delete[] pd.val;
	delete[] pd.tmp;
}

void free_prop_data_deep(RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd)
{
	delete[] pd.deep_dump_times;

	delete[] pd.deep_evals;
	delete[] pd.deep_avg_rho;

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

	// Add here saving regular characteristics
}

void dump_prop_data_deep(RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd)
{
	cout << "Saving specs to files" << endl << endl;

	string fn = rp.path;

	fn = rp.path + "evals" + file_name_suffix(cp, 4);
	cout << "Saving evals to file:" << endl << fn << endl << endl;
	save_double_data(fn, pd.deep_evals, pd.total_num_dumps * md.size, 16, false);

	for (int st_id_1 = 0; st_id_1 < md.size; st_id_1++)
	{
		for (int st_id_2 = 0; st_id_2 < md.size; st_id_2++)
		{
			pd.deep_avg_rho[st_id_1 * md.size + st_id_2].real /= double(pd.total_num_dumps);
			pd.deep_avg_rho[st_id_1 * md.size + st_id_2].imag /= double(pd.total_num_dumps);
		}
	}
	
	fn = rp.path + "rho_avg" + file_name_suffix(cp, 4);
	cout << "Saving times to file:" << endl << fn << endl << endl;
	save_complex_data(fn, pd.deep_avg_rho, md.size * md.size, 16, false);
}
