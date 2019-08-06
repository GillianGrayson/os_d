#include "data.h"
#include "Model.h"
#include "calcODE.h"
#include "f_basis_init.h"

void init_main_data(RunParam &rp, ConfigParam &cp, MainData &md)
{
	md.size = cp.N + 1;
	md.T = 2 * 3.14159265358979323846;
	md.step = md.T / double(cp.num_steps);
	cp.diss_gamma /= cp.N;
}

void delete_main_data(MainData &md)
{
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
		fn = "rho" + file_name_suffix(cp, 4);
		cout << "Saving rho to file:" << endl << fn << endl << endl;
		save_sparse_complex_mtx(fn, model->Rho, 16, false);
	}

	freeModel(model);
	delete_main_data(md);
}

void f_basis_prop_floquet(RunParam &rp, ConfigParam &cp, MainData &md)
{
	string fn;

	init_main_data(rp, cp, md);

	Model * model;
	model = createModel(md.size - 1, cp);

	f_basis_init(model, rp, cp, md);

	PropData pd;
	init_prop_data_floquet(rp, cp, md, pd);

	calcODE_floquet(model, rp, cp, md, pd);

	dump_prop_data_floquet(rp, cp, md, pd);

	free_prop_data_floquet(rp, cp, md, pd);

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
			pd.dump_periods[dump_id + 1] = (dump_id + 1) * dump_shift;
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

	string fn = "periods" + file_name_suffix(cp, 4);
	cout << "Saving times to file:" << endl << fn << endl << endl;
	save_int_data(fn, pd.dump_periods, pd.total_num_dumps, false);

	int N_mat = md.size * md.size - 1;

	pd.k1 = new double[N_mat];
	pd.k2 = new double[N_mat];
	pd.k3 = new double[N_mat];
	pd.k4 = new double[N_mat];
	pd.val = new double[N_mat];
	pd.tmp = new double[N_mat];
	pd.tmp_drv = new double[N_mat];
}

void init_prop_data_floquet(RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd)
{
	int N_mat = md.size * md.size - 1;

	pd.k1 = new double[N_mat];
	pd.k2 = new double[N_mat];
	pd.k3 = new double[N_mat];
	pd.k4 = new double[N_mat];
	pd.val = new double[N_mat];
	pd.tmp = new double[N_mat];
	pd.tmp_drv = new double[N_mat];

	pd.floquet = new MKL_Complex16[md.size * md.size * md.size * md.size];
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
	delete[] pd.tmp_drv;
}

void free_prop_data_floquet(RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd)
{
	delete[] pd.k1;
	delete[] pd.k2;
	delete[] pd.k3;
	delete[] pd.k4;
	delete[] pd.val;
	delete[] pd.tmp; 
	delete[] pd.tmp_drv;

	delete[] pd.floquet;
}

void dump_prop_data_std(RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd)
{
	cout << "Saving specs to files" << endl << endl;

	// Add here saving regular characteristics
}

void dump_prop_data_floquet(RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd)
{
	int size_xtd = md.size * md.size;
	MKL_Complex16 * evals = new MKL_Complex16[size_xtd];

	int info = LAPACKE_zgeev(
		LAPACK_ROW_MAJOR,
		'N',
		'N',
		size_xtd,
		pd.floquet,
		size_xtd,
		evals,
		NULL,
		size_xtd,
		NULL,
		size_xtd
	);

	/* Check for convergence */
	if (info > 0) {
		printf("The algorithm failed to compute eigenvalues.\n");
		exit(1);
	}

	cout << "Saving specs to files" << endl << endl;

	string fn;

	fn = "floquet_evals" + file_name_suffix(cp, 4);
	cout << "Saving floquet_evals to file:" << endl << fn << endl << endl;
	save_complex_data(fn, evals, size_xtd, 16, false);

	delete[] evals;
}

