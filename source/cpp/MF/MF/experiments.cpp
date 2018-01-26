#include "experiments.h"

void basic_exp(RunParam &rp, ConfigParam &cp)
{
	for (int U_id = 0; U_id < rp.U_num; U_id++)
	{
		cp.U = rp.U_start + double(U_id) * rp.U_shift;

		for (int seed = rp.seed_start; seed < rp.seed_start + rp.seed_num; seed++)
		{
			cp.seed = seed;

			cout << "U = " << cp.U << " seed = " << seed << endl;
			string fn_data = rp.path + "data" + file_name_suffix(rp, cp, 4);

			MainData dt;

			init_main_data(cp, dt);
			init_cond(rp, cp, dt);
			init_evo_data(cp, dt);

			int_trans_proc(cp, dt);

			for (int per_id = 0; per_id < cp.np; per_id++)
			{
				int_period(cp, dt, cp.npt + per_id);
			}

			write_2d_double_data(fn_data, dt.data_evo, dt.size, dt.num_dumps, 16, 0);

			delete_main_data(dt);
			delete_evo_data(dt);
		}
	}
}

void lpn_fin_exp(RunParam &rp, ConfigParam &cp)
{
	for (int U_id = 0; U_id < rp.U_num; U_id++)
	{
		cp.U = rp.U_start + double(U_id) * rp.U_shift;

		for (int seed = rp.seed_start; seed < rp.seed_start + rp.seed_num; seed++)
		{
			cp.seed = seed;

			cout << "U = " << cp.U << " seed = " << seed << endl;
			string fn_exps_lpn = rp.path + "exps_lpn" + file_name_suffix(rp, cp, 4);

			MainData md;

			init_main_data(cp, md);
			init_lpn_data(cp, md);

			double * exps_lpn = new double[md.num_lpn];
			for (int lpn_id = 0; lpn_id < md.num_lpn; lpn_id++)
			{
				exps_lpn[lpn_id] = 0.0;
			}

			init_cond(rp, cp, md);
			init_cond_lpn(cp, md);
			init_evo_data(cp, md);

			int_trans_proc(cp, md);

			for (int per_id = 0; per_id < cp.np; per_id++)
			{
				int_period_lpn(cp, md, per_id);
			}

			for (int lpn_id = 0; lpn_id < md.num_lpn; lpn_id++)
			{
				exps_lpn[lpn_id] = md.exps_lpn[lpn_id];
			}

			write_double_data(fn_exps_lpn, exps_lpn, md.num_lpn, 16, 0);

			delete_main_data(md);
			delete_lpn_data(md);
			delete_evo_data(md);

			delete[] exps_lpn;
		}
	}
}

void basic_and_lpn_fin_exp(RunParam &rp, ConfigParam &cp)
{
	for (int U_id = 0; U_id < rp.U_num; U_id++)
	{
		cp.U = rp.U_start + double(U_id) * rp.U_shift;

		for (int seed = rp.seed_start; seed < rp.seed_start + rp.seed_num; seed++)
		{
			cp.seed = seed;

			cout << "U = " << cp.U << " seed = " << seed << endl;
			string fn_exps_lpn = rp.path + "exps_lpn" + file_name_suffix(rp, cp, 4);
			string fn_data = rp.path + "data" + file_name_suffix(rp, cp, 4);

			MainData md;

			init_main_data(cp, md);
			init_lpn_data(cp, md);

			double * exps_lpn = new double[md.num_lpn];
			for (int lpn_id = 0; lpn_id < md.num_lpn; lpn_id++)
			{
				exps_lpn[lpn_id] = 0.0;
			}

			init_cond(rp, cp, md);
			init_cond_lpn(cp, md);
			init_evo_data(cp, md);

			int_trans_proc(cp, md);

			for (int per_id = 0; per_id < cp.np; per_id++)
			{
				int_period_lpn(cp, md, per_id);
			}

			for (int lpn_id = 0; lpn_id < md.num_lpn; lpn_id++)
			{
				exps_lpn[lpn_id] = md.exps_lpn[lpn_id];
			}

			write_2d_double_data(fn_data, md.data_evo, md.size, md.num_dumps, 16, 0);
			write_double_data(fn_exps_lpn, exps_lpn, md.num_lpn, 16, 0);

			delete_main_data(md);
			delete_lpn_data(md);
			delete_evo_data(md);

			delete[] exps_lpn;
		}
	}
}

void cd_exp(RunParam &rp, ConfigParam &cp)
{
	double * cd_eps = new double[rp.cd_eps_num];
	
	for (int U_id = 0; U_id < rp.U_num; U_id++)
	{
		cp.U = rp.U_start + double(U_id) * rp.U_shift;

		for (int seed = rp.seed_start; seed < rp.seed_start + rp.seed_num; seed++)
		{
			cp.seed = seed;

			for (int eps_id = 0; eps_id < rp.cd_eps_num; eps_id++)
			{
				cd_eps[eps_id] = rp.cd_eps * pow(10.0, (1.0 / double(rp.cd_eps_ndpd)) * double(eps_id));

				cp.cd_eps = cd_eps[eps_id];

				cout << "U = " << cp.U << " seed = " << cp.seed << " eps = " << cp.cd_eps << endl;

				MainData md;

				init_main_data(cp, md);
				init_cd_data(cp, md);

				print_cd_info(md);

				init_cond(rp, cp, md);
				init_evo_data(cp, md);

				int_trans_proc(cp, md);

				for (int per_id = 0; per_id < cp.np; per_id++)
				{
					int_period_cd(cp, md, per_id);
				}

				calc_ci(cp, md);

				string fn_ci = rp.path + "ci" + file_name_suffix(rp, cp, 4);
				write_double_data(fn_ci, &md.cd_ci, 1, 16, 0);

				delete_main_data(md);
				delete_cd_data(md);
				delete_evo_data(md);
			}
		}
	}

	delete[] cd_eps;
}

void cd_d_exp(RunParam &rp, ConfigParam &cp)
{
	double * cd_eps = new double[rp.cd_eps_num];

	for (int U_id = 0; U_id < rp.U_num; U_id++)
	{
		cp.U = rp.U_start + double(U_id) * rp.U_shift;

		for (int seed = rp.seed_start; seed < rp.seed_start + rp.seed_num; seed++)
		{
			cp.seed = seed;

			for (int eps_id = 0; eps_id < rp.cd_eps_num; eps_id++)
			{
				cd_eps[eps_id] = rp.cd_eps * pow(10.0, (1.0 / double(rp.cd_eps_ndpd)) * double(eps_id));

				cp.cd_eps = cd_eps[eps_id];

				cout << "U = " << cp.U << " seed = " << cp.seed << " eps = " << cp.cd_eps << endl;

				MainData md;

				init_main_data(cp, md);
				init_cd_d_data(cp, md);

				print_cd_info(md);

				init_cond(rp, cp, md);
				init_evo_data(cp, md);

				int_trans_proc(cp, md);

				for (int per_id = 0; per_id < cp.np; per_id++)
				{
					int_period_cd_d(cp, md, per_id);
				}

				calc_ci(cp, md);

				string fn_ci = rp.path + "ci" + file_name_suffix(rp, cp, 4);
				write_double_data(fn_ci, &md.cd_ci, 1, 16, 0);

				delete_main_data(md);
				delete_cd_d_data(md);
				delete_evo_data(md);
			}
		}
	}

	delete[] cd_eps;
}

void cd_sd_exp(RunParam &rp, ConfigParam &cp)
{
	double * cd_eps = new double[rp.cd_eps_num];

	for (int U_id = 0; U_id < rp.U_num; U_id++)
	{
		cp.U = rp.U_start + double(U_id) * rp.U_shift;

		for (int seed = rp.seed_start; seed < rp.seed_start + rp.seed_num; seed++)
		{
			cp.seed = seed;

			for (int eps_id = 0; eps_id < rp.cd_eps_num; eps_id++)
			{
				cd_eps[eps_id] = rp.cd_eps * pow(10.0, (1.0 / double(rp.cd_eps_ndpd)) * double(eps_id));

				cp.cd_eps = cd_eps[eps_id];

				cout << "U = " << cp.U << " seed = " << cp.seed << " eps = " << cp.cd_eps << endl;

				MainData md;

				init_main_data(cp, md);
				init_cd_sd_data(cp, md);

				print_cd_info(md);

				init_cond(rp, cp, md);
				init_evo_data(cp, md);

				int_trans_proc(cp, md);

				for (int per_id = 0; per_id < cp.np; per_id++)
				{
					int_period_cd_sd(cp, md, per_id);
				}

				calc_ci(cp, md);

				string fn_ci = rp.path + "ci" + file_name_suffix(rp, cp, 4);
				write_double_data(fn_ci, &md.cd_ci, 1, 16, 0);

				delete_main_data(md);
				delete_cd_sd_data(md);
				delete_evo_data(md);
			}
		}
	}

	delete[] cd_eps;
}

