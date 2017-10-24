#include "experiments.h"

void only_data_exp(RunParam &rp, ConfigParam &cp)
{
	for (int U_id = 0; U_id < rp.U_num; U_id++)
	{
		cp.U = rp.U_start + double(U_id) * rp.U_shift;

		for (int seed = rp.seed_start; seed < rp.seed_start + rp.seed_num; seed++)
		{
			cp.seed = seed;

			cout << "U = " << cp.U << " seed = " << seed << endl;
			string fn_data = rp.path + "nu" + file_name_suffix(cp, 4);
			string fn_time = rp.path + "time" + file_name_suffix(cp, 4);

			double * times = new double[cp.np + 1];
			double * nus = new double[cp.np + 1];

			MainData dt;

			init_main_data(cp, dt);
			init_cond(rp, cp, dt);

			int_trans_proc(cp, dt);
			times[0] = dt.time;
			nus[0] = dt.data[0];

			for (int per_id = 0; per_id < cp.np; per_id++)
			{
				int_period(cp, dt, cp.npt + per_id);
				times[per_id + 1] = dt.time;
				nus[per_id + 1] = dt.data[0];
			}

			write_double_data(fn_data, nus, cp.np + 1, 16, 0);
			write_double_data(fn_time, times, cp.np + 1, 16, 0);

			delete_main_data(dt);

			delete[] times;
			delete[] nus;
		}
	}
}

void lpn_exp(RunParam &rp, ConfigParam &cp)
{
	for (int U_id = 0; U_id < rp.U_num; U_id++)
	{
		cp.U = rp.U_start + double(U_id) * rp.U_shift;

		for (int seed = rp.seed_start; seed < rp.seed_start + rp.seed_num; seed++)
		{
			cp.seed = seed;

			cout << "U = " << cp.U << " seed = " << seed << endl;
			string fn_data = rp.path + "nu" + file_name_suffix(cp, 4);
			string fn_time = rp.path + "time" + file_name_suffix(cp, 4);
			string fn_exps_lpn = rp.path + "exps_lpn" + file_name_suffix(cp, 4);

			MainData md;

			init_main_data(cp, md);
			init_lpn_data(cp, md);

			double * times = new double[cp.np + 1];
			double * nus = new double[cp.np + 1];
			double ** exps_lpn = new double *[md.num_lpn];
			for (int lpn_id = 0; lpn_id < md.num_lpn; lpn_id++)
			{
				exps_lpn[lpn_id] = new double[cp.np + 1];
			}

			init_cond(rp, cp, md);
			init_cond_lpn(cp, md);

			int_trans_proc(cp, md);

			times[0] = md.time;
			nus[0] = md.data[0];
			for (int lpn_id = 0; lpn_id < md.num_lpn; lpn_id++)
			{
				exps_lpn[lpn_id][0] = 0.0;
			}

			for (int per_id = 0; per_id < cp.np; per_id++)
			{
				int_period_lpn(cp, md, per_id);

				times[per_id + 1] = md.time;
				nus[per_id + 1] = md.data[0];
				for (int lpn_id = 0; lpn_id < md.num_lpn; lpn_id++)
				{
					exps_lpn[lpn_id][per_id + 1] = md.exps_lpn[lpn_id];
				}
			}

			write_double_data(fn_data, nus, cp.np + 1, 16, 0);
			write_double_data(fn_time, times, cp.np + 1, 16, 0);
			write_2d_double_data(fn_exps_lpn, exps_lpn, md.num_lpn, cp.np + 1, 16, 0);

			delete_main_data(md);
			delete_lpn_data(md);

			delete[] times;
			delete[] nus;
			for (int lpn_id = 0; lpn_id < md.num_lpn; lpn_id++)
			{
				delete[] exps_lpn[lpn_id];
			}
			delete[] exps_lpn;
		}
	}
}