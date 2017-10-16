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

			Data dt;

			init_main_data(cp, dt);
			init_cond(cp, dt);

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