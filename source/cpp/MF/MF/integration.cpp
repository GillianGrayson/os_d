#include "integration.h"

void right_part(ConfigParam &cp, double * ks, double * x, double time)
{
	double driving = sin(cp.omega * time + cp.phase);

	if (cp.mt == 1)
	{
		driving *= cp.A;
	}
	else
	{
		driving = cp.A * (2.0 * (driving > 0.0) - 1.0);
	}

	ks[0] = 2.0 * cp.J * sin(x[1]) + 4.0 * cp.gamma * cos(x[1]) * cos(x[0]);
	ks[1] = 2.0 * cp.J * cos(x[0]) * cos(x[1]) / sin(x[0]) - 2.0 * cp.E - 2.0 * driving + cp.U * cos(x[0]) - 4.0 * cp.gamma * sin(x[1]) / sin(x[0]);
	ks[2] = 1.0;
}

void right_part_lpn(ConfigParam &cp, double * ks_lpn, double * x_lpn, double time, double * x)
{
	double F1_x = (-4.0 * cp.gamma * cos(x[1]) * sin(x[0]));
	double F1_y = (2.0 * cp.J * cos(x[1]) - 4.0 * cp.gamma * cos(x[0]) * sin(x[1]));
	double F1_t = 0.0;
	double F2_x = (-2.0 * cp.J * cos(x[1]) / pow(sin(x[0]), 2.0) - cp.U * sin(x[0]) + 4.0 * cp.gamma * sin(x[1]) * cos(x[0]) / pow(sin(x[0]), 2.0));
	double F2_y = (-2.0 * cp.J * cos(x[0]) / sin(x[0]) * sin(x[1]) - 4.0 * cp.gamma * cos(x[1]) / sin(x[0]));
	double F2_t = 0.0;
	if (cp.mt == 1)
	{
		F2_t = -2.0 * cp.A * cp.omega *  cos(cp.omega * time + cp.phase);
	}
	double F3_x = 0.0;
	double F3_y = 0.0;
	double F3_t = 0.0;

	ks_lpn[0] = F1_x * x_lpn[0] + F1_y * x_lpn[1] + F1_t * x_lpn[2];
	ks_lpn[1] = F2_x * x_lpn[0] + F2_y * x_lpn[1] + F2_t * x_lpn[2];
	ks_lpn[2] = F3_x * x_lpn[0] + F3_y * x_lpn[1] + F3_t * x_lpn[2];
}

void upd_arg(int size, double * x_arg, double * x, double * ks, double coeff)
{
	for (int st_id = 0; st_id < size; st_id++)
	{
		x_arg[st_id] = x[st_id] + coeff * ks[st_id];
	}
}

void rk_final(int size, double * x, double * k1s, double * k2s, double * k3s, double * k4s, double step)
{
	for (int st_id = 0; st_id < size; st_id++)
	{
		x[st_id] += (k1s[st_id] + 2.0 * k2s[st_id] + 2.0 * k3s[st_id] + k4s[st_id]) * step / 6.0;
	}
}

void rk_step(ConfigParam &cp, MainData &md)
{
	right_part(cp, md.k1s, md.data, md.time);
	upd_arg(md.size, md.args, md.data, md.k1s, md.step * 0.5);

	md.time += md.step * 0.5;

	right_part(cp, md.k2s, md.args, md.time);
	upd_arg(md.size, md.args, md.data, md.k2s, md.step * 0.5);

	right_part(cp, md.k3s, md.args, md.time);
	upd_arg(md.size, md.args, md.data, md.k3s, md.step * 1.0);

	md.time += md.step * 0.5;

	right_part(cp, md.k4s, md.args, md.time);

	rk_final(md.size, md.data, md.k1s, md.k2s, md.k3s, md.k4s, md.step);
}

void rk_step_lpn(ConfigParam &cp, MainData &md)
{
	// ========= K1 ==========
	right_part(cp, md.k1s, md.data, md.time);
	for (int lpn_id = 0; lpn_id < md.num_lpn; lpn_id++)
	{
		right_part_lpn(cp, md.k1s_lpn[lpn_id], md.data_lpn[lpn_id], md.time, md.data);
	}

	upd_arg(md.size, md.args, md.data, md.k1s, md.step * 0.5);
	for (int lpn_id = 0; lpn_id < md.num_lpn; lpn_id++)
	{
		upd_arg(md.size, md.args_lpn[lpn_id], md.data_lpn[lpn_id], md.k1s_lpn[lpn_id], md.step * 0.5);
	}

	md.time += md.step * 0.5;

	// ========= K2 ==========
	right_part(cp, md.k2s, md.args, md.time);
	for (int lpn_id = 0; lpn_id < md.num_lpn; lpn_id++)
	{
		right_part_lpn(cp, md.k2s_lpn[lpn_id], md.args_lpn[lpn_id], md.time, md.args);
	}

	upd_arg(md.size, md.args, md.data, md.k2s, md.step * 0.5);
	for (int lpn_id = 0; lpn_id < md.num_lpn; lpn_id++)
	{
		upd_arg(md.size, md.args_lpn[lpn_id], md.data_lpn[lpn_id], md.k2s_lpn[lpn_id], md.step * 0.5);
	}

	// ========= K3 ==========
	right_part(cp, md.k3s, md.args, md.time);
	for (int lpn_id = 0; lpn_id < md.num_lpn; lpn_id++)
	{
		right_part_lpn(cp, md.k3s_lpn[lpn_id], md.args_lpn[lpn_id], md.time, md.args);
	}

	upd_arg(md.size, md.args, md.data, md.k3s, md.step * 1.0);
	for (int lpn_id = 0; lpn_id < md.num_lpn; lpn_id++)
	{
		upd_arg(md.size, md.args_lpn[lpn_id], md.data_lpn[lpn_id], md.k3s_lpn[lpn_id], md.step * 1.0);
	}

	md.time += md.step * 0.5;

	// ========= K4 ==========
	right_part(cp, md.k4s, md.args, md.time);
	for (int lpn_id = 0; lpn_id < md.num_lpn; lpn_id++)
	{
		right_part_lpn(cp, md.k4s_lpn[lpn_id], md.args_lpn[lpn_id], md.time, md.args);
	}

	rk_final(md.size, md.data, md.k1s, md.k2s, md.k3s, md.k4s, md.step);
	for (int lpn_id = 0; lpn_id < md.num_lpn; lpn_id++)
	{
		rk_final(md.size, md.data_lpn[lpn_id], md.k1s_lpn[lpn_id], md.k2s_lpn[lpn_id], md.k3s_lpn[lpn_id], md.k4s_lpn[lpn_id], md.step);
	}

	gsorth_lpn(md);
}

void int_period(ConfigParam &cp, MainData &md, int per_id)
{
	for (int step_id = 0; step_id < cp.num_steps; step_id++)
	{
		md.time = per_id * cp.T + step_id * md.step;
		md.data[2] = md.time;
		rk_step(cp, md);
	}
}

void int_period_lpn(ConfigParam &cp, MainData &md, int per_id)
{
	for (int step_id = 0; step_id < cp.num_steps; step_id++)
	{
		md.time = per_id * cp.T + step_id * md.step;
		md.data[2] = md.time;
		rk_step_lpn(cp, md);
	}
}

void int_period_cd(ConfigParam &cp, MainData &md, int per_id)
{
	int curr_point_id = 0;

	for (int step_id = 0; step_id < cp.num_steps; step_id++)
	{
		if (step_id == md.cd_ti[curr_point_id])
		{
			md.cd_obs = cos(md.data[0]) + 1.0;
			md.cd_rd[per_id][curr_point_id] = md.cd_obs;

			curr_point_id++;
		}

		md.time = per_id * cp.T + step_id * md.step;
		md.data[2] = md.time;
		rk_step(cp, md);
	}
}

void int_period_cd_d(ConfigParam &cp, MainData &md, int per_id)
{
	int curr_point_id = 0;
	int global_point_id = 0;

	for (int step_id = 0; step_id < cp.num_steps; step_id++)
	{
		global_point_id = per_id * cp.num_steps + step_id;

		for (int cd_st_id = 0; cd_st_id < md.cd_size; cd_st_id++)
		{
			curr_point_id = global_point_id - cd_st_id;

			if (curr_point_id >= 0 && curr_point_id < md.cd_M)
			{
				md.cd_obs = cos(md.data[0]) + 1.0;

				md.cd_rd[curr_point_id][cd_st_id] = md.cd_obs;
			}
		}

		md.time = per_id * cp.T + step_id * md.step;
		md.data[2] = md.time;
		rk_step(cp, md);
	}
}

void int_period_cd_sd(ConfigParam &cp, MainData &md, int per_id)
{
	int dump_point_id = 0;
	int curr_point_id = 0;
	int global_point_id = 0;

	for (int step_id = 0; step_id < cp.num_steps; step_id++)
	{
		if (step_id % md.cd_ps == 0)
		{
			global_point_id = per_id * cp.cd_nppp + dump_point_id;

			dump_point_id++;

			for (int cd_st_id = 0; cd_st_id < md.cd_size; cd_st_id++)
			{
				curr_point_id = global_point_id - cd_st_id;

				if (curr_point_id >= 0 && curr_point_id < md.cd_M)
				{
					md.cd_obs = cos(md.data[0]) + 1.0;

					md.cd_rd[curr_point_id][cd_st_id] = md.cd_obs;
				}
			}
		}

		md.time = per_id * cp.T + step_id * md.step;
		md.data[2] = md.time;
		rk_step(cp, md);
	}
}

void int_trans_proc(ConfigParam &cp, MainData &md)
{
	for (int per_id = 0; per_id < cp.npt; per_id++)
	{
		int_period(cp, md, per_id);
	}
}