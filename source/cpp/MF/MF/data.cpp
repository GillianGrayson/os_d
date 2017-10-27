#include "data.h"

void init_main_data(ConfigParam &cp, MainData &md)
{
	md.size = 3;

	md.step = (2.0 * PI) / (cp.omega * cp.num_steps);

	md.time = 0;

	md.data = new double[md.size];
	md.args = new double[md.size];
	md.k1s = new  double[md.size];
	md.k2s = new  double[md.size];
	md.k3s = new  double[md.size];
	md.k4s = new  double[md.size];

	for (int v_id = 0; v_id < md.size; v_id++)
	{
		md.data[v_id] = 0.0;
		md.args[v_id] = 0.0;
		md.k1s[v_id] = 0.0;
		md.k2s[v_id] = 0.0;
		md.k3s[v_id] = 0.0;
		md.k4s[v_id] = 0.0;
	}
}

void init_lpn_data(ConfigParam &cp, MainData &md)
{
	md.num_lpn = md.size;

	md.data_lpn = new double *[md.num_lpn];
	md.args_lpn = new double *[md.num_lpn];
	md.k1s_lpn = new double *[md.num_lpn];
	md.k2s_lpn = new double *[md.num_lpn];
	md.k3s_lpn = new double *[md.num_lpn];
	md.k4s_lpn = new double *[md.num_lpn];

	md.norms_lpn = new double[md.num_lpn];
	md.exps_lpn = new double[md.num_lpn];
	md.rvm_lpn = new double[md.num_lpn];

	for (int lpn_id = 0; lpn_id < md.num_lpn; lpn_id++)
	{
		md.data_lpn[lpn_id] = new double[md.size];
		md.args_lpn[lpn_id] = new double[md.size];
		md.k1s_lpn[lpn_id] = new double[md.size];
		md.k2s_lpn[lpn_id] = new double[md.size];
		md.k3s_lpn[lpn_id] = new double[md.size];
		md.k4s_lpn[lpn_id] = new double[md.size];

		for (int st_id = 0; st_id < md.size; st_id++)
		{
			md.data_lpn[lpn_id][st_id] = 0.0;
			md.args_lpn[lpn_id][st_id] = 0.0;
			md.k1s_lpn[lpn_id][st_id] = 0.0;
			md.k2s_lpn[lpn_id][st_id] = 0.0;
			md.k3s_lpn[lpn_id][st_id] = 0.0;
			md.k4s_lpn[lpn_id][st_id] = 0.0;
		}

		md.norms_lpn[lpn_id] = 0.0;
		md.exps_lpn[lpn_id] = 0.0;
		md.rvm_lpn[lpn_id] = 0.0;
	}
}

void init_cd_data(ConfigParam &cp, MainData &md)
{
	md.cd_ps = cp.num_steps / cp.cd_dim;
	md.cd_M = cp.np;

	md.cd_size = cp.cd_dim;

	md.cd_obs = 0.0;

	md.cd_ti = new int[md.cd_size];
	for (int cd_st_id = 0; cd_st_id < md.cd_size; cd_st_id++)
	{
		md.cd_ti[cd_st_id] = cd_st_id * md.cd_ps;
	}

	md.cd_rd = new double *[md.cd_M];
	for (int cd_p_id = 0; cd_p_id < md.cd_M; cd_p_id++)
	{
		md.cd_rd[cd_p_id] = new double[md.cd_size];

		for (int cd_st_id = 0; cd_st_id < md.cd_size; cd_st_id++)
		{
			md.cd_rd[cd_p_id][cd_st_id] = 0.0;
		}
	}
}

void init_cd_d_data(ConfigParam &cp, MainData &md)
{
	md.cd_ps = 1;
	md.cd_size = cp.cd_dim;
	md.cd_M = cp.np * cp.num_steps - md.cd_size + 1;

	md.cd_obs = 0.0;

	md.cd_rd = new double *[md.cd_M];
	for (int cd_p_id = 0; cd_p_id < md.cd_M; cd_p_id++)
	{
		md.cd_rd[cd_p_id] = new double[md.cd_size];

		for (int cd_st_id = 0; cd_st_id < md.cd_size; cd_st_id++)
		{
			md.cd_rd[cd_p_id][cd_st_id] = 0.0;
		}
	}
}

void delete_main_data(MainData &md)
{
	delete_data(md.data);
	delete_data(md.args);
	delete_data(md.k1s);
	delete_data(md.k2s);
	delete_data(md.k3s);
	delete_data(md.k4s);
}

void delete_lpn_data(MainData &md)
{
	for (int lpn_id = 0; lpn_id < md.num_lpn; lpn_id++)
	{
		delete[] md.data_lpn[lpn_id];
		delete[] md.args_lpn[lpn_id];
		delete[] md.k1s_lpn[lpn_id];
		delete[] md.k2s_lpn[lpn_id];
		delete[] md.k3s_lpn[lpn_id];
		delete[] md.k4s_lpn[lpn_id];
	}

	delete[] md.data_lpn;
	delete[] md.args_lpn;
	delete[] md.k1s_lpn;
	delete[] md.k2s_lpn;
	delete[] md.k3s_lpn;
	delete[] md.k4s_lpn;

	delete[] md.norms_lpn;
	delete[] md.exps_lpn;
	delete[] md.rvm_lpn;
}

void delete_cd_data(MainData &md)
{
	delete[] md.cd_ti;
	
	for (int cd_p_id = 0; cd_p_id < md.cd_size; cd_p_id++)
	{
		delete[] md.cd_rd[cd_p_id];
	}
	delete[] md.cd_rd;
}

void delete_cd_d_data(MainData &md)
{
	for (int cd_p_id = 0; cd_p_id < md.cd_size; cd_p_id++)
	{
		delete[] md.cd_rd[cd_p_id];
	}
	delete[] md.cd_rd;
}

void init_cond(RunParam &rp, ConfigParam &cp, MainData &md)
{
	double left_border = -PI;
	double right_border = PI;
	int max_num_seeds = rp.max_num_seeds;

	VSLStreamStatePtr stream;
	vslNewStream(&stream, VSL_BRNG_MCG31, 77778888);
	vslLeapfrogStream(stream, cp.seed, max_num_seeds);
	vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, md.size, md.data, left_border, right_border);

	md.data[2] = 0.0;
}

void init_cond_lpn(ConfigParam &cp, MainData &md)
{
	for (int lpn_st_id = 0; lpn_st_id < md.num_lpn; lpn_st_id++)
	{
		md.data_lpn[lpn_st_id][lpn_st_id] = 1.0;
	}
}

void calc_norm_lpn(MainData &md, int lpn_id)
{
	double sum = 0.0;

	for (int st_id = 0; st_id < md.size; st_id++)
	{
		sum += (md.data_lpn[lpn_id][st_id] * md.data_lpn[lpn_id][st_id]);
	}

	md.norms_lpn[lpn_id] = sqrt(sum);
}

void normalization_lpn(MainData &md, int lpn_id)
{
	for (int st_id = 0; st_id < md.size; st_id++)
	{
		md.data_lpn[lpn_id][st_id] = md.data_lpn[lpn_id][st_id] / md.norms_lpn[lpn_id];
	}
}

void scalar_mult_lpn(MainData &md, double * mults, int lpn_id, int lpn_id_tmp)
{
	for (int st_id = 0; st_id < md.size; st_id++)
	{
		mults[lpn_id_tmp] += md.data_lpn[lpn_id][st_id] * md.data_lpn[lpn_id_tmp][st_id];
	}
}

void sub_lpn(MainData &md, double* mults, int lpn_id, int lpn_id_tmp)
{
	for (int st_id = 0; st_id < md.size; st_id++)
	{
		md.data_lpn[lpn_id][st_id] -= mults[lpn_id_tmp] * md.data_lpn[lpn_id_tmp][st_id];
	}
}

void gsorth_lpn(MainData &md)
{
	/* Compute first norm */
	int lpn_id = 0;
	calc_norm_lpn(md, lpn_id);
	/* Renorm - find first vector in new basis */
	normalization_lpn(md, lpn_id);

	/* Find  others basis vectors */
	for (int lpn_id = 1; lpn_id < md.num_lpn; lpn_id++)
	{
		/* Firstly find all nessesary scalar mults */
		double *mults = new double[md.num_lpn];

		for (int lpn_id_tmp = 0; lpn_id_tmp < lpn_id; lpn_id_tmp++)
		{
			mults[lpn_id_tmp] = 0.0;
			scalar_mult_lpn(md, mults, lpn_id, lpn_id_tmp);
		}

		/* Compute unnormed basis vector using scalar mults*/
		for (int lpn_id_tmp = 0; lpn_id_tmp < lpn_id; lpn_id_tmp++)
		{
			sub_lpn(md, mults, lpn_id, lpn_id_tmp);
		}

		delete[] mults;

		/* Find norm of new basis vector */
		calc_norm_lpn(md, lpn_id);

		/* Renorm new vector */
		normalization_lpn(md, lpn_id);
	}

	/* Calculate Lyapunov exps */

	for (int lpn_id = 0; lpn_id < md.num_lpn; lpn_id++)
	{
		md.rvm_lpn[lpn_id] += log(md.norms_lpn[lpn_id]);
		md.exps_lpn[lpn_id] = md.rvm_lpn[lpn_id] / md.time;
	}
}

void calc_ci(ConfigParam &cp, MainData &md)
{
	double * curr_diff = new double[md.cd_size];
	double curr_norm = 0.0;

	double integral = 0.0;

	for (int cd_p_id_1 = 0; cd_p_id_1 < md.cd_M; cd_p_id_1++)
	{
		for (int cd_p_id_2 = 0; cd_p_id_2 < md.cd_M; cd_p_id_2++)
		{
			if (cd_p_id_1 != cd_p_id_2)
			{
				for (int cd_st_id = 0; cd_st_id < md.cd_size; cd_st_id++)
				{
					curr_diff[cd_st_id] = md.cd_rd[cd_p_id_1][cd_st_id] - md.cd_rd[cd_p_id_2][cd_st_id];
				}

				curr_norm = calc_norm(curr_diff, md.cd_size);

				if (curr_norm < cp.cd_eps)
				{
					integral += 1.0;
				}
			}
		}
	}

	integral /= (double(md.cd_M) * double(md.cd_M - 1));

	md.cd_ci = integral;

	delete curr_diff;
}

