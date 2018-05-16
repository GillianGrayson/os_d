#include "calcODE.h"
#include "calcRho.h"
#include "characteristics.h"

void init_conditions(dcomplex *mtx, RunParam &rp, ConfigParam &cp, MainData &md)
{
	for (int s_id_1 = 0; s_id_1 < md.size; s_id_1++)
	{
		for (int s_id_2 = 0; s_id_2 < md.size; s_id_2++)
		{
			mtx[s_id_1 * md.size + s_id_2].re = 0.0;
			mtx[s_id_1 * md.size + s_id_2].im = 0.0;
		}
	}

	if (cp.int_ist == 0)
	{
		mtx[cp.int_isi * md.size + cp.int_isi].re = 1.0;
		mtx[cp.int_isi * md.size + cp.int_isi].im = 0.0;
	}
	else if (cp.int_ist == 1)
	{
		for (int s_id_1 = 0; s_id_1 < md.size; s_id_1++)
		{
			mtx[s_id_1 * md.size + s_id_1].re = 1.0 / md.size;
			mtx[s_id_1 * md.size + s_id_1].im = 0.0;
		}
	}
	else
	{
		stringstream msg;
		msg << "wrong int_ist value: " << cp.int_ist << endl;
		Error(msg.str());
	}
}

void multMatVec(crsMatrix *mat, dcomplex * x, dcomplex * res)
{
	int i, j, s, f;
	for (i = 0; i < mat->N; i++)
	{
		s = mat->RowIndex[i];
		f = mat->RowIndex[i + 1];
		res[i].re = 0.0;
		res[i].im = 0.0;
		for (j = s; j < f; j++)
		{
			dcomplex v1 = mat->Value[j];
			dcomplex v2 = x[mat->Col[j]];
			res[i].re += v1.re * v2.re;// - v1.im * v2.im;
			res[i].im = 0.0;//+= v1.re * v2.im + v1.im * v2.re;
		}
	}
}

void calcVectValue_t0(double h, Model * m, dcomplex *x, dcomplex * res, dcomplex * tmp)
{
	int i;
	int N_mat = m->N_mat;
	crsMatrix * G_0_s = m->G_0_s;
	dcomplex  * Ks = m->Ks;

	multMatVec(G_0_s, x, tmp);

	for (i = 0; i < N_mat; i++)
	{
		res[i].re = (tmp[i].re - Ks[i].re) * h;
		res[i].im = 0.0;
	}
}

void calcVectValue_t1(double h, Model * m, dcomplex *x, dcomplex * res, dcomplex * tmp)
{
	int i;
	int N_mat = m->N_mat;
	crsMatrix * G_1_s = m->G_1_s;
	dcomplex  * Ks = m->Ks;

	multMatVec(G_1_s, x, tmp);

	for (i = 0; i < N_mat; i++)
	{
		res[i].re = (tmp[i].re - Ks[i].re) * h;
		res[i].im = 0.0;
	}
}

void initRho(Model *m, RunParam &rp, ConfigParam &cp, MainData &md)
{
	int N = m->N;
	int N_mat = m->N_mat;
	dcomplex  * RhoF = m->RhoF;
	dcomplex * psi = new dcomplex[(N + 1)*(N + 1)];
	init_conditions(psi, rp, cp, md);

	if (rp.debug == 1)
	{
		string fn = rp.path + "rho_ini" + file_name_suffix(cp, 4);
		save_complex_data(fn, (MKL_Complex16 *)psi, (N + 1)*(N + 1), 16, false);
	}

	for (int i = 0; i < N_mat; i++)
	{
		RhoF[i].re = 0.0;
		RhoF[i].im = 0.0;
	}
	int k = 0;
	double val = 1.0 / sqrt(2.0);
	for (int i = 0; i < N + 1; i++)
	{
		for (int j = i + 1; j < N + 1; j++)
		{
			RhoF[k].re += psi[(i)* (N + 1) + (j)].re * val;
			RhoF[k].re += psi[(j)* (N + 1) + (i)].re * val;
			k++;

			RhoF[k].re += psi[(i)* (N + 1) + (j)].im * (-val);
			RhoF[k].re += psi[(j)* (N + 1) + (i)].im * (+val);
			k++;
		}
	}

	for (int i = 0; i < N; i++)
	{
		val = 1.0 / sqrt((double)((i + 1)* (i + 2)));
		for (int j = 0; j <= i; j++)
		{
			RhoF[k].re += psi[j * (N + 1) + j].re *val;
		}
		RhoF[k].re -= psi[(i + 1) * (N + 1) + (i + 1)].re * val * (i + 1);
		k++;
	}

	delete[] psi;
}

dcomplex calcDiffIter(Model *m)
{
	dcomplex max_diff, diff;
	max_diff.re = 0.0;
	max_diff.im = 0.0;
	for (int i = 0; i < m->N_mat; i++)
	{
		diff.re = fabs(m->prevRhoF[i].re - m->RhoF[i].re);
		diff.im = fabs(m->prevRhoF[i].im - m->RhoF[i].im);
		if (max_diff.re < diff.re)max_diff.re = diff.re;
		if (max_diff.im < diff.im)max_diff.im = diff.im;
	}

	return max_diff;
}

void calcODE_trans(Model *m, RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd)
{
	int i;
	int N_mat = m->N_mat;
	dcomplex * RhoF = m->RhoF;

	double eps = 1.0e-10;

	for (i = 0; i < N_mat; i++)
	{
		m->prevRhoF[i] = RhoF[i];
	}

	dcomplex * k1 = pd.k1;
	dcomplex * k2 = pd.k2;
	dcomplex * k3 = pd.k3;
	dcomplex * k4 = pd.k4;
	dcomplex * val = pd.val;
	dcomplex * tmp = pd.tmp;

	double step_t_0 = md.step_t_0;
	double step_t_1 = md.step_t_1;


	for (int period = 1; period <= cp.num_periods_trans; period++)
	{
		for (int t_0_step_id = 0; t_0_step_id < cp.num_steps_t_0; t_0_step_id++)
		{
			calcVectValue_t0(step_t_0, m, RhoF, k1, tmp);
			for (i = 0; i < N_mat; i++)
			{
				val[i].re = RhoF[i].re + k1[i].re / 2.0;
				val[i].im = RhoF[i].im + k1[i].im / 2.0;
			}
			calcVectValue_t0(step_t_0, m, val, k2, tmp);
			for (i = 0; i < N_mat; i++)
			{
				val[i].re = RhoF[i].re + k2[i].re / 2.0;
				val[i].im = RhoF[i].im + k2[i].im / 2.0;
			}
			calcVectValue_t0(step_t_0, m, val, k3, tmp);
			for (i = 0; i < N_mat; i++)
			{
				val[i].re = RhoF[i].re + k3[i].re;
				val[i].im = RhoF[i].im + k3[i].im;
			}
			calcVectValue_t0(step_t_0, m, val, k4, tmp);

			for (i = 0; i < N_mat; i++)
			{
				RhoF[i].re += (k1[i].re + 2.0 * k2[i].re + 2.0 * k3[i].re + k4[i].re) / 6.0;
				RhoF[i].im += (k1[i].im + 2.0 * k2[i].im + 2.0 * k3[i].im + k4[i].im) / 6.0;
			}
		}

		for (int t_1_step_id = 0; t_1_step_id < cp.num_steps_t_1; t_1_step_id++)
		{
			calcVectValue_t1(step_t_1, m, RhoF, k1, tmp);
			for (i = 0; i < N_mat; i++)
			{
				val[i].re = RhoF[i].re + k1[i].re / 2.0;
				val[i].im = RhoF[i].im + k1[i].im / 2.0;
			}
			calcVectValue_t1(step_t_1, m, val, k2, tmp);
			for (i = 0; i < N_mat; i++)
			{
				val[i].re = RhoF[i].re + k2[i].re / 2.0;
				val[i].im = RhoF[i].im + k2[i].im / 2.0;
			}
			calcVectValue_t1(step_t_1, m, val, k3, tmp);
			for (i = 0; i < N_mat; i++)
			{
				val[i].re = RhoF[i].re + k3[i].re;
				val[i].im = RhoF[i].im + k3[i].im;
			}
			calcVectValue_t1(step_t_1, m, val, k4, tmp);

			for (i = 0; i < N_mat; i++)
			{
				RhoF[i].re += (k1[i].re + 2.0 * k2[i].re + 2.0 * k3[i].re + k4[i].re) / 6.0;
				RhoF[i].im += (k1[i].im + 2.0 * k2[i].im + 2.0 * k3[i].im + k4[i].im) / 6.0;
			}
		}
	}
}

void calcODE_std(Model *m, RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd)
{
	int i;
	int N_mat = m->N_mat;
	dcomplex * RhoF = m->RhoF;

	double eps = 1.0e-10;

	for (i = 0; i < N_mat; i++)
	{
		m->prevRhoF[i] = RhoF[i];
	}

	dcomplex * k1 = pd.k1;
	dcomplex * k2 = pd.k2;
	dcomplex * k3 = pd.k3;
	dcomplex * k4 = pd.k4;
	dcomplex * val = pd.val;
	dcomplex * tmp = pd.tmp;

	double step_t_0 = md.step_t_0;
	double step_t_1 = md.step_t_1;

	double curr_time = 0.0;
	int dump_id = 0;

	calcRho(m);

	characteristics_std(m, rp, cp, md, pd, dump_id);
	dump_id++;

	for (int period = 1; period <= cp.num_periods_obser; period++)
	{
		for (int t_0_step_id = 0; t_0_step_id < cp.num_steps_t_0; t_0_step_id++)
		{
			calcVectValue_t0(step_t_0, m, RhoF, k1, tmp);
			for (i = 0; i < N_mat; i++)
			{
				val[i].re = RhoF[i].re + k1[i].re / 2.0;
				val[i].im = RhoF[i].im + k1[i].im / 2.0;
			}
			calcVectValue_t0(step_t_0, m, val, k2, tmp);
			for (i = 0; i < N_mat; i++)
			{
				val[i].re = RhoF[i].re + k2[i].re / 2.0;
				val[i].im = RhoF[i].im + k2[i].im / 2.0;
			}
			calcVectValue_t0(step_t_0, m, val, k3, tmp);
			for (i = 0; i < N_mat; i++)
			{
				val[i].re = RhoF[i].re + k3[i].re;
				val[i].im = RhoF[i].im + k3[i].im;
			}
			calcVectValue_t0(step_t_0, m, val, k4, tmp);

			for (i = 0; i < N_mat; i++)
			{
				RhoF[i].re += (k1[i].re + 2.0 * k2[i].re + 2.0 * k3[i].re + k4[i].re) / 6.0;
				RhoF[i].im += (k1[i].im + 2.0 * k2[i].im + 2.0 * k3[i].im + k4[i].im) / 6.0;
			}
		}

		for (int t_1_step_id = 0; t_1_step_id < cp.num_steps_t_1; t_1_step_id++)
		{
			calcVectValue_t1(step_t_1, m, RhoF, k1, tmp);
			for (i = 0; i < N_mat; i++)
			{
				val[i].re = RhoF[i].re + k1[i].re / 2.0;
				val[i].im = RhoF[i].im + k1[i].im / 2.0;
			}
			calcVectValue_t1(step_t_1, m, val, k2, tmp);
			for (i = 0; i < N_mat; i++)
			{
				val[i].re = RhoF[i].re + k2[i].re / 2.0;
				val[i].im = RhoF[i].im + k2[i].im / 2.0;
			}
			calcVectValue_t1(step_t_1, m, val, k3, tmp);
			for (i = 0; i < N_mat; i++)
			{
				val[i].re = RhoF[i].re + k3[i].re;
				val[i].im = RhoF[i].im + k3[i].im;
			}
			calcVectValue_t1(step_t_1, m, val, k4, tmp);

			for (i = 0; i < N_mat; i++)
			{
				RhoF[i].re += (k1[i].re + 2.0 * k2[i].re + 2.0 * k3[i].re + k4[i].re) / 6.0;
				RhoF[i].im += (k1[i].im + 2.0 * k2[i].im + 2.0 * k3[i].im + k4[i].im) / 6.0;
			}
		}

		curr_time = period * md.T;

		if (period == pd.dump_periods[dump_id])
		{

			if (rp.ipp == 1)
			{
				cout << endl << "Dump period: " << period << endl;
			}

			calcRho(m);
			characteristics_std(m, rp, cp, md, pd, dump_id);

			dump_id++;
		}
	}
}

void calcODE_deep(Model *m, RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd)
{
	string fn;

	int dupm_step_t_0 = cp.num_steps_t_0 / cp.int_dn;
	int dupm_step_t_1 = cp.num_steps_t_1 / cp.int_dn;

	int i;
	int N_mat = m->N_mat;
	dcomplex * RhoF = m->RhoF;

	double eps = 1.0e-10;

	for (i = 0; i < N_mat; i++)
	{
		m->prevRhoF[i] = RhoF[i];
	}

	dcomplex * k1 = pd.k1;
	dcomplex * k2 = pd.k2;
	dcomplex * k3 = pd.k3;
	dcomplex * k4 = pd.k4;
	dcomplex * val = pd.val;
	dcomplex * tmp = pd.tmp;

	double step_t_0 = md.step_t_0;
	double step_t_1 = md.step_t_1;

	double curr_time = 0.0;
	int dump_id = 0;

	calcRho(m);

	if (rp.debug == 1)
	{
		fn = rp.path + "rho" + "_" + to_string(dump_id) + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, m->Rho, 16, false);
	}

	characteristics_deep(m, rp, cp, md, pd, dump_id);
	dump_id++;

	for (int period = 1; period <= cp.num_periods_obser; period++)
	{
		for (int t_0_step_id = 0; t_0_step_id < cp.num_steps_t_0; t_0_step_id++)
		{
			calcVectValue_t0(step_t_0, m, RhoF, k1, tmp);
			for (i = 0; i < N_mat; i++)
			{
				val[i].re = RhoF[i].re + k1[i].re / 2.0;
				val[i].im = RhoF[i].im + k1[i].im / 2.0;
			}
			calcVectValue_t0(step_t_0, m, val, k2, tmp);
			for (i = 0; i < N_mat; i++)
			{
				val[i].re = RhoF[i].re + k2[i].re / 2.0;
				val[i].im = RhoF[i].im + k2[i].im / 2.0;
			}
			calcVectValue_t0(step_t_0, m, val, k3, tmp);
			for (i = 0; i < N_mat; i++)
			{
				val[i].re = RhoF[i].re + k3[i].re;
				val[i].im = RhoF[i].im + k3[i].im;
			}
			calcVectValue_t0(step_t_0, m, val, k4, tmp);

			for (i = 0; i < N_mat; i++)
			{
				RhoF[i].re += (k1[i].re + 2.0 * k2[i].re + 2.0 * k3[i].re + k4[i].re) / 6.0;
				RhoF[i].im += (k1[i].im + 2.0 * k2[i].im + 2.0 * k3[i].im + k4[i].im) / 6.0;
			}

			if (rp.debug == 1)
			{
				calcRho(m);

				fn = rp.path + "rho" + "_" + to_string(dump_id) + file_name_suffix(cp, 4);
				save_sparse_complex_mtx(fn, m->Rho, 16, false);
			}

			if (t_0_step_id % dupm_step_t_0 == 0)
			{
				calcRho(m);
				characteristics_deep(m, rp, cp, md, pd, dump_id);
				dump_id++;
			}
		}

		for (int t_1_step_id = 0; t_1_step_id < cp.num_steps_t_1; t_1_step_id++)
		{
			calcVectValue_t1(step_t_1, m, RhoF, k1, tmp);
			for (i = 0; i < N_mat; i++)
			{
				val[i].re = RhoF[i].re + k1[i].re / 2.0;
				val[i].im = RhoF[i].im + k1[i].im / 2.0;
			}
			calcVectValue_t1(step_t_1, m, val, k2, tmp);
			for (i = 0; i < N_mat; i++)
			{
				val[i].re = RhoF[i].re + k2[i].re / 2.0;
				val[i].im = RhoF[i].im + k2[i].im / 2.0;
			}
			calcVectValue_t1(step_t_1, m, val, k3, tmp);
			for (i = 0; i < N_mat; i++)
			{
				val[i].re = RhoF[i].re + k3[i].re;
				val[i].im = RhoF[i].im + k3[i].im;
			}
			calcVectValue_t1(step_t_1, m, val, k4, tmp);

			for (i = 0; i < N_mat; i++)
			{
				RhoF[i].re += (k1[i].re + 2.0 * k2[i].re + 2.0 * k3[i].re + k4[i].re) / 6.0;
				RhoF[i].im += (k1[i].im + 2.0 * k2[i].im + 2.0 * k3[i].im + k4[i].im) / 6.0;
			}

			if (rp.debug == 1)
			{
				calcRho(m);

				fn = rp.path + "rho" + "_" + to_string(dump_id) + file_name_suffix(cp, 4);
				save_sparse_complex_mtx(fn, m->Rho, 16, false);
			}

			if (t_1_step_id % dupm_step_t_1 == 0)
			{
				calcRho(m);
				characteristics_deep(m, rp, cp, md, pd, dump_id);
				dump_id++;
			}
		}

		curr_time = period * md.T;

		if (rp.ipp == 1)
		{
			cout << endl << "Dump period: " << period << endl;
		}
	}
}