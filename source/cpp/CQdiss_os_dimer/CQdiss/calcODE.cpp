#include "calcODE.h"
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
void init_conditions_floquet(dcomplex *mtx, RunParam &rp, ConfigParam &cp, MainData &md, int state_id)
{
	for (int s_id_1 = 0; s_id_1 < md.size; s_id_1++)
	{
		for (int s_id_2 = 0; s_id_2 < md.size; s_id_2++)
		{
			mtx[s_id_1 * md.size + s_id_2].re = 0.0;
			mtx[s_id_1 * md.size + s_id_2].im = 0.0;
		}
	}

	mtx[state_id].re = 1.0;
	mtx[state_id].im = 0.0;
}

void multMatVec_real(crsMatrix *mat, double * x, double * res)
{
	char trans = 'n';
	double *value = (double *)(mat->Value);
	mkl_dcsrgemv(&trans, &(mat->N), value, mat->RowIndex, mat->Col, x, res);
}
void multMatVec_complex(crsMatrix *mat, dcomplex * x, dcomplex * res)
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
			res[i].re += v1.re * v2.re;
			res[i].im = 0.0;
		}
	}
}

void calcVectValue_real(
	double t,
	double h,
	Model* m,
	double* x,
	double* res,
	double* tmp1,
	double* tmp2
)
{
	int i;
	int N_mat = m->N_mat;
	crsMatrix * Gs = m->Gs;
	crsMatrix * Qs_drv = m->Qs_drv;
	double    * Ks = (double *)(m->Ks);

	double drv_ampl = -m->conf.drv_ampl;
	double drv_freq = m->conf.drv_freq;
	double drv_phase = m->conf.drv_phase;

	int drv_type = m->conf.drv_type;

	multMatVec_real(Gs, x, tmp1);
	multMatVec_real(Qs_drv, x, tmp2);

	for (i = 0; i < N_mat; i++)
	{
		if (drv_type == 1)
		{
			res[i] = (tmp1[i] + drv_ampl * sin(drv_freq * t + drv_phase) * tmp2[i] - Ks[i]) * h;
		}
		else if (drv_type == 0)
		{
			if (sin(drv_freq * t + drv_phase) < 0.0)
			{
				res[i] = (tmp1[i] + drv_ampl * (-1.0) * tmp2[i] - Ks[i]) * h;
			}
			else
			{
				res[i] = (tmp1[i] + drv_ampl * (+1.0) * tmp2[i] - Ks[i]) * h;
			}
		}
	}
}

void initRhoODE(Model *m, RunParam &rp, ConfigParam &cp, MainData &md)
{
	int N = m->N;
	int N_mat = m->N_mat;
	dcomplex  * RhoF = m->RhoF;
	dcomplex * psi = new dcomplex[(N + 1)*(N + 1)];
	init_conditions(psi, rp, cp, md);

	if (rp.debug == 1)
	{
		string fn = "rho_ini" + file_name_suffix(cp, 4);
		save_complex_data(fn, (MKL_Complex16 *)psi, (N + 1) * (N + 1), 16, false);
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
			RhoF[k].re += psi[(i) * (N + 1) + (j)].re * val;
			RhoF[k].re += psi[(j) * (N + 1) + (i)].re * val;
			k++;

			RhoF[k].re += psi[(i) * (N + 1) + (j)].im * (-val);
			RhoF[k].re += psi[(j) * (N + 1) + (i)].im * (+val);
			k++;
		}
	}

	for (int i = 0; i < N; i++)
	{
		val = 1.0 / sqrt((double)((i + 1) * (i + 2)));
		for (int j = 0; j <= i; j++)
		{
			RhoF[k].re += psi[j * (N + 1) + j].re * val;
		}
		RhoF[k].re -= psi[(i + 1) * (N + 1) + (i + 1)].re * val * (i + 1);
		k++;
	}

	delete[] psi;
}

void initRhoODE_floquet(Model *m, RunParam &rp, ConfigParam &cp, MainData &md, int state_id)
{
	int N = m->N;
	int N_mat = m->N_mat;
	dcomplex  * RhoF = m->RhoF;
	dcomplex * psi = new dcomplex[(N + 1)*(N + 1)];
	init_conditions_floquet(psi, rp, cp, md, state_id);

	if (rp.debug == 1)
	{
		string fn = "rho_ini_floquet_" + std::to_string(state_id) + file_name_suffix(cp, 4);
		save_complex_data(fn, (MKL_Complex16 *)psi, (N + 1) * (N + 1), 16, false);
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
			RhoF[k].re += psi[(i) * (N + 1) + (j)].re * val;
			RhoF[k].re += psi[(j) * (N + 1) + (i)].re * val;
			k++;

			RhoF[k].re += psi[(i) * (N + 1) + (j)].im * (-val);
			RhoF[k].re += psi[(j) * (N + 1) + (i)].im * (+val);
			k++;
		}
	}

	for (int i = 0; i < N; i++)
	{
		val = 1.0 / sqrt((double)((i + 1) * (i + 2)));
		for (int j = 0; j <= i; j++)
		{
			RhoF[k].re += psi[j * (N + 1) + j].re * val;
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
		diff.re = abs(m->prevRhoF[i].re - m->RhoF[i].re);
		diff.im = abs(m->prevRhoF[i].im - m->RhoF[i].im);
		if (max_diff.re < diff.re)max_diff.re = diff.re;
		if (max_diff.im < diff.im)max_diff.im = diff.im;
	}

	return max_diff;
}

void calcODE_trans(Model *m, RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd)
{
	int i;
	int N_mat = m->N_mat;
	double * RhoF = (double *)(m->RhoF);

	double eps = 1.0e-10;

	double* k1 = pd.k1;
	double* k2 = pd.k2;
	double* k3 = pd.k3;
	double* k4 = pd.k4;
	double* val = pd.val;
	double* tmp = pd.tmp;
	double* tmp_drv = pd.tmp_drv;

	double step = md.step;

	double time = 0;
	dcomplex diff_it;
	before(m);
	for (int period = 0; period < cp.num_periods_trans; period++)
	{
		for (i = 0; i < N_mat; i++)
		{
			m->prevRhoF[i].re = m->RhoF[i].re;
			m->prevRhoF[i].im = m->RhoF[i].im;
		}

		for (int step_id = 0; step_id < cp.num_steps; step_id++)
		{
			time = period * md.T + step_id * md.step;

			calcVectValue_real(time, step, m, RhoF, k1, tmp, tmp_drv);
			for (i = 0; i < N_mat; i++)
			{
				val[i] = RhoF[i] + k1[i] / 2.0;
			}
			calcVectValue_real(time + step / 2.0, step, m, val, k2, tmp, tmp_drv);
			for (i = 0; i < N_mat; i++)
			{
				val[i] = RhoF[i] + k2[i] / 2.0;
			}
			calcVectValue_real(time + step / 2.0, step, m, val, k3, tmp, tmp_drv);
			for (i = 0; i < N_mat; i++)
			{
				val[i] = RhoF[i] + k3[i];
			}
			calcVectValue_real(time + step, step, m, val, k4, tmp, tmp_drv);

			for (i = 0; i < N_mat; i++)
			{
				RhoF[i] += (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
			}
		}

		dcomplex diff_it = calcDiffIter(m);
		printf("diff on %d is %0.16le %0.16le \n", period + 1, diff_it.re, diff_it.im);
	}
	after(m);
}
void calcODE_std(Model *m, RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd)
{
	int i;
	int N_mat = m->N_mat;
	double * RhoF = (double *)(m->RhoF);
	double * prevRhoF = (double *)(m->prevRhoF);

	double eps = 1.0e-10;

	for (i = 0; i < N_mat; i++)
	{
		prevRhoF[i] = RhoF[i];
	}

	double * k1 = pd.k1;
	double * k2 = pd.k2;
	double * k3 = pd.k3;
	double * k4 = pd.k4;
	double * val = pd.val;
	double * tmp = pd.tmp;
	double* tmp_drv = pd.tmp_drv;

	double step = md.step;

	double time = 0;

	int dump_id = 0;
	calcRho(m);
	characteristics_std(m, rp, cp, md, pd, dump_id);
	dump_id++;

	before(m);
	for (int period = 0; period < cp.num_periods_obser; period++)
	{
		for (int step_id = 0; step_id < cp.num_steps; step_id++)
		{
			time = period * md.T + step_id * md.step;

			calcVectValue_real(time, step, m, RhoF, k1, tmp, tmp_drv);
			for (i = 0; i < N_mat; i++)
			{
				val[i] = RhoF[i] + k1[i] / 2.0;
			}
			calcVectValue_real(time + step / 2.0, step, m, val, k2, tmp, tmp_drv);
			for (i = 0; i < N_mat; i++)
			{
				val[i] = RhoF[i] + k2[i] / 2.0;
			}
			calcVectValue_real(time + step / 2.0, step, m, val, k3, tmp, tmp_drv);
			for (i = 0; i < N_mat; i++)
			{
				val[i] = RhoF[i] + k3[i];
			}
			calcVectValue_real(time + step, step, m, val, k4, tmp, tmp_drv);

			for (i = 0; i < N_mat; i++)
			{
				RhoF[i] += (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
			}
		}

		time = (period + 1) * md.T;

		if ((period + 1) == pd.dump_periods[dump_id])
		{
			if (rp.ipp == 1)
			{
				cout << endl << "Dump period: " << (period + 1) << endl;
			}

			after(m);
			calcRho(m);
			characteristics_std(m, rp, cp, md, pd, dump_id);
			before(m);
			dump_id++;
		}
	}

	after(m);
}
void calcODE_floquet(Model *m, RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd)
{
	int i;
	int N_mat = m->N_mat;
	double * RhoF = (double *)(m->RhoF);
	double * prevRhoF = (double *)(m->prevRhoF);

	double * k1 = pd.k1;
	double * k2 = pd.k2;
	double * k3 = pd.k3;
	double * k4 = pd.k4;
	double * val = pd.val;
	double * tmp = pd.tmp;
	double * tmp_drv = pd.tmp_drv;

	double step = md.step;

	int size = md.size;
	int size_xtd = size * size;

	for (int fl_id = 0; fl_id < size_xtd; fl_id++)
	{
		cout << "floquet: " << fl_id << endl;

		initRhoODE_floquet(m, rp, cp, md, fl_id);
		double time = 0;

		before(m);
		for (int step_id = 0; step_id < cp.num_steps; step_id++)
		{
			time = step_id * md.step;

			calcVectValue_real(time, step, m, RhoF, k1, tmp, tmp_drv);
			for (i = 0; i < N_mat; i++)
			{
				val[i] = RhoF[i] + k1[i] / 2.0;
			}
			calcVectValue_real(time + step / 2.0, step, m, val, k2, tmp, tmp_drv);
			for (i = 0; i < N_mat; i++)
			{
				val[i] = RhoF[i] + k2[i] / 2.0;
			}
			calcVectValue_real(time + step / 2.0, step, m, val, k3, tmp, tmp_drv);
			for (i = 0; i < N_mat; i++)
			{
				val[i] = RhoF[i] + k3[i];
			}
			calcVectValue_real(time + step, step, m, val, k4, tmp, tmp_drv);

			for (i = 0; i < N_mat; i++)
			{
				RhoF[i] += (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
			}
		}
		after(m);
		calcRho(m);

		for (int i = 0; i < size; i++)
		{
			int s = m->Rho->RowIndex[i];
			int f = m->Rho->RowIndex[i + 1];
			for (int k = s; k < f; k++)
			{
				int j = m->Rho->Col[k];
				int index = (i * size + j) * size_xtd + fl_id;
				pd.floquet[index].real = m->Rho->Value[k].re;
				pd.floquet[index].imag = m->Rho->Value[k].im;
			}
		}
	}
}

void before(Model *m)
{
	complex_to_real(m->Gs->Value, m->Gs->NZ);
	complex_to_real(m->Qs_drv->Value, m->Qs_drv->NZ);
	complex_to_real(m->Ks, m->N_mat);
	complex_to_real(m->RhoF, m->N_mat);
	toOneBase(*(m->Gs));
	toOneBase(*(m->Qs_drv));
}
void after(Model *m)
{
	toZeroBase(*(m->Gs));
	toZeroBase(*(m->Qs_drv));
	real_to_complex(m->Gs->Value, m->Gs->NZ);
	real_to_complex(m->Qs_drv->Value, m->Qs_drv->NZ);
	real_to_complex(m->Ks, m->N_mat);
	real_to_complex(m->RhoF, m->N_mat);
}

void complex_to_real(dcomplex *mat, int N)
{
	double *value = (double *)(mat);
	for (int i = 0; i <N; i++)
	{
		value[i] = mat[i].re;
	}
}
void real_to_complex(dcomplex *mat, int N)
{
	double *value = (double *)(mat);
	for (int i = N - 1; i >= 0; i--)
	{
		mat[i].re = value[i];
		mat[i].im = 0;
	}
}
