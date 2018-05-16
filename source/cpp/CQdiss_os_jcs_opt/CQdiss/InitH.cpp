#include "initH.h"
#include <math.h>
#include <algorithm>
#include <mkl.h>

void to_F_basis_for_zeros(crsMatrix * Mat, crsMatrix * vec)
{
	int N = Mat->N;
	int N_mat = N * N - 1;
	int cnt = 0;
	int *mask = new int[N_mat];
	int *col = new int[N_mat];
	dcomplex *value = new dcomplex[N_mat];
	for (int i = 0; i < N_mat; i++)
	{
		mask[i] = 0;
		value[i].re = 0.0;
		value[i].im = 0.0;
	}

	dcomplex sum;
	sum.re = 0.0;
	sum.im = 0.0;
	for (int i = 0; i < N; i++)
	{
		int start = Mat->RowIndex[i];
		int finish = Mat->RowIndex[i + 1];
		for (int j = start; j < finish; j++)
		{
			int k = Mat->Col[j];

			if (k != i)
			{
				int ii = i;
				int jj = k;
				int z = -1;
				if (ii > jj) {
					ii = k;
					jj = i;
					z = 1;
				}

				int index = ((N - 1 + N - ii) * ii) + jj - ii - 1;
				if (mask[index] != 1)
				{
					col[cnt + 0] = index + 0;
					col[cnt + 1] = index + 1;
					cnt += 2;
				}

				mask[index] = 1;
				mask[index + 1] = 1;

				value[index].re += Mat->Value[j].re / sqrt(2.0);
				value[index].im += Mat->Value[j].im / sqrt(2.0);

				value[index + 1].re += z * Mat->Value[j].im / sqrt(2.0);
				value[index + 1].im += -z * Mat->Value[j].re / sqrt(2.0);

			}
			else
			{
				if (k != 0)
				{
					int index = N * (N - 1) + k - 1;
					mask[index] = 1;
					double value_d = 1.0 / sqrt((double)((k + 0)* (k + 1)));
					value[index].re = sum.re * value_d;
					value[index].re -= Mat->Value[j].re * (k + 0) * value_d;
					value[index].im = sum.im * value_d;
					value[index].im -= Mat->Value[j].im * (k + 0) * value_d;
					col[cnt] = index;
					cnt++;
				}

				sum.re += Mat->Value[j].re;
				sum.im += Mat->Value[j].im;
			}
		}
	}

	std::sort(col, col + cnt);

	vec->NZ = cnt;
	for (int i = 0; i < vec->N + 1; i++)
	{
		vec->RowIndex[i] = 0;
	}
	for (int i = 0; i < cnt; i++)
	{
		vec->Col[i] = col[i];
		vec->Value[i] = value[col[i]];
		vec->RowIndex[col[i] + 1] ++;
	}
	for (int i = 0; i < vec->N; i++)
	{
		vec->RowIndex[i + 1] = vec->RowIndex[i] + vec->RowIndex[i + 1];
	}

	delete[] mask;
	delete[] col;
	delete[] value;
}
void to_F_basis(crsMatrix * Mat, crsMatrix * vec)
{
	int N = Mat->N;
	int N_mat = N * N - 1;
	int cnt = 0;
//	int cnt_diag = 0;
	int *mask = new int[N_mat];
	int *col = new int[N_mat];
	dcomplex *value = new dcomplex[N_mat];
	for (int i = 0; i < N_mat; i++)
	{
		mask[i] = 0;
		value[i].re = 0.0;
		value[i].im = 0.0;
	}

	dcomplex sum;
	sum.re = 0.0;
	sum.im = 0.0;
	for (int i = 0; i < N; i++)
	{
		int diag_ex = 0;
		int start = Mat->RowIndex[i];
		int finish = Mat->RowIndex[i + 1];
		for (int j = start; j < finish; j++)
		{
			int k = Mat->Col[j];

			if (k != i)
			{
				int ii = i;
				int jj = k;
				int z = -1;
				if (ii > jj) {
					ii = k;
					jj = i;
					z = 1;
				}

				int index = ((N - 1 + N - ii) * ii) + (jj - ii - 1) * 2;
				if (mask[index] != 1)
				{
					col[cnt + 0] = index + 0;
					col[cnt + 1] = index + 1;
					cnt += 2;
				}

				mask[index] = 1;
				mask[index + 1] = 1;

				value[index].re += Mat->Value[j].re / sqrt(2.0);
				value[index].im += Mat->Value[j].im / sqrt(2.0);

				value[index + 1].re += z * Mat->Value[j].im / sqrt(2.0);
				value[index + 1].im += -z * Mat->Value[j].re / sqrt(2.0);

			}
			else
			{
				diag_ex = 1;
				if (k != 0)
				{
					int index = N * (N - 1) + k - 1;
					mask[index] = 1;
					double value_d = 1.0 / sqrt((double)((k + 0)* (k + 1)));
					value[index].re = sum.re * value_d;
					value[index].re -= Mat->Value[j].re * (k + 0) * value_d;
					value[index].im = sum.im * value_d;
					value[index].im -= Mat->Value[j].im * (k + 0) * value_d;
					col[cnt] = index;
					cnt++;
				}

				sum.re += Mat->Value[j].re;
				sum.im += Mat->Value[j].im;
			}
		}
		if (diag_ex == 0)
		{
			if (i != 0)
			{
				int index = N * (N - 1) + i - 1;
				mask[index] = 1;
				double value_d = 1.0 / sqrt((double)((i + 0)* (i + 1)));
				value[index].re = sum.re * value_d;
				value[index].im = sum.im * value_d;
				col[cnt] = index;
				cnt++;
			}
		}
	}

	for (int i = 0; i < cnt; i++)
	{
		if ((value[col[i]].re == 0.0) && (value[col[i]].im == 0.0))
		{
			cnt--;
			col[i] = col[cnt];
			i--;
		}
	}
	std::sort(col, col + cnt);

	vec->NZ = cnt;
	for (int i = 0; i < vec->N + 1; i++)
	{
		vec->RowIndex[i] = 0;
	}
	for (int i = 0; i < cnt; i++)
	{
		vec->Col[i] = col[i];
		vec->Value[i] = value[col[i]];
		vec->RowIndex[col[i] + 1] ++;
	}
	for (int i = 0; i < vec->N; i++)
	{
		vec->RowIndex[i + 1] = vec->RowIndex[i] + vec->RowIndex[i + 1];
	}

	delete[] mask;
	delete[] col;
	delete[] value;
}

crsMatrix * create_H_0_matrix(Model * m, RunParam &rp, ConfigParam &cp, MainData &md)
{
	int N = m->N;

	crsMatrix * a_std = create_a_std_matrix(m);
	crsMatrix * a_std_copy = new crsMatrix(*a_std);
	crsMatrix * a_dag = create_a_dag_matrix(m);
	crsMatrix * a_dag_copy = new crsMatrix(*a_dag);

	if (rp.debug == 1)
	{
		string fn = rp.path + "a_std" + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, a_std, 16, false);

		fn = rp.path + "a_dag" + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, a_dag, 16, false);
	}

	crsMatrix * tmp_0 = new crsMatrix;
	crsMatrix * tmp_1 = new crsMatrix;
	crsMatrix * tmp_2 = new crsMatrix;

	SparseMKLMult(*a_dag, *a_dag_copy, *tmp_0);
	SparseMKLMult(*a_std, *a_std_copy, *tmp_1);
	SparseMKLMult(*tmp_0, *tmp_1, *tmp_2);

	double coeff = 0.5 * 1.0 / pow(cp.prm_alpha, 3);

	scalar_mult(tmp_2, coeff);

	crsMatrix * H_0 = new crsMatrix(*tmp_2);

	if (rp.debug == 1)
	{
		string fn = rp.path + "H_0" + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, H_0, 16, false);

		fn = rp.path + "a_std" + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, a_std, 16, false);

		fn = rp.path + "a_dag" + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, a_dag, 16, false);
	}

	delete a_std;
	delete a_std_copy;
	delete a_dag;
	delete a_dag_copy;

	delete tmp_0;
	delete tmp_1;
	delete tmp_2;

	return H_0;
}
void init_h_0_vector(Model * m, RunParam &rp, ConfigParam &cp, MainData &md)
{
	crsMatrix * H_0 = create_H_0_matrix(m, rp, cp, md), *res;
	m->H_0 = H_0;
	int N = m->N;
	crsMatrix * h_0 = m->h_0;

	to_F_basis(H_0, h_0);
}

crsMatrix * create_H_1_matrix(Model * m, RunParam &rp, ConfigParam &cp, MainData &md)
{
	int N = m->N;

	crsMatrix * a_std = create_a_std_matrix(m);
	crsMatrix * a_dag = create_a_dag_matrix(m);

	crsMatrix * tmp_0 = new crsMatrix();

	dcomplex beta = { -1.0, 0.0 };

	SparseMKLAdd(*a_dag, beta, *a_std, *tmp_0);

	int NZ = tmp_0->NZ;
	dcomplex * Value = tmp_0->Value;
	double real = 0.0;
	double imag = 0.0;
	for (int nz_id = 0; nz_id < NZ; nz_id++)
	{
		real = Value[nz_id].re;
		imag = Value[nz_id].im;

		Value[nz_id].re = -imag;
		Value[nz_id].im = real;
	}

	crsMatrix * H_1 = new crsMatrix(*tmp_0);

	if (rp.debug == 1)
	{
		string fn = rp.path + "H_1" + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, H_1, 16, false);
	}

	delete a_std;
	delete a_dag;

	delete tmp_0;

	return H_1;
}
void init_h_1_vector(Model * m, RunParam &rp, ConfigParam &cp, MainData &md)
{
	crsMatrix * H_1 = create_H_1_matrix(m, rp, cp, md), *res;
	m->H_1 = H_1;
	int N = m->N;

	crsMatrix * h_1 = m->h_1;

	to_F_basis(H_1, h_1);
}

void init_H0(Model * m, RunParam &rp, ConfigParam &cp, MainData &md)
{
	int N_mat = m->N_mat;
	crsMatrix * H0 = create_H_0_matrix(m, rp, cp, md);
	m->H0 = H0;
}
void init_H1(Model * m, RunParam &rp, ConfigParam &cp, MainData &md)
{
	int N_mat = m->N_mat;
	crsMatrix * H1 = create_H_1_matrix(m, rp, cp, md);
	m->H1 = H1;
}


