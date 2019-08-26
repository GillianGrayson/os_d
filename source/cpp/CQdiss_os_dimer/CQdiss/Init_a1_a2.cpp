#include "Init_a1_a2.h"
#include "initH.h"
#include <math.h>


crsMatrix * createA1mat(int N)
{
	crsMatrix * mat;
	mat = new crsMatrix(N + 1, N + 1);
	mat->RowIndex[0] = 0;

	for (int i = 0; i < mat->N; i++)
	{
		mat->Col[i] = i;
		mat->Value[i].re = 2.0 * i - N;
		mat->RowIndex[i + 1] = mat->RowIndex[i] + 1;
	}

	return mat;
}

crsMatrix * createA2mat(int N)
{
	int i, j;
	crsMatrix * mat;
	mat = new crsMatrix(N + 1, N * 2);
	mat->RowIndex[0] = 0;
	mat->Col[0] = 1;
	mat->RowIndex[1] = 1;
	for (i = 1; i < N; i++)
		mat->RowIndex[i + 1] = mat->RowIndex[i] + 2;
	mat->RowIndex[i + 1] = mat->RowIndex[i] + 1;
	double val;
	val = N;
	mat->Value[0].im = sqrt(val);
	for (i = 1; i < N; i++)
	{
		j = mat->RowIndex[i];
		mat->Col[j] = i - 1;
		mat->Col[j + 1] = i + 1;
		val = (N - i + 1) * (i + 0);
		mat->Value[j].im = -sqrt(val);
		val = (N - i + 0)  * (i + 1);
		mat->Value[j + 1].im = sqrt(val);
	}
	j = mat->RowIndex[i];
	mat->Col[j] = i - 1;
	val = (N - i + 1)  * (i + 0);
	mat->Value[j].im = -sqrt(val);

	return mat;
}

void init_diss(Model * m, RunParam &rp, ConfigParam &cp, MainData &md)
{
	int N = m->N;

	crsMatrix * A1 = createA1mat(N);
	crsMatrix * A2 = createA2mat(N);
	if (rp.debug == 1)
	{
		string fn = "A1" + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, A1, 16, false);
		fn = "A2" + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, A2, 16, false);
	}

	int N_mat = m->N_mat;

	crsMatrix * a1_mat = new crsMatrix(N_mat, N_mat);
	crsMatrix * a2_mat = new crsMatrix(N_mat, N_mat);

	to_F_basis(A1, a1_mat);
	to_F_basis(A2, a2_mat);

	if (rp.debug == 1)
	{
		string fn = "A1F" + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, a1_mat, 16, false);
		fn = "A2F" + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, a2_mat, 16, false);
	}

	crsMatrix * a1_i_a2_mat = new crsMatrix(N_mat, a2_mat->NZ + a1_mat->NZ);
	int k = 0;
	int k1, k2;
	for (int i = 0; i < N_mat; i++)
	{
		a1_i_a2_mat->RowIndex[i] = k;
		int c1 = a1_mat->RowIndex[i + 1] - a1_mat->RowIndex[i];
		int c2 = a2_mat->RowIndex[i + 1] - a2_mat->RowIndex[i];
		if (c1 > 0)
		{
			a1_i_a2_mat->Value[k] = a1_mat->Value[a1_mat->RowIndex[i]];
		}
		if (c2 > 0)
		{
			a1_i_a2_mat->Value[k].re += a2_mat->Value[a2_mat->RowIndex[i]].im;
			a1_i_a2_mat->Value[k].im -= a2_mat->Value[a2_mat->RowIndex[i]].re;
		}
		if (c1 + c2 > 0)
		{
			a1_i_a2_mat->Col[k] = i;
			k++;
		}
	}
	a1_i_a2_mat->RowIndex[N_mat] = k;
	a1_i_a2_mat->NZ = k;

	m->l_mat = a1_i_a2_mat;
	delete a1_mat;
	delete a2_mat;

	if (rp.debug == 1)
	{
		string fn = "diss_mtx" + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, a1_i_a2_mat, 16, false);
	}

	delete A1;
	delete A2;
}


