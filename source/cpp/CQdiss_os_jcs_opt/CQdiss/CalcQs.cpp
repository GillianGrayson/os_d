#include "CalcQs.h"
#include "sortTensor.h"
#include <string.h>
#include "coef_coord.h"
#include "f_d_sym.h"
#include "f_d_ch.h"

ulli countSelect(Tensor_Coordinates * f_ijk, crsMatrix *hMat, ulli from, ulli to)
{
	unsigned int cnt = 0;
	for (unsigned int i = from; i < to; i++)
	{
		int j = f_ijk->coord1[i];
		cnt += hMat->RowIndex[j + 1] - hMat->RowIndex[j];
	}
	return cnt;
}
void dataSelect(int N_mat, Tensor_Coordinates * f_ijk, crsMatrix *hMat, unsigned int from, unsigned int to, Tensor_Coordinates * res)
{
	unsigned int k = 0;
	unsigned int cnt = 0;
	for (unsigned int i = from; i < to; i++)
	{
		int j = f_ijk->coord1[i];
		cnt = hMat->RowIndex[j + 1] - hMat->RowIndex[j];
		if (cnt > 0)
		{
			dcomplex val1 = f_ijk->data[i];
			val1.im = val1.im;
			dcomplex val2 = hMat->Value[hMat->RowIndex[j]];
			res->coord1[k] = f_ijk->coord1[i];
			res->coord2[k] = f_ijk->coord2[i];
			res->coord3[k] = f_ijk->coord3[i];
			res->data[k].re = val1.re * val2.re - val1.im * val2.im;
			res->data[k].im = val1.re * val2.im + val1.re * val2.im;

			res->hash[k] = 0;
			unsigned long long int step = N_mat;
			res->hash[k] += res->coord2[k] + step * res->coord3[k];
			k++;
		}
	}
}

void calc_CooQs(int N_mat, Model *m, Tensor_Coordinates * f_ijk, crsMatrix *hMat, crsMatrix *&res)
{
	Tensor_Coordinates * select = new Tensor_Coordinates;
	int *hash = new int[N_mat];
	int *rowi = new int[N_mat + 1];
	dcomplex *hash_calc = new dcomplex[N_mat];
	ulli cnt;

	cnt = fijk_coord_sym(hMat, m->N + 1);
	if (cnt == 0){
		res = new crsMatrix(N_mat, 1);
		res->Value[0].re = 0;
		res->Value[0].im = 0;
		res->Col[0] = 0;
		res->NZ = 0;
		for (int i = 0; i <= N_mat; i++)
		{
			res->RowIndex[i] = 0;
		}
		return;
	}
	fijk_coord_ch(select, hMat, cnt + 1, m->N + 1);

	for (ulli i = 0; i < select->k; i++)
	{
		int j = select->coord1[i];
		dcomplex val1 = select->data[i];
		dcomplex val2 = hMat->Value[hMat->RowIndex[j]];
		select->data[i].re = val1.re * val2.re - val1.im * val2.im;
		select->data[i].im = val1.re * val2.im + val1.re * val2.im;

		ulli step = N_mat;
		select->hash[i] = select->coord2[i] + step * select->coord3[i];
	}

	sort_matrix(select);


	for (int i = 0; i < N_mat; i++)
	{
		rowi[i] = -1;
	}
	for (ulli i = 0; i < select->k; i++)
	{
		if (rowi[select->coord3[i]] == -1)
		{
			rowi[select->coord3[i]] = i;
		}
	}
	rowi[N_mat] = select->k;
	for (int i = N_mat - 1; i >0; i--)
	{
		if (rowi[i] == -1)
		{
			rowi[i] = rowi[i + 1];
		}
	}
	rowi[N_mat] = select->k;

	res = new crsMatrix(N_mat, 1);

	int countNotZero = 0;
	for (int j = 0; j < N_mat; j++)
	{
		res->RowIndex[j] = countNotZero;
		memset(hash, 0, sizeof(int) * N_mat);
		int start, finish;
		start = rowi[j];
		finish = rowi[j + 1];
		if (start >= 0)
		{
			for (int k = start; k < finish; k++)
			{
				if (hash[select->coord2[k]] == 0)
				{
					hash[select->coord2[k]] = 1;
					countNotZero++;
				}
			}
		}
	}
	res->RowIndex[N_mat] = countNotZero;
	res->setNZ(countNotZero);

	countNotZero = 0;
	for (int j = 0; j < N_mat; j++)
	{
		memset(hash, 0, sizeof(int) * N_mat);
		int start, finish;
		start = rowi[j];
		finish = rowi[j + 1];
		if (start >= 0)
		{
			for (int k = start; k < finish; k++)
			{
				if (hash[select->coord2[k]] == 0)
				{
					hash[select->coord2[k]] = 1;
					res->Col[countNotZero] = select->coord2[k];
					countNotZero++;
					hash_calc[select->coord2[k]] = select->data[k];
				}
				else
				{
					hash_calc[select->coord2[k]].re += select->data[k].re;
					hash_calc[select->coord2[k]].im += select->data[k].im;
				}
			}
		}
		start = res->RowIndex[j];
		finish = res->RowIndex[j + 1];
		if (start >= 0)
		{
			for (int k = start; k < finish; k++)
			{
				res->Value[k] = hash_calc[res->Col[k]];
			}
		}
	}

	crsMatrix * R = new crsMatrix(*res);
	Transpose(*res, *R, false);
	Transpose(*R, *res, false);

	delete R;
	delete[] hash;
	delete[] hash_calc;
	delete[] rowi;

	free_matrix(select);
	delete select;
}

void calc_Q_0(Model * m)
{
	int N_mat = m->N_mat;
	crsMatrix * h_0 = m->h_0;

	crsMatrix * Q_0;

	if (m->f_ijk == NULL)
	{
		m->f_ijk = new Tensor_Coordinates;
		fijk_coord(m->f_ijk, m->N + 1);
	}

	calc_CooQs(N_mat, m, m->f_ijk, h_0, Q_0);

	m->Q_0 = Q_0;
}
void calc_Q_1(Model * m)
{
	int N_mat = m->N_mat;
	crsMatrix * h_1 = m->h_1;

	crsMatrix * Q_1;

	if (m->f_ijk == NULL)
	{
		m->f_ijk = new Tensor_Coordinates;
		fijk_coord(m->f_ijk, m->N + 1);
	}

	calc_CooQs(N_mat, m, m->f_ijk, h_1, Q_1);

	m->Q_1 = Q_1;
}