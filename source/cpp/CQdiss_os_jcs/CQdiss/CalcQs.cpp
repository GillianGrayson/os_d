#include "CalcQs.h"

void calc_Q_0(Model * m)
{
	int N_mat = m->N_mat;
	dcomplex * h_0 = m->h_0;

	crsMatrix ** f_mat = m->f_H_mat;
	crsMatrix * Q_0 = new crsMatrix(m->N_mat, 1);
	crsMatrix * resSum;

	Q_0->Col[0] = 0;
	Q_0->RowIndex[0] = 0;
	for (int i = 1; i <= N_mat; i++)
	{
		Q_0->RowIndex[i] = 0;
	}

	for (int i = 0; i < N_mat; i++)
	{
		if ((h_0[i].re != 0.0) || (h_0[i].im != 0.0))
		{
			resSum = new crsMatrix;
			SparseMKLAdd(*Q_0, h_0[i], *(f_mat[i]), *resSum);
			delete Q_0;
			Q_0 = resSum;
		}
	}

	m->Q_0 = Q_0;
}

void calc_Q_1(Model * m)
{
	int N_mat = m->N_mat;
	dcomplex * h_1 = m->h_1;

	crsMatrix ** f_mat = m->f_H_mat;
	crsMatrix * Q_1 = new crsMatrix(m->N_mat, 1);
	crsMatrix * resSum;

	Q_1->Col[0] = 0;
	Q_1->RowIndex[0] = 0;
	for (int i = 1; i <= N_mat; i++)
	{
		Q_1->RowIndex[i] = 0;
	}

	for (int i = 0; i < N_mat; i++)
	{
		if ((h_1[i].re != 0.0) || (h_1[i].im != 0.0))
		{
			resSum = new crsMatrix;
			SparseMKLAdd(*Q_1, h_1[i], *(f_mat[i]), *resSum);
			delete Q_1;
			Q_1 = resSum;
		}
	}

	m->Q_1 = Q_1;
}
