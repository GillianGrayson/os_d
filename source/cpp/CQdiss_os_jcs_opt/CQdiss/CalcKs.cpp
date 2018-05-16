#include "CalcKs.h"
#include "coef_coord.h"
#include "sortTensor.h"
#include <string.h>

void calcKs(Model *m)
{
	int N = m->N;
	int N_mat = m->N_mat;
	dcomplex  *Ks = m->Ks;

	crsMatrix * l_mat = m->l_mat;

	m->f_ijk = new Tensor_Coordinates;
	fijk_coord(m->f_ijk, m->N + 1);
	Tensor_Coordinates * f_ijk = m->f_ijk;

	for (int i = 0; i < N_mat; i++)
	{
		Ks[i].re = 0.0;
		Ks[i].im = 0.0;
	}

	ulli uN = l_mat->N;
	ulli uN2 = uN * uN;
	for (ulli i = 0; i < f_ijk->k; i++)
	{
		int ii = f_ijk->coord1[i];
		int ji = f_ijk->coord2[i];
		int ki = f_ijk->coord3[i];
		f_ijk->hash[i] = ji * uN + ii;
	}
	sort_matrix(f_ijk);

	dcomplex * hash = new dcomplex[l_mat->N];

	for (ulli i = 0; i < f_ijk->k; i++)
	{
		int sji = f_ijk->coord2[i];

		memset(hash, 0, sizeof(dcomplex) * l_mat->N);
		int start = l_mat->RowIndex[sji];
		int finish = l_mat->RowIndex[sji + 1];
		for (int ind = start; ind < finish; ind++)
		{
			hash[l_mat->Col[ind]] = l_mat->Value[ind];
		}

		while ((i < f_ijk->k) && (f_ijk->coord2[i] == sji))
		{
			int ii = f_ijk->coord1[i];
			int ki = f_ijk->coord3[i];
			int ji = f_ijk->coord2[i];

			dcomplex v1 = hash[ii];
			dcomplex v2 = f_ijk->data[i];
			dcomplex tmp;
			tmp.re = v1.re * v2.re - v1.im * v2.im;
			tmp.im = v1.re * v2.im + v1.im * v2.re;

			Ks[ki].re += tmp.re;
			Ks[ki].im += tmp.im;
			i++;
		}
		i--;
	}

	dcomplex val;
	for (int i = 0; i < N_mat; i++)
	{
      Ks[i].re *= (1.0) / (N + 1);
      Ks[i].im *= (-1.0) / (N + 1);
		val = Ks[i];
		Ks[i].re = val.im;
		Ks[i].im = val.re;
	}
}
