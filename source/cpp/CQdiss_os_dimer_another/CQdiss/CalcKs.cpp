#include "CalcKs.h"
#include "coef_coord.h"
#include "sortTensor.h"
#include <string.h>
#include "f_d_sym.h"
#include "f_d_ch.h"

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
		Ks[i].re *= (-m->conf.diss_gamma) / (N + 1);
		Ks[i].im *= (m->conf.diss_gamma) / (N + 1);
		val = Ks[i];
		Ks[i].re = val.im;
		Ks[i].im = val.re;
	}
}


void calcKs_dimer(Model *m)
{
	int N = m->N;
	int N_mat = m->N_mat;
	crsMatrix *Ks_tmp;
	crsMatrix *FsT;
	dcomplex  *Ks = m->Ks;
	crsMatrix *AsT;

	crsMatrix * l_mat = m->l_mat;

	ulli cnt = fijk_coord_sym(m->l_mat, m->N + 1);
	m->f_ijk = new Tensor_Coordinates;
	fijk_coord_ch(m->f_ijk, m->l_mat, cnt + 1, m->N + 1);
	Tensor_Coordinates * f_ijk = m->f_ijk;

	for (int i = 0; i < N_mat; i++)
	{
		Ks[i].re = 0.0;
		Ks[i].im = 0.0;
	}

	for (ulli i = 0; i < f_ijk->k; i++)
	{
		int ii = f_ijk->coord1[i];
		int ji = f_ijk->coord2[i];
		int ki = f_ijk->coord3[i];
		int c1 = l_mat->RowIndex[ii + 1] - l_mat->RowIndex[ii];
		int c2 = l_mat->RowIndex[ji + 1] - l_mat->RowIndex[ji];
		if (c1 * c2 > 0)
		{
			dcomplex v1 = l_mat->Value[l_mat->RowIndex[ii]];
			dcomplex v2 = l_mat->Value[l_mat->RowIndex[ji]];
			v2.im = -v2.im;
			dcomplex v3 = f_ijk->data[i];
			dcomplex tmp;
			tmp.re = v1.re * v2.re - v1.im * v2.im;
			tmp.im = v1.re * v2.im + v1.im * v2.re;

			Ks[ki].re += tmp.re * v3.re - tmp.im * v3.im;
			Ks[ki].im += tmp.re * v3.im + tmp.im * v3.re;
		}
	}

	dcomplex val;
	for (int i = 0; i < N_mat; i++)
	{
		Ks[i].re *= (-m->conf.diss_gamma) / (N + 1);
		Ks[i].im *= (m->conf.diss_gamma) / (N + 1);
		val = Ks[i];
		Ks[i].re = val.im;
		Ks[i].im = val.re;
	}
}
