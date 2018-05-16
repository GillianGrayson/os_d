#include "Init_a1_a2.h"
#include <math.h>

crsMatrix * create_A1_diss1_matrix(Model * m, int diss_id, RunParam &rp, ConfigParam &cp)
{
	int N = m->N;

	crsMatrix * a_std = create_a_std_matrix(m);

	return a_std;
}

crsMatrix * create_A2_diss1_matrix(Model * m, int diss_id, RunParam &rp, ConfigParam &cp)
{
	int N = m->N;

	crsMatrix * mat;

	mat = new crsMatrix((N + 1), (N + 1));

	int * RowIndex = mat->RowIndex;
	int * Col = mat->Col;
	dcomplex * Value = mat->Value;

	for (int state_id_1 = 0; state_id_1 < N + 1; state_id_1++)
	{
		Col[state_id_1] = state_id_1;
		RowIndex[state_id_1] = state_id_1;
		Value[state_id_1].re = 0.0;
	}
	RowIndex[N + 1] = N + 1;

	return mat;
}

vector<pair<int, dcomplex> > * create_a_std_matrix(crsMatrix * a1_mat, crsMatrix * a2_mat, int N_mat)
{
	vector<pair<int, dcomplex> > * a_std;
	a_std = new vector<pair<int, dcomplex> >[N_mat];

	for (int i = 0; i < N_mat; i++)
	{
		for (int j = 0; j < N_mat; j++)
		{
			dcomplex v1, v2, v3, v4, res1, res2, res;
			int c1, c2, c3, c4;
			c1 = a1_mat->RowIndex[i + 1] - a1_mat->RowIndex[i];
			c2 = a1_mat->RowIndex[j + 1] - a1_mat->RowIndex[j];
			c3 = a2_mat->RowIndex[i + 1] - a2_mat->RowIndex[i];
			c4 = a2_mat->RowIndex[j + 1] - a2_mat->RowIndex[j];

			if ((c1 + c3) * (c2 + c4) > 0)
			{

				v1.re = a1_mat->Value[a1_mat->RowIndex[i]].re * c1;
				v3.re = a2_mat->Value[a2_mat->RowIndex[i]].re * c3;
				v2.re = a1_mat->Value[a1_mat->RowIndex[j]].re * c2;
				v4.re = a2_mat->Value[a2_mat->RowIndex[j]].re * c4;

				v1.im = a1_mat->Value[a1_mat->RowIndex[i]].im * c1;
				v3.im = a2_mat->Value[a2_mat->RowIndex[i]].im * c3;
				v2.im = a1_mat->Value[a1_mat->RowIndex[j]].im * c2;
				v4.im = a2_mat->Value[a2_mat->RowIndex[j]].im * c4;

				res1.re = v1.re + v3.im;
				res1.im = v1.im - v3.re;

				res2.re = v2.re + v4.im;
				res2.im = v2.im - v4.re;

				res2.im = -res2.im;

				res.re = res1.re * res2.re - res1.im * res2.im;
				res.im = res1.re * res2.im + res1.im * res2.re;
				a_std[i].push_back(make_pair(j, res));
			}
		}
	}

	return a_std;
}

void init_diss_1(Model * m, RunParam &rp, ConfigParam &cp, MainData &md)
{
	int N = m->N;
	FMatrixs *Fs = m->Fs;

	int N_mat = m->N_mat;

	crsMatrix * result_a_matrix = NULL;

	double time = omp_get_wtime();
	double init_time = time;

	crsMatrix * A1 = create_A1_diss1_matrix(m, 0, rp, cp);
	crsMatrix * A2 = create_A2_diss1_matrix(m, 0, rp, cp);

	crsMatrix * a1_mat = new crsMatrix(N_mat, N_mat);
	crsMatrix * a2_mat = new crsMatrix(N_mat, N_mat);

	cout << "Dissipation" << endl;

	int k = 0;

	crsMatrix * res;
	int cnt;
	a1_mat->RowIndex[0] = 0;
	for (int i = 0; i < N_mat; i++)
	{
		res = new crsMatrix;
		SparseMKLMult(*A1, *(Fs->F[i + 1]), *res);
		cnt = trace_struct(*res);
		if (cnt > 0)
		{
			a1_mat->Value[k] = trace(*res);
			a1_mat->Col[k] = i;
			k++;
		}
		a1_mat->RowIndex[i + 1] = k;
		delete res;
	}
	a1_mat->NZ = k;

	k = 0;
	a2_mat->RowIndex[0] = 0;
	for (int i = 0; i < N_mat; i++)
	{
		res = new crsMatrix;
		SparseMKLMult(*A2, *(Fs->F[i + 1]), *res);
		cnt = trace_struct(*res);
		if (cnt > 0)
		{
			a2_mat->Value[k] = trace(*res);
			a2_mat->Col[k] = i;
			k++;
		}
		a2_mat->RowIndex[i + 1] = k;
		delete res;
	}
	a2_mat->NZ = k;

	vector<pair<int, dcomplex> > * a_std = create_a_std_matrix(a1_mat, a2_mat, N_mat);

	delete a1_mat;
	delete a2_mat;

	delete A1;
	delete A2;

	crsMatrix * a1_i_a2_mat = NULL;
	a1_i_a2_mat = stdToCrs(a_std, N_mat);

	delete[] a_std;

	scalar_mult(a1_i_a2_mat, cp.g);
	result_a_matrix = new crsMatrix(*a1_i_a2_mat);

	m->a_mat = new crsMatrix(*a1_i_a2_mat);

	delete a1_i_a2_mat;

	time = omp_get_wtime() - init_time;
	cout << "time of a_" << "0 : " << time << endl << endl;
}
