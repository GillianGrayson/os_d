#include "initH.h"
#include <math.h>
#include <mkl.h>

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
	FMatrixs *Fs = m->Fs;
	dcomplex * h_0 = m->h_0;

	int i = 0;
	for (i = 0; i < (N + 1) * (N + 1) - 1; i++)
	{
		res = new crsMatrix;
		SparseMKLMult(*H_0, *(Fs->F[i + 1]), *res);
		h_0[i] = trace(*res);

		delete res;
	}
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
	FMatrixs *Fs = m->Fs;
	dcomplex * h_1 = m->h_1;

	int i = 0;
	for (i = 0; i < (N + 1) * (N + 1) - 1; i++)
	{
		res = new crsMatrix;
		SparseMKLMult(*H_1, *(Fs->F[i + 1]), *res);
		h_1[i] = trace(*res);

		delete res;
	}
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