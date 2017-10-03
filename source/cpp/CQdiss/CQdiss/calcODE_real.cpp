#include "calcODE.h"
#include "genMatrix.h"
#include <math.h>
#include <stdlib.h>
#include <mkl.h>

void complex_to_real(dcomplex *mat, int N)
{
	double *value = (double *)(mat);
	for (int i = 0; i < N; i++)
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

void multMatVec_real(crsMatrix *mat, double * x, double * res)
{

	char trans = 'n';

	double *value = (double *)(mat->Value);

	mkl_dcsrgemv(&trans, &(mat->N), value, mat->RowIndex, mat->Col, x, res);
}

void calcVectValue_real(double t, double h,
	Model * m, double *x, double * res,
	double * tmp1, double * tmp2)
{
	int i;
	int N_mat = m->N_mat;
	crsMatrix * Gs = m->Gs;
	crsMatrix * QEs = m->QEs;
	double    * Ks = (double *)(m->Ks);
	double T = m->conf.T;

	double A0 = m->conf.A0;
	double w = m->conf.w;

	multMatVec_real(Gs, x, tmp1);
	multMatVec_real(QEs, x, tmp2);

	for (i = 0; i < N_mat; i++)
	{
		//res[i] = (tmp1[i] + A0 * sin(w * t) * tmp2[i] - Ks[i]) * h;
		if(sin(w * t) < 0.0)
		{
		  res[i] = (tmp1[i] + A0 * (-1.0) * tmp2[i] - Ks[i]) * h;
		}
		else
		{
		  res[i] = (tmp1[i] + A0 * (+1.0) * tmp2[i] - Ks[i]) * h;
		}
	}
}

void calcODE_real(Model *m, double h, int cntItr, double t)
{
	double time;
	int itr, i;
	int N_mat = m->N_mat;
	double * RhoF = (double *)(m->RhoF);

	double * k1 = new double[N_mat];
	double * k2 = new double[N_mat];
	double * k3 = new double[N_mat];
	double * k4 = new double[N_mat];
	double * val = new double[N_mat];
	double * tmp1 = new double[N_mat];
	double * tmp2 = new double[N_mat];
	double * tmp3 = new double[N_mat];

	for (itr = 0; itr < cntItr; itr++)
	{
		time = t + h * itr;
		calcVectValue_real(time, h, m, RhoF, k1, tmp1, tmp2);
		for (i = 0; i < N_mat; i++)
		{
			val[i] = RhoF[i] + k1[i] / 2.0;
		}
		calcVectValue_real(time + h / 2.0, h, m, val, k2, tmp1, tmp2);
		for (i = 0; i < N_mat; i++)
		{
			val[i] = RhoF[i] + k2[i] / 2.0;
		}
		calcVectValue_real(time + h / 2.0, h, m, val, k3, tmp1, tmp2);
		for (i = 0; i < N_mat; i++)
		{
			val[i] = RhoF[i] + k3[i];
		}
		calcVectValue_real(time + h, h, m, val, k4, tmp1, tmp2);

		for (i = 0; i < N_mat; i++)
		{
			RhoF[i] += (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
		}
	}

	delete[] k1;
	delete[] k2;
	delete[] k3;
	delete[] k4;
	delete[] val;
	delete[] tmp1;
	delete[] tmp2;
	delete[] tmp3;
}

