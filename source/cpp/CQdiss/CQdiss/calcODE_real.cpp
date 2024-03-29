#include "calcODE.h"
#include "genMatrix.h"
#include <math.h>
#include <stdlib.h>
#include <mkl.h>
#include "calcRho.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>

using namespace std;

void save_double_data(string file_name, double * data, int size, int precision, bool append)
{
	if (append)
	{
		ofstream ofs = ofstream(file_name, ios::app);

		if (ofs.is_open())
		{
			ofs << setprecision(precision) << scientific;
			for (int i = 0; i < size; i++)
			{
				ofs << data[i] << endl;
			}

			ofs.close();
		}
	}
	else
	{
		ofstream ofs = ofstream(file_name);

		if (ofs.is_open())
		{
			ofs << setprecision(precision) << scientific;
			for (int i = 0; i < size; i++)
			{
				ofs << data[i] << endl;
			}

			ofs.close();
		}
	}
}

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

	int driving_type = m->conf.driving_type;

	multMatVec_real(Gs, x, tmp1);
	multMatVec_real(QEs, x, tmp2);

	for (i = 0; i < N_mat; i++)
	{
		if (driving_type == 1)
		{
			res[i] = (tmp1[i] + A0 * sin(w * t) * tmp2[i] - Ks[i]) * h;
		}
		else if (driving_type == 0)
		{
			if (sin(w * t) < 0.0)
			{
				res[i] = (tmp1[i] + A0 * (-1.0) * tmp2[i] - Ks[i]) * h;
			}
			else
			{
				res[i] = (tmp1[i] + A0 * (+1.0) * tmp2[i] - Ks[i]) * h;
			}
		}
	}
}

void calcVectValue_mult(double t, double h,
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

	int driving_type = m->conf.driving_type;

	multMatVec_real(Gs, x, tmp1);
	multMatVec_real(QEs, x, tmp2);

	for (i = 0; i < N_mat; i++)
	{
		if (driving_type == 1)
		{
			res[i] = (tmp1[i] + A0 * sin(w * t) * tmp2[i]) * h;
		}
		else if (driving_type == 0)
		{
			if (sin(w * t) < 0.0)
			{
				res[i] = (tmp1[i] + A0 * (-1.0) * tmp2[i]) * h;
			}
			else
			{
				res[i] = (tmp1[i] + A0 * (+1.0) * tmp2[i]) * h;
			}
		}
	}
}

void calcODE_real(Model *m, double h, int cntItr, double t)
{
	double time;
	int itr, i;
	int N_mat = m->N_mat;
	double * RhoF = (double *)(m->RhoF);

	for (i = 0; i < N_mat; i++)
	{
		m->prevRhoF[i].re = m->RhoF[i].re;
		m->prevRhoF[i].im = m->RhoF[i].im;
	}

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


		if (m->conf.deep_dump == 2)
		{
			if (itr % 100 == 0)
			{
				real_to_complex(m->RhoF, m->N_mat);

				calcRho_fill(m);

				string fn = "diag_rho_deep_evo.txt";
				ofstream ofs = ofstream(fn, ios::app);

				if (ofs.is_open())
				{
					ofs << setprecision(16) << scientific;

					for (int i = 0; i < m->Rho->N; i++)
					{
						for (int k = m->Rho->RowIndex[i]; k < m->Rho->RowIndex[i + 1]; k++)
						{
							int j = m->Rho->Col[k];
							if (i == j)
							{
								ofs << m->Rho->Value[k].re << endl;
							}
						}
					}

					ofs.close();
				}

				complex_to_real(m->RhoF, m->N_mat);
			}
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


void calcODE_mult(Model *m, double h, int cntItr, double t)
{
	double time;
	int itr, i;
	int N_mat = m->N_mat;
	double * RhoF = (double *)(m->RhoF);

	//for (i = 0; i < N_mat; i++)
	//{
	//	m->prevRhoF[i].re = m->RhoF[i].re;
	//	m->prevRhoF[i].im = m->RhoF[i].im;
	//}

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
