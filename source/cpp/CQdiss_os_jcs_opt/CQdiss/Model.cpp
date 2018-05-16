#include "Model.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

Model * createModel(int N, ConfigParam conf)
{
	int i;
	Model * model = new Model;

	model->N = N;
	model->N_mat = (N + 1) * (N + 1) - 1;
	model->conf = conf;

	model->h_0 = new crsMatrix(model->N_mat, model->N_mat);
	model->h_1 = new crsMatrix(model->N_mat, model->N_mat);

	model->H_0 = new crsMatrix(model->N_mat, model->N_mat);
	model->H_1 = new crsMatrix(model->N_mat, model->N_mat);

	model->H0 = NULL;
	model->H1 = NULL;

	model->Q_0 = NULL;
	model->Q_1 = NULL;

	model->Ks = new dcomplex[model->N_mat];
	model->Rs = NULL;
	model->G_0_s = NULL;
	model->G_1_s = NULL;

	model->prevRhoF = new dcomplex[model->N_mat];
	model->RhoF = new dcomplex[model->N_mat];
	memset(model->RhoF, 0, sizeof(dcomplex) * model->N_mat);
	memset(model->prevRhoF, 0, sizeof(dcomplex) * model->N_mat);
	for (i = 0; i < model->N_mat; i++)
	{
		model->RhoF[i].re = (double)rand() / (double)RAND_MAX;
		model->RhoF[i].im = 0.0;
	}

	model->Rho = NULL;

	model->f_ijk = NULL;

	model->l_mat = NULL;

	return model;
}
void freeModel(Model * model)
{
	if (model->h_0 != NULL)
	{
		delete model->h_0;
	}
	if (model->h_1 != NULL)
	{
		delete model->h_1;
	}

	if (model->H_0 != NULL)
	{
		delete model->H_0;
	}
	if (model->H_1 != NULL)
	{
		delete model->H_1;
	}

	if (model->H0 != NULL)
	{
		delete model->H0;
	}
	if (model->H1 != NULL)
	{
		delete model->H1;
	}

	if (model->Q_0 != NULL)
	{
		delete model->Q_0;
	}
	if (model->Q_1 != NULL)
	{
		delete model->Q_1;
	}

	if (model->Ks != NULL)
	{
		delete[] model->Ks;
	}
	if (model->Rs != NULL)
	{
		delete model->Rs;
	}
	if (model->G_0_s != NULL)
	{
		delete model->G_0_s;
	}
	if (model->G_1_s != NULL)
	{
		delete model->G_1_s;
	}

	if (model->RhoF != NULL)
	{
		delete[] model->RhoF;
	}

	if (model->prevRhoF != NULL)
	{
		delete[] model->prevRhoF;
	}

	if (model->Rho != NULL)
	{
		delete model->Rho;
	}

	if (model->l_mat != NULL)
	{
		delete model->l_mat;
	}

	if (model->f_ijk != NULL)
	{
		free_matrix(model->f_ijk);
		delete model->f_ijk;
	}

	delete model;
}

crsMatrix * create_a_std_matrix(Model * m)
{
	int N = m->N;

	crsMatrix * a_mtx = new crsMatrix(N + 1, N);

	for (int i = 0; i < N; i++)
	{
		a_mtx->RowIndex[i] = i;
		a_mtx->Col[i] = i + 1;
		a_mtx->Value[i].re = sqrt(double(i + 1));
	}
	a_mtx->RowIndex[N] = N;
	a_mtx->RowIndex[N + 1] = N;

	return a_mtx;
}

crsMatrix * create_a_dag_matrix(Model * m)
{
	int N = m->N;

	crsMatrix * a_mtx = new crsMatrix(N + 1, N);

	a_mtx->RowIndex[0] = 0;
	for (int i = 0; i < N; i++)
	{
		a_mtx->RowIndex[i + 1] = i;
		a_mtx->Col[i] = i;
		a_mtx->Value[i].re = sqrt(double(i + 1));
	}
	a_mtx->RowIndex[N + 1] = N;

	return a_mtx;
}