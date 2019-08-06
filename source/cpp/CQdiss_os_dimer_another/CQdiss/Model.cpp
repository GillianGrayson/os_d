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

	model->h_base = new crsMatrix(model->N_mat, model->N_mat);
	model->h_drv = new crsMatrix(model->N_mat, model->N_mat);

	model->H_base = new crsMatrix(model->N_mat, model->N_mat);
	model->H_drv = new crsMatrix(model->N_mat, model->N_mat);

	model->Qs_base = NULL;
	model->Qs_drv = NULL;

	model->Ks = new dcomplex[model->N_mat];
	model->Rs = NULL;
	model->Gs = NULL;

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
	if (model->h_base != NULL)
	{
		delete model->h_base;
	}
	if (model->h_drv != NULL)
	{
		delete model->h_drv;
	}

	if (model->H_base != NULL)
	{
		delete model->H_base;
	}
	if (model->H_drv != NULL)
	{
		delete model->H_drv;
	}

	if (model->Qs_base != NULL)
	{
		delete model->Qs_base;
	}
	if (model->Qs_drv != NULL)
	{
		delete model->Qs_drv;
	}

	if (model->Ks != NULL)
	{
		delete[] model->Ks;
	}
	if (model->Rs != NULL)
	{
		delete model->Rs;
	}
	if (model->Gs != NULL)
	{
		delete model->Gs;
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
