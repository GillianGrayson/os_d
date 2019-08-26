#include "Model.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//создание модели
Model * createModel(int N, ConfigParam conf)
{
  int i;
  Model * model = new Model;

  model -> memlog = fopen("memlog.txt", "w");

  model->N = N;
  model->N_mat = (N+1) * (N+1) - 1;
  model->conf = conf;

  //model->Fs = new FMatrixs;
  //createFMatrixs(model->Fs, N);

  model->h = new crsMatrix(model->N_mat, model->N_mat);
  model->he = new crsMatrix(model->N_mat, model->N_mat);

//  model->f_mat = NULL;
//  model->f_H_mat = NULL;
//  model->d_mat = NULL;

//  model->a_mat = NULL;

  model->Qs = NULL;
  model->QEs = NULL;
  model->Ks = new dcomplex[model->N_mat];
  model->Rs = NULL;
  model->Gs = NULL;

  model->prevRhoF = new dcomplex[model->N_mat];
  model->RhoF = new dcomplex[model->N_mat];
  memset(model->RhoF, 0, sizeof(dcomplex) * model->N_mat);
  memset(model->prevRhoF, 0, sizeof(dcomplex) * model->N_mat);
  for(i = 0; i < model->N_mat; i++)
  {
    model->RhoF[i].re = (double)rand() / (double)RAND_MAX;
    model->RhoF[i].im = 0.0;
  }
  
  model->f_ijk = NULL;
//  model->d_ijk = NULL;
  model->Rho = NULL;

  model->l_mat = NULL;

  return model;
}

void createFMatrixs(FMatrixs * Fs, int N)
{
  Fs->countF = (2 + N) * N + 1;
  Fs->F = new crsMatrix *[Fs->countF];
  for(int i = 0; i < Fs->countF; i++)
  {
    Fs->F[i] = NULL;
  }
}

//освобождение памяти из под модели
void freeModel(Model * model)
{
  fclose(model->memlog);
  //freeFMatrixs(model->Fs);
  //delete model->Fs;
  //delete[] model->h;
  delete model->h;
  delete model->he;

//  if(model->f_mat != NULL)
//  {
//    int N = model->N;
//    for(int i = 0; i < (N + 1) * (N + 1) - 1; i++)
//    {
//      delete model->f_mat[i];
//    }
//    delete[] model->f_mat;
//  }

//  if(model->f_H_mat != NULL)
//  {
//    int N = model->N;
//    for(int i = 0; i < (N + 1) * (N + 1) - 1; i++)
//    {
//      delete model->f_H_mat[i];
//    }
//    delete[] model->f_H_mat;
//  }

//  if(model->d_mat != NULL)
//  {
//    int N = model->N;
//    for(int i = 0; i < (N + 1) * (N + 1) - 1; i++)
//    {
//      delete model->d_mat[i];
//    }
//    delete[] model->d_mat;
//  }

//  if(model->a_mat != NULL)
//  {
//    delete model->a_mat;
//  }

  if(model->Qs != NULL)
  {
    delete model->Qs;
  }

  if(model->QEs != NULL)
  {
    delete model->QEs;
  }


  if(model->Ks != NULL)
  {
    delete []model->Ks;
  }
  if(model->Rs != NULL)
  {
    delete model->Rs;
  }

  if(model->Gs != NULL)
  {
    delete model->Gs;
  }

  if(model->RhoF != NULL)
  {
    delete[] model->RhoF;
  }

  if(model->prevRhoF != NULL)
  {
    delete[] model->prevRhoF;
  }  
  
  if(model->Rho != NULL)
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

void freeFMatrixs(FMatrixs * Fs)
{
  for(int i = 0; i < Fs->countF; i++)
  {
    if(Fs->F[i] != NULL)
    {
      delete Fs->F[i];
    }
  }
  delete[] Fs->F;
}