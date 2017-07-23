#include "CalcRho.h"
#include "math.h"

void calcRho(Model *m)
{
  FMatrixs  * Fs = m->Fs;
  crsMatrix * Rho = new crsMatrix(*(Fs->F[0]));
  crsMatrix * Rho_tmp = new crsMatrix(*(Fs->F[0]));
  crsMatrix * tmp;
  dcomplex  * RhoF = m->RhoF;
  int N_mat = m->N_mat;
  int N = m->N;
 
  int i;
  for(i = 0; i < Rho->NZ; i++)
  {
    Rho->Value[i].re /= N + 1;
    Rho->Value[i].im /= N + 1;
  }

  for(int i = 0; i < N_mat; i++)
  {
    SparseMKLAdd(*Rho, RhoF[i], *(Fs->F[i + 1]), *Rho_tmp, true); 
    tmp = Rho; Rho = Rho_tmp; Rho_tmp = tmp;
  }
  delete Rho_tmp;
  
  if (m->Rho != NULL) delete m->Rho;
  m->Rho = Rho;
}

void calcRho_fill(Model *m)
{
	int N_mat = m->N_mat;
	int N = m->N;
	dcomplex  * RhoF = m->RhoF;
	crsMatrix * mRho = new crsMatrix(N + 1, (N + 1) * (N + 1));
	dcomplex  * Rho = mRho->Value;

	int i, j;
	for (i = 0; i < (N + 1); i++)
	{
		mRho->RowIndex[i] = i * (N + 1);
		for (j = 0; j < (N + 1); j++)
		{
			mRho->Col[i * (N + 1) + j] = j;
		}
	}
	mRho->RowIndex[i] = i * (N + 1);
	
	
	for (i = 0; i < (N + 1) * (N + 1); i++)
	{
		Rho[i].re = 0.0;
		Rho[i].im = 0.0;
	}
	for (i = 0; i < (N + 1); i++)
	{
		Rho[i * (N + 1) + i].re = 1.0 / (N + 1);
	}

	int k = 0;
	double val = 1.0 / sqrt(2.0);
	for (i = 0; i < N + 1; i++)
	{
		for (j = i + 1; j < N + 1; j++)
		{
			Rho[(i) * (N + 1) + (j)].re += RhoF[k].re * val;
			Rho[(j) * (N + 1) + (i)].re += RhoF[k].re * val;
			k++;

			Rho[(i) * (N + 1) + (j)].im += RhoF[k].re * (-val);
			Rho[(j) * (N + 1) + (i)].im += RhoF[k].re * (+val);
			k++;
		}
	}

	for (i = 0; i < N; i++)
	{
		val = 1.0 / sqrt((double)((i + 1)* (i + 2)));
		for (j = 0; j <= i; j++)
		{
			Rho[j * (N + 1) + j].re += RhoF[k].re * val;
		}
		Rho[j * (N + 1) + j].re -= RhoF[k].re * val * (i + 1);
		k++;
	}



	if (m->Rho != NULL) delete m->Rho;
	m->Rho = mRho;
}

