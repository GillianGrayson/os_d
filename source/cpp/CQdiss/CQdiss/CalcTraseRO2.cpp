#include "CalcTraseRO2.h"
#include "Matrix_op.h"
#include <stdio.h>
#include <math.h>

void calcTraseRO2(Model *m)
{
  crsMatrix * Rho = m->Rho;
  crsMatrix * Rho2 = new crsMatrix;
  
  int i, j;
  toOneBase(*Rho);

  SparseMKLMultOne(*Rho, *Rho, *Rho2);

  toZeroBase(*Rho);
  toZeroBase(*Rho2);
  
  dcomplex tr= trace(*Rho2);
  
  delete Rho2;
  
  FILE * f = fopen("purity_final.txt", "w");
  
  fprintf(f, "%0.16le %0.16le\n", tr.re, tr.im);
  
  fclose(f);
}

void calc_negativity_final(Model *m)
{
	crsMatrix * Rho = m->Rho;

	double result = 0.0;

	for (int i = 0; i < Rho->N; i++)
	{
		for (int k = Rho->RowIndex[i]; k < Rho->RowIndex[i + 1]; k++)
		{
			int j = Rho->Col[k];
			if (i != j)
			{
				result += sqrt(Rho->Value[k].re * Rho->Value[k].re + Rho->Value[k].im * Rho->Value[k].im);
			}
		}
	}

	result *= 0.5;

	FILE * f = fopen("negativity_final.txt", "w");

	fprintf(f, "%0.16le\n", result);

	fclose(f);
}

double calc_purity(Model *m)
{
	crsMatrix * Rho = m->Rho;
	crsMatrix * Rho2 = new crsMatrix;

	int i, j;
	toOneBase(*Rho);

	SparseMKLMultOne(*Rho, *Rho, *Rho2);

	toZeroBase(*Rho);
	toZeroBase(*Rho2);

	dcomplex tr = trace(*Rho2);

	delete Rho2;

	return tr.re;
}

double calc_negativity(Model *m)
{
	crsMatrix * Rho = m->Rho;

	double result = 0.0;

	for (int i = 0; i < Rho->N; i++)
	{
		for (int k = Rho->RowIndex[i]; k < Rho->RowIndex[i + 1]; k++)
		{
			int j = Rho->Col[k];
			if (i != j)
			{
				result += sqrt(Rho->Value[k].re * Rho->Value[k].re + Rho->Value[k].im * Rho->Value[k].im);
			}
		}
	}

	result *= 0.5;

	return result;
}
