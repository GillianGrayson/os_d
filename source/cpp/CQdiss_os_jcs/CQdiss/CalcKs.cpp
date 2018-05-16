#include "CalcKs.h"

void calcKs(Model *m)
{
	int N = m->N;
	int N_mat = m->N_mat;
	crsMatrix *Ks_tmp;
	crsMatrix *FsT;
	dcomplex  *Ks = m->Ks;
	crsMatrix *As = m->a_mat;
	crsMatrix **Fs = m->f_mat;
	crsMatrix *AsT;

	for(int i = 0; i < N_mat; i++)
	{
		Ks[i].re = 0.0;
		Ks[i].im = 0.0;
	}

	for(int i = 0; i < N_mat; i++)
	{
		AsT = As;

		FsT = Fs[i];
		FsT = new crsMatrix(*(Fs[i]));
		Transpose(*(Fs[i]), *FsT, false);

		for(int j = 0; j < N_mat; j++)
		{
			int ii,jj, m1, m2;
			ii = As->RowIndex[i];
			m1 = As->RowIndex[i + 1];
			jj = FsT->RowIndex[j];
			m2 = FsT->RowIndex[j + 1];

			while((ii < m1)&&(jj < m2))
			{
				if(AsT->Col[ii] < FsT->Col[jj])
				{  
					ii++;
				} 
				else
				{
					if(AsT->Col[ii] > FsT->Col[jj])
					{
						jj++;
					} 
					else
					{
						dcomplex as, fs;
						as = AsT->Value[ii];
						fs = FsT->Value[jj];
						Ks[j].re += as.re * fs.re - as.im * fs.im;
						Ks[j].im += as.re * fs.im + as.im * fs.re;

						ii++;
						jj++;
					}
				}
			}

		}
		delete FsT;
	}

	dcomplex val;
	for(int i = 0; i < N_mat; i++)
	{
		Ks[i].re *= (-1.0) / (N + 1);
		Ks[i].im *= (1.0) / (N + 1);
		val = Ks[i];
		Ks[i].re = val.im;
		Ks[i].im = val.re;
	}
}