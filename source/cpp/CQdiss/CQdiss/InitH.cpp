#include "initH.h"
#include <math.h>
#include <algorithm>

crsMatrix * createHmatrix(Model * m, bool sub_trace)
{
	int N = m->N;
	double E0 = m->conf.E0;
	double U = m->conf.U;
	double J = m->conf.J;

	crsMatrix * H = new crsMatrix(N + 1, (N + 1) * 3 - 2);
	int * RowIndex = H->RowIndex;
	int * Col = H->Col;
	dcomplex * Value = H->Value;

	int i, j;
	RowIndex[0] = 0;
	RowIndex[1] = 2;
	for (i = 2; i < N + 1; i++)
	{
		RowIndex[i] = RowIndex[i - 1] + 3;
	}
	RowIndex[i] = RowIndex[i - 1] + 2;

	Col[0] = 0;
	Col[1] = 1;

	for (i = 1; i < N; i++)
	{
		j = RowIndex[i];
		Col[j + 0] = i - 1;
		Col[j + 1] = i;
		Col[j + 2] = i + 1;
	}
	j = RowIndex[i];
	Col[j + 0] = i - 1;
	Col[j + 1] = i;

	double * h_diag1 = new double[N + 1];
	double * h_diag2 = new double[N + 1];
	double tr = 0.0;

	for (i = 0; i < N + 1; i++)
	{
		h_diag1[i] = i * 2 - N;
	}

	for (i = 0; i < N + 1; i++)
	{
		h_diag2[i] = (h_diag1[i] * E0 + h_diag1[i] * h_diag1[i] * U / N);
		tr += h_diag2[i];
	}
	tr /= (N + 1);
	if (sub_trace)
	{
		for (i = 0; i < N + 1; i++)
		{
			h_diag2[i] -= tr;
		}
	}

	Value[0].re = h_diag2[0];
	Value[1].re = sqrt((double)N) * J;

	for (i = 1; i < N; i++)
	{
		j = RowIndex[i];
		Value[j + 0].re = J * sqrt((double)((N - i + 1)  * (i + 0)));
		Value[j + 1].re = h_diag2[i];
		Value[j + 2].re = J * sqrt((double)((N - i + 0)  * (i + 1)));
	}
	j = RowIndex[i];
	Value[j + 0].re = J * sqrt((double)N);
	Value[j + 1].re = h_diag2[i];

	//printMatrixVal(H);

	delete[] h_diag1;
	delete[] h_diag2;

	return H;
}

void to_F_basis(crsMatrix * Mat, crsMatrix * vec)
{
	int N = Mat->N;
	int N_mat = N * N - 1;
	int cnt = 0;
	int *mask = new int[N_mat];
	int *col = new int[N_mat];
	dcomplex *value = new dcomplex[N_mat];
	for (int i = 0; i < N_mat; i++)
	{
		mask[i] = 0;
		value[i].re = 0.0;
		value[i].im = 0.0;
	}

	dcomplex sum;
	sum.re = 0.0;
	sum.im = 0.0;
	for (int i = 0; i < N; i++)
	{
		int start = Mat->RowIndex[i];
		int finish = Mat->RowIndex[i + 1];
		for (int j = start; j < finish; j++)
		{
			int k = Mat->Col[j];

			if (k != i)
			{
				int ii = i;
				int jj = k;
				int z = -1;
				if (ii > jj) {
					ii = k;
					jj = i;
					z = 1;
				}

				int index = ((N - 1 + N - ii) * ii) + jj - ii - 1;
				if (mask[index] != 1)
				{
					col[cnt + 0] = index + 0;
					col[cnt + 1] = index + 1;
					cnt += 2;
				}

				mask[index] = 1;
				mask[index + 1] = 1;

				value[index].re += Mat->Value[j].re / sqrt(2.0);
				value[index].im += Mat->Value[j].im / sqrt(2.0);

				value[index + 1].re += z * Mat->Value[j].im / sqrt(2.0);
				value[index + 1].im += -z * Mat->Value[j].re / sqrt(2.0);

			}
			else
			{
				if (k != 0)
				{
					int index = N * (N - 1) + k - 1;
					mask[index] = 1;
					double value_d = 1.0 / sqrt((double)((k + 0)* (k + 1)));
					value[index].re = sum.re * value_d;
					value[index].re -= Mat->Value[j].re * (k + 0) * value_d;
					value[index].im = sum.im * value_d;
					value[index].im -= Mat->Value[j].im * (k + 0) * value_d;
					col[cnt] = index;
					cnt++;
				}

				sum.re += Mat->Value[j].re;
				sum.im += Mat->Value[j].im;
			}
		}
	}

	for (int i = 0; i < cnt; i++)
	{
		if ((value[col[i]].re == 0.0) && (value[col[i]].im == 0.0))
		{
			cnt--;
			col[i] = col[cnt];
			i--;
		}
	}
	std::sort(col, col + cnt);

	vec->NZ = cnt;
	for (int i = 0; i < vec->N + 1; i++)
	{
		vec->RowIndex[i] = 0;
	}
	for (int i = 0; i < cnt; i++)
	{
		vec->Col[i] = col[i];
		vec->Value[i] = value[col[i]];
		vec->RowIndex[col[i] + 1] ++;
	}
	for (int i = 0; i < vec->N; i++)
	{
		vec->RowIndex[i + 1] = vec->RowIndex[i] + vec->RowIndex[i + 1];
	}

	delete[] mask;
	delete[] col;
	delete[] value;
}

void init_h_vector_opt(Model * m)
{
	crsMatrix * H = createHmatrix(m), *res;
	int N = m->N;
	crsMatrix * h = m->h;

	to_F_basis(H, h);

	delete H;
}


void init_h_vector(Model * m)
{
	crsMatrix * H = createHmatrix(m), *res;
	int N = m->N;
	FMatrixs *Fs = m->Fs;
	crsMatrix * h = m->h;

	int cnt = 0;
	int k = 0;
	int i = 0;
	for (i = 0; i < (N + 1) * (N + 1) - 1; i++)
	{
		h->RowIndex[i] = k;
		res = new crsMatrix;
		SparseMKLMult(*H, *(Fs->F[i + 1]), *res);
		cnt = trace_struct(*res);
		if (cnt > 0)
		{
			h->Col[k] = i;
			h->Value[k] = trace(*res);
			if ((h->Value[k].re != 0.0) || (h->Value[k].im != 0.0))
			{
				k++;
			}
		}
		delete res;
	}
	h->RowIndex[i] = k;
	h->NZ = k;

	//  printMatrixVal(h);

	delete H;
}

crsMatrix * createHe_matrix(Model * m)
{
	int N = m->N;

	crsMatrix * H = new crsMatrix(N + 1, N + 1);
	int * RowIndex = H->RowIndex;
	int * Col = H->Col;
	dcomplex * Value = H->Value;

	int i, j;
	for (i = 0; i < N + 1; i++)
	{
		RowIndex[i] = i;
		Col[i] = i;
		Value[i].re = i * 2 - N;
	}
	RowIndex[i] = i;

	return H;
}

void init_he_vector_opt(Model * m)
{
	crsMatrix * H = createHe_matrix(m), *res;
	int N = m->N;
	FMatrixs *Fs = m->Fs;
	crsMatrix * he = m->he;

	to_F_basis(H, he);

	delete H;
}
void init_he_vector(Model * m)
{
	crsMatrix * H = createHe_matrix(m), *res;
	int N = m->N;
	FMatrixs *Fs = m->Fs;
	crsMatrix * he = m->he;

	int cnt = 0;
	int k = 0;
	int i = 0;
	for (i = 0; i < (N + 1) * (N + 1) - 1; i++)
	{
		he->RowIndex[i] = k;
		res = new crsMatrix;
		SparseMKLMult(*H, *(Fs->F[i + 1]), *res);
		he->Col[k] = i;
		he->Value[k] = trace(*res);
		if ((he->Value[k].re != 0.0) || (he->Value[k].im != 0.0))
		{
			k++;
		}
		delete res;
	}
	he->RowIndex[i] = k;
	he->NZ = k;


	delete H;
}

void print_h_vector(Model * m)
{
	int N = m->N;
	int i;
	for (i = 0; i < (N + 1) * (N + 1) - 1; i++)
	{
		//   printf("(%5.2lf;%5.2lf)", m->h[i].re, m->h[i].im);
	}
}
