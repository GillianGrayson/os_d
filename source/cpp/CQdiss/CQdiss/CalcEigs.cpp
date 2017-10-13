#include "CalcGs.h"
#include <string.h>
#include <mkl.h>
#include <omp.h>
#include <stdio.h>
#include <math.h> 

void calcEig(Model *m)
{
	crsMatrix * Gs = m->Gs;
	int N_mat = Gs->N;
	dcomplex * fillMat = new dcomplex[N_mat * N_mat];
	dcomplex * DP = new dcomplex[N_mat];
	int i, j, k, s, f;

	memset(fillMat, 0, sizeof(dcomplex) * N_mat * N_mat);
	for (i = 0; i < N_mat; i++)
	{
		s = Gs->RowIndex[i];
		f = Gs->RowIndex[i + 1];
		for (k = s; k < f; k++)
		{
			j = Gs->Col[k];
			fillMat[i * N_mat + j].re = Gs->Value[k].re;
			fillMat[i * N_mat + j].im = Gs->Value[k].im;
		}
	}

	double l_time = omp_get_wtime();
	int info;
	info = LAPACKE_zgeev(LAPACK_ROW_MAJOR, /*'V'*/'N', /*'V'*/'N', N_mat, (MKL_Complex16 *)fillMat, N_mat,
		(MKL_Complex16 *)DP, NULL/* (MKL_Complex16 *)&LeftSigmaP[0][0]*/,
		N_mat, NULL/*(MKL_Complex16 *)&SigmaP[0][0]*/, N_mat);
	/* Check for convergence */
	if (info > 0) {
		printf("The algorithm failed to compute eigenvalues.\n");
		exit(1);
	}
	l_time = omp_get_wtime() - l_time;
	printf("LAPACKE_zgeev time:%lf\n", l_time);

	FILE * file = fopen("eigs.txt", "w");

	for (i = 0; i < N_mat; i++)
	{
		fprintf(file, "%1.16lf %1.16lf\n", DP[i].re, DP[i].im);
	}

	fclose(file);

	delete[] DP;
	delete[] fillMat;
}


void check_rho_evals(Model * m)
{
	crsMatrix * Rho = m->Rho;
	int N = m->N;
	dcomplex * rho_dense = new dcomplex[(N+1) * (N+1)];
	dcomplex * evals = new dcomplex[(N + 1)];

	for (int st_id_1 = 0; st_id_1 < (N + 1); st_id_1++)
	{
		for (int st_id_2 = 0; st_id_2 < (N + 1); st_id_2++)
		{
			rho_dense[st_id_1 * (N + 1) + st_id_2].re = 0.0;
			rho_dense[st_id_1 * (N + 1) + st_id_2].im = 0.0;
		}

		evals[st_id_1].re = 0.0;
		evals[st_id_1].im = 0.0;
	}

	int s = 0;
	int f = 0;
	int k = 0;
	int j = 0;
	for (int i = 0; i < (N + 1); i++)
	{
		s = Rho->RowIndex[i];
		f = Rho->RowIndex[i + 1];
		for (k = s; k < f; k++)
		{
			j = Rho->Col[k];
			rho_dense[i * (N+1) + j].re = Rho->Value[k].re;
			rho_dense[i * (N+1) + j].im = Rho->Value[k].im;
		}
	}

	int info;
	info = LAPACKE_zgeev(
		LAPACK_ROW_MAJOR,
		'N',
		'N',
		(N + 1),
		(MKL_Complex16 *)rho_dense,
		(N + 1),
		(MKL_Complex16 *)evals,
		NULL,
		(N + 1),
		NULL,
		(N + 1)
	);

	/* Check for convergence */
	if (info > 0) {
		printf("The algorithm failed to compute eigenvalues.\n");
		exit(1);
	}


	double sum = 0.0;
	double eps_imag = 1.0e-14;
	double eps_norm = 1.0e-14;
	int is_evals_ok = 1;
	for (int eval_id = 0; eval_id < (N + 1); eval_id++)
	{
		if (fabs(evals[eval_id].im) > eps_imag)
		{
			printf("Imag part of eval too big. Index: %d\n", eval_id);
			is_evals_ok = 0;
		}

		if ((evals[eval_id].re < 0.0) && (evals[eval_id].re > 1.0))
		{
			printf("Wrong eval value. Index: %d\n", eval_id);
			is_evals_ok = 0;
		}

		sum += evals[eval_id].re;
	}

	if (fabs(sum - 1.0) > eps_norm)
	{
		printf("Sum of evals is not 1.0. Diff: %0.16le \n", sum - 1.0);
		is_evals_ok = 0;
	}

	if (is_evals_ok == 1)
	{
		printf("Evals is ok. Diff: %0.16le \n", sum - 1.0);
	}

	delete[] rho_dense;
	delete[] evals;
}