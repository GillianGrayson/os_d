#include "CalcODE_rk.h"
#include <math.h>
#include <string.h>
#include <mkl.h>

dcomplex * initRhoODE_rk(Model *m)
{
	int i;
	int N = m->N + 1;
	dcomplex * rho;

	rho = new dcomplex[N * N];
	for (int i = 0; i < N * N; i++)
	{
		rho[i].re = 0.0;
		rho[i].im = 0.0;
	}

	for (i = 0; i < N; i++)
	{
		rho[i * N + i].re = 1.0 / N;//sqrt(1.0 / N / 2.0);
		rho[i * N + i].im = 0.0;//sqrt(1.0 / N / 2.0);
	}

	return rho;
}

void calcVectValue_rk(double t, double h,
	Model * m, dcomplex * H, dcomplex * He, dcomplex * L,
	dcomplex *x, dcomplex * res,
	dcomplex * tmp1, dcomplex * tmp2, dcomplex * tmp3);

void CRStoMAT(crsMatrix * sm, dcomplex *m);

void calcODE_rk(Model *m, crsMatrix * H, crsMatrix * He, crsMatrix * L,
	dcomplex *rho, double h, int cntItr, double t)
{
	double time;
	int itr, i;
	int N = m->N+1;
	int Np2 = N * N;
	dcomplex  * RhoF = rho;

	dcomplex * H_mat = new dcomplex[Np2];
	dcomplex * L_mat = new dcomplex[Np2];
	dcomplex * He_mat = new dcomplex[Np2];
	
	dcomplex * k1 = new dcomplex[Np2];
	dcomplex * k2   = new dcomplex[Np2];
	dcomplex * k3   = new dcomplex[Np2];
	dcomplex * k4   = new dcomplex[Np2];
	dcomplex * val  = new dcomplex[Np2];
	dcomplex * tmp1 = new dcomplex[Np2];
	dcomplex * tmp2 = new dcomplex[Np2];
	dcomplex * tmp3 = new dcomplex[Np2];

	CRStoMAT(H, H_mat);
	CRStoMAT(L, L_mat);
	CRStoMAT(He, He_mat);

	for (itr = 0; itr < cntItr; itr++)
	{
		time = t + h * itr;
		calcVectValue_rk(time, h, m, H_mat, He_mat, L_mat, RhoF, k1, tmp1, tmp2, tmp3);
		for (i = 0; i < Np2; i++)
		{
			val[i].re = RhoF[i].re + k1[i].re / 2.0;
			val[i].im = RhoF[i].im + k1[i].im / 2.0;
		}
		calcVectValue_rk(time + h / 2.0, h, m, H_mat, He_mat, L_mat, val, k2, tmp1, tmp2, tmp3);
		for (i = 0; i < Np2; i++)
		{
			val[i].re = RhoF[i].re + k2[i].re / 2.0;
			val[i].im = RhoF[i].im + k2[i].im / 2.0;
		}
		calcVectValue_rk(time + h / 2.0, h, m, H_mat, He_mat, L_mat, val, k3, tmp1, tmp2, tmp3);
		for (i = 0; i < Np2; i++)
		{
			val[i].re = RhoF[i].re + k3[i].re;
			val[i].im = RhoF[i].im + k3[i].im;
		}
		calcVectValue_rk(time + h, h, m, H_mat, He_mat, L_mat, val, k4, tmp1, tmp2, tmp3);

		for (i = 0; i < Np2; i++)
		{
			RhoF[i].re += (k1[i].re + 2.0 * k2[i].re + 2.0 * k3[i].re + k4[i].re) / 6.0;
			RhoF[i].im += (k1[i].im + 2.0 * k2[i].im + 2.0 * k3[i].im + k4[i].im) / 6.0;
		}
	}

	delete[] H_mat;
	delete[] He_mat;
	delete[] L_mat;

	delete[] k1;
	delete[] k2;
	delete[] k3;
	delete[] k4;
	delete[] val;
	delete[] tmp1;
	delete[] tmp2;
	delete[] tmp3;
}

void CRStoMAT(crsMatrix * sm, dcomplex *m)
{
	int i, j, c;
	memset(m, 0.0, sizeof(dcomplex) * sm->N * sm->N);
	for (i = 0; i < sm->N; i++)
	{
		for (j = sm->RowIndex[i]; j < sm->RowIndex[i + 1]; j++)
		{
			c = sm->Col[j];
			m[i * sm->N + c] = sm->Value[j];
		}
	}
}

void calcVectValue_rk(double t, double h,
	Model * m, dcomplex * H, dcomplex * He, dcomplex * L,
	dcomplex *x, dcomplex * res,
	dcomplex * tmp1, dcomplex * tmp2, dcomplex * tmp3)
{
	int i;
	int N = m->N + 1;
	int Np2 = N * N;
	double g = m->conf.g / 2.0;

	dcomplex alpha; 
	dcomplex beta;
	alpha.re = 1.0;
	alpha.im = 0.0;
	beta.re = 0.0;
	beta.im = 0.0;

	double A0 = m->conf.A0;
	double w = m->conf.w;

	for (i = 0; i < Np2; i++)
	{
		//res[i].re = (tmp1[i].re + A0 * sin(w * t) * tmp2[i].re - Ks[i].re) * h;
		if (sin(w * t) < 0.0)
		{
			tmp1[i].re = (H[i].re + A0 * (-1.0) * He[i].re);
			tmp1[i].im = (H[i].im + A0 * (-1.0) * He[i].im);
		}
		else
		{
			tmp1[i].re = (H[i].re + A0 * (+1.0) * He[i].re);
			tmp1[i].im = (H[i].im + A0 * (+1.0) * He[i].im);
		}
		//res[i].im = (tmp1[i].im + A0 * cos(w * t) * tmp2[i].im - Ks[i].im) * h;
		//res[i].re = (tmp1[i].re + A0 * sin(w * t) * tmp2[i].re - Ks[i].re) * h;
	}

	//[H, rho]
	cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
		N, N, N, &alpha, tmp1, N, x, N, &beta, tmp2, N);
	cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
		N, N, N, &alpha, x, N, tmp1, N, &beta, tmp3, N);
	//-i [H, rho] * h
	for (i = 0; i < Np2; i++)
	{
		res[i].re = +(tmp2[i].im - tmp3[i].im) * h;
		res[i].im = -(tmp2[i].re - tmp3[i].re) * h;
	}
//	printMatrixVal_com(x, N);
//	printMatrixVal_com(tmp1, N);
//	printMatrixVal_com(tmp2, N);
//	printMatrixVal_com(tmp3, N);
//	printMatrixVal_com(res, N);

	// rho L_c
	cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans,
		N, N, N, &alpha, x, N, L, N, &beta, tmp1, N);
	//[L, rho L_c]
	cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
		N, N, N, &alpha, L, N, tmp1, N, &beta, tmp2, N);
	cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
		N, N, N, &alpha, tmp1, N, L, N, &beta, tmp3, N);
//	printMatrixVal_com(x, N);
//	printMatrixVal_com(L, N);
//	printMatrixVal_com(tmp1, N);
//	printMatrixVal_com(tmp2, N);
//	printMatrixVal_com(tmp3, N);


	//g / 2 * [L, rho L_c] * h
	for (i = 0; i < Np2; i++)
	{
		res[i].re += (2.0 * tmp2[i].re - tmp3[i].re) * g * h;
		res[i].im += (2.0 * tmp2[i].im - tmp3[i].im) * g * h;
	}

//	printMatrixVal_com(res, N);

	// L rho
	cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
		N, N, N, &alpha, L, N, x, N, &beta, tmp1, N);
//	//[L rho, L_c]
//	cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasConjTrans,
//		N, N, N, &alpha, tmp1, N, L, N, &beta, tmp2, N);
	cblas_zgemm(CblasRowMajor, CblasConjTrans, CblasNoTrans,
		N, N, N, &alpha, L, N, tmp1, N, &beta, tmp3, N);
	for (i = 0; i < Np2; i++)
	{
		res[i].re += (/*tmp2[i].re*/ - tmp3[i].re) * g * h;
		res[i].im += (/*tmp2[i].im*/ - tmp3[i].im) * g * h;
	}

//	printMatrixVal_com(res, N);
}


