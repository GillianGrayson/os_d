#include <omp.h>

#include <stdio.h>
#include <stdlib.h>
#include "read_config.h"
#include "Model.h"
#include "InitFs.h"
#include "initH.h"
#include "Init_a1_a2.h"
#include "CalcODE_rk.h"


#include "InitFs.h"
#include "Init_f_d.h"
#include "TranspFs.h"

#include "CalcQs.h"
#include "CalcKs.h"
#include "CalcRs.h"
#include "CalcGs.h"

#include "linSolv.h"
#include "calcODE.h"
#include "calcODE_real.h"

#include "CalcRho.h"

bool ismap = false;
bool ismap_print = false;
double map_size = 0.0;

double trace(dcomplex * mat, int N)
{
	double sum;
	sum = 0.0;
	for (int i = 0; i < N; i++)
	{
		sum += mat[i * N + i].re; //* mat[i * N + i].re + mat[i * N + i].im * mat[i * N + i].im;
	}
	return sum;
}

int main(int argc, char ** argv)
{
//################################# init RK #####################################
	printf("current path (%s)\n\n", argv[0]);
	if (argc > 1)
		omp_set_num_threads(atoi(argv[1]));

	double time, all_time;
	ConfigParam param;
	read_config(param, "config.txt");

	Model * model;
	model = createModel(param.N, param);

	crsMatrix * H = createHmatrix(model, false);
	crsMatrix * He = createHe_matrix(model);
//	printMatrixVal(H);
//	printMatrixVal(He);

	crsMatrix * a1_mat = createA1mat(model->N);
	crsMatrix * a2_mat = createA2mat(model->N);
	crsMatrix * L = new crsMatrix(a1_mat->N, a1_mat->NZ + a2_mat->NZ);
	dcomplex alfa;
	alfa.re = 1.0;
	alfa.im = 0.0;
	dcomplex betta;
	betta.re = 0.0;
	betta.im = -1.0;
	SparseMKLAdd(*a1_mat, betta, *a2_mat, *L);
	
//	crsMatrix * Lt = new crsMatrix(*L);
//	Transpose(*L, *Lt);
//	printMatrixVal(a1_mat);
//	printMatrixVal(a2_mat);
//	printMatrixVal(L);
//	printMatrixVal(Lt);
//	crsMatrix *A = new crsMatrix;
//	SparseMKLMult(*L, *Lt, *A, true);
//	printMatrixVal(A);


	dcomplex * rho = initRhoODE_rk(model);
	dcomplex * rho_rk_p = new dcomplex[(model->N + 1)*(model->N + 1)];
	dcomplex * rho_f_p = new dcomplex[(model->N + 1)*(model->N + 1)];
	dcomplex * rho_f_c = new dcomplex[(model->N + 1)*(model->N + 1)];
	double v = trace(rho, model->N + 1);
	printf("(%1.2lf)", v);

// ############################################ init F basis ######################################

	time = omp_get_wtime();
	// init_h_vector(model);
	init_h_vector_opt(model);
	time = omp_get_wtime() - time;
	printf("init_h_vector: %2.4lf\n", time);
	fflush(stdout);

	time = omp_get_wtime();
	init_he_vector_opt(model);
	time = omp_get_wtime() - time;
	printf("init_he_vector: %2.4lf\n", time);
	fflush(stdout);

	time = omp_get_wtime();
	//  init_a1_a2(model);
	init_a1_a2_opt(model);
	time = omp_get_wtime() - time;
	printf("init_a1_a2: %2.4lf\n", time);
	fflush(stdout);

	time = omp_get_wtime();
	calcQs(model);
	time = omp_get_wtime() - time;
	printf("calcQs: %2.4lf\n", time);
	fflush(stdout);

	time = omp_get_wtime();
	calcQEs(model);
	time = omp_get_wtime() - time;
	printf("calcQEs: %2.4lf\n", time);
	fflush(stdout);

	time = omp_get_wtime();
	calcKs(model);
	time = omp_get_wtime() - time;
	printf("calcKs: %2.4lf\n", time);
	fflush(stdout);

	time = omp_get_wtime();
	calcRs(model);
	time = omp_get_wtime() - time;
	printf("calcRs: %2.4lf\n", time);
	fflush(stdout);

	time = omp_get_wtime();
	calcGs(model);
	time = omp_get_wtime() - time;
	printf("calcGs: %2.4lf\n", time);
	fflush(stdout);

	time = omp_get_wtime();
	initRhoODE(model);
	//    printVectorVal(model->RhoF, model->N_mat);
	time = omp_get_wtime() - time;
	printf("initRhoODE: %2.4lf\n", time);

	complex_to_real(model->Gs->Value, model->Gs->NZ);
	complex_to_real(model->QEs->Value, model->QEs->NZ);
	complex_to_real(model->Ks, model->N_mat);
	complex_to_real(model->RhoF, model->N_mat);
	toOneBase(*(model->Gs));
	toOneBase(*(model->QEs));


	for (int i = 0; i < (model->N + 1)*(model->N + 1); i++)
	{
		rho_rk_p[i] = rho[i];
		rho_f_p[i].re = rho_f_p[i].im = 0;
	}
	calcRho_fill(model);
	for (int i = 0; i < model->Rho->N; i++)
	{
		for (int j = model->Rho->RowIndex[i]; j < model->Rho->RowIndex[i + 1]; j++)
		{
			int c = model->Rho->Col[j];
			rho_f_p[i * (model->N + 1) + c].re = model->Rho->Value[j].re;
			rho_f_p[i * (model->N + 1) + c].im = model->Rho->Value[j].im;
		}
	}
	FILE  *f_rk, *f_ff, *f_diff;
	
	f_rk = fopen("rk_diff_o.txt", "w");
	f_ff = fopen("ff_diff_o.txt", "w");
	f_diff = fopen("diff_o.txt", "w");

	for (int itr = 0; itr < model->conf.N_T; itr++)
	{
		calcODE_rk(model, H, He, L, rho, model->conf.h,
			model->conf.NSTEP, itr * model->conf.T);
		calcODE_real(model, model->conf.h,
			model->conf.NSTEP, itr * model->conf.T);

		real_to_complex(model->RhoF, model->N_mat);
//		calcRho(model);
		calcRho_fill(model);
		complex_to_real(model->RhoF, model->N_mat);

		v = trace(rho, model->N + 1);
		//for (int i = 0; i < (model->N + 1) * (model->N + 1); i++)
		//{
		//	rho[i].re /= sqrt(v);
		//	rho[i].im /= sqrt(v);
		//}
		//printMatrixVal_com(rho, model->N + 1);
		//printMatrixVal(model->Rho);
		
		for (int i = 0; i < (model->N + 1)*(model->N + 1); i++)
		{
			rho_f_c[i].re = rho_f_c[i].im = 0;
		}
		for (int i = 0; i < model->Rho->N; i++)
		{
			for (int j = model->Rho->RowIndex[i]; j < model->Rho->RowIndex[i + 1]; j++)
			{
				int c = model->Rho->Col[j];
				rho_f_c[i * (model->N + 1) + c].re = model->Rho->Value[j].re;
				rho_f_c[i * (model->N + 1) + c].im = model->Rho->Value[j].im;
			}
		}

//		printMatrixVal_com(rho_f_c, model->N + 1);
//		printMatrixVal_com(rho, model->N + 1);
		
		dcomplex max_val;
		max_val.re = 0.0;
		max_val.im = 0.0;
		double mod = max_val.re * max_val.re + max_val.im * max_val.im;

		dcomplex max_diff;
		max_diff.re = max_diff.im = 0.0;
		dcomplex max_diff_rk;
		max_diff_rk.re = max_diff_rk.im = 0.0;
		dcomplex max_diff_f;
		max_diff_f.re = max_diff_f.im = 0.0;
		for (int i = 0; i < (model->N + 1)*(model->N + 1); i++)
		{
			if (max_diff.re < fabs(rho_f_c[i].re - rho[i].re)) max_diff.re = fabs(rho_f_c[i].re - rho[i].re);
			if (max_diff.im < fabs(rho_f_c[i].im - rho[i].im)) max_diff.im = fabs(rho_f_c[i].im - rho[i].im);
			if (max_diff_rk.re < fabs(rho_rk_p[i].re - rho[i].re)) max_diff_rk.re = fabs(rho_rk_p[i].re - rho[i].re);
			if (max_diff_rk.im < fabs(rho_rk_p[i].im - rho[i].im)) max_diff_rk.im = fabs(rho_rk_p[i].im - rho[i].im);
			if (max_diff_f.re < fabs(rho_f_p[i].re - rho_f_c[i].re)) max_diff_f.re = fabs(rho_f_p[i].re - rho_f_c[i].re);
			if (max_diff_f.im < fabs(rho_f_p[i].im - rho_f_c[i].im)) max_diff_f.im = fabs(rho_f_p[i].im - rho_f_c[i].im);
			
			if (mod < (rho[i].re * rho[i].re + rho[i].im * rho[i].im))
			{
				max_val.re = rho[i].re;
				max_val.im = rho[i].im;
				mod = max_val.re * max_val.re + max_val.im * max_val.im;
			}
		}

		mod = sqrt(mod);
		printf("%1.8lf \n", mod);

		dcomplex tmp;
		for (int i = 0; i < (model->N + 1); i++)
		{
			int ind = (model->N + 1) * i + i;
			tmp.re = rho_rk_p[ind].re - rho[ind].re;
			tmp.im = rho_rk_p[ind].im - rho[ind].im;
			fprintf(f_rk, "%1.8lf ", sqrt(tmp.re * tmp.re + tmp.im * tmp.im) / mod);

			tmp.re = rho_f_p[ind].re - rho_f_c[ind].re;
			tmp.im = rho_f_p[ind].im - rho_f_c[ind].im;
			fprintf(f_ff, "%1.8lf ", sqrt(tmp.re * tmp.re + tmp.im * tmp.im) / mod);

			tmp.re = fabs(rho_f_c[ind].re - rho[ind].re);
			tmp.im = fabs(rho_f_c[ind].im - rho[ind].im);
			fprintf(f_diff, "%1.16lf ", fmax(tmp.re, tmp.im)/ mod);

		}
		fprintf(f_rk, "\n"); fflush(f_rk);
		fprintf(f_ff, "\n"); fflush(f_ff);
		fprintf(f_diff, "\n"); fflush(f_diff);

		for (int i = 0; i < (model->N + 1)*(model->N + 1); i++)
		{
			rho_rk_p[i] = rho[i];
			rho_f_p[i] = rho_f_c[i];
		}

		printf("tr %1.6lf diff:rk(%1.8lf, %1.8lf)f(%1.8lf, %1.8lf)f_vs_rk(%1.16lf, %1.16lf)\n", v, 
			max_diff_rk.re, max_diff_rk.im, max_diff_f.re, max_diff_f.im, max_diff.re, max_diff.im);
	}

	fclose(f_rk);
	fclose(f_ff);
	fclose(f_diff);


	toZeroBase(*(model->Gs));
	toZeroBase(*(model->QEs));
	real_to_complex(model->Gs->Value, model->Gs->NZ);
	real_to_complex(model->QEs->Value, model->QEs->NZ);
	real_to_complex(model->Ks, model->N_mat);
	real_to_complex(model->RhoF, model->N_mat);


	//for (int i = 0; i < (model->N + 1) * (model->N + 1); i++)
	//{
	//	rho[i].re /= sqrt(v);
	//	rho[i].im /= sqrt(v);
	//}
	//v = trace(rho, model->N + 1);
	//printf("(%1.2lf)", v);
	//printMatrixVal_com(rho, model->N + 1);

	saveMatrixVal_com("rho.txt", rho, model->N + 1);
	saveAbsMatrixVal_com("absRho.txt", rho, model->N + 1);

	delete a1_mat;
	delete a2_mat;
	delete L;
	freeModel(model);
	return 0;
}