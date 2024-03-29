#include <omp.h>
#include <mkl.h>
#include <math.h>

#include "read_config.h"
#include "Model.h"
#include "InitFs.h"
#include "initH.h"
#include "Init_a1_a2.h"
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

#include "CalcEigs.h"
#include "CalcTraseRO2.h"

#include "testODE.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>

using namespace std;

double number_of_allocs = 0.0;

std::map<void *, double> mem_map;

FILE * mem_time;

bool ismap = false;
bool ismap_print = false;
double map_size = 0.0;

bool program = true;
bool isEnd = false;
void* operator new(std::size_t size) throw(std::bad_alloc){
	number_of_allocs += size / 1024.0 / 1024.0;
	void *p = malloc(size);

	if (mem_time != NULL)
	{
		if (ismap)
		{
			map_size += size / 1024.0 / 1024.0;
		}
		if (ismap_print)
		{
			//fprintf(mem_time, "%lf \n", map_size);
			//fprintf(mem_time, "%lf \n", -map_size);
			fflush(mem_time);
			ismap = false;
			ismap_print = false;
		}
		if (program && !ismap)
		{
			program = false;
			std::pair<void *, double> pr(p, size / 1024.0 / 1024.0);
			program = false;
			mem_map.insert(pr);

			//fprintf(mem_time, "%lf \n", mem_map[p]);
			fflush(mem_time);
		}
	}
	program = true;
	if (!p) throw std::bad_alloc();
	return p;
}
void* operator new   (std::size_t size, const std::nothrow_t&) throw(){
	number_of_allocs += size / 1024.0 / 1024.0;
	void *p = malloc(size);

	if (mem_time != NULL)
	{
		if (ismap)
		{
			map_size += size / 1024.0 / 1024.0;
		}
		if (ismap_print)
		{
			//fprintf(mem_time, "%lf \n", map_size);
			//fprintf(mem_time, "%lf \n", -map_size);
			fflush(mem_time);
			ismap = false;
			ismap_print = false;
		}
		if (program && !ismap)
		{
			program = false;
			std::pair<void *, double> pr(p, size / 1024.0 / 1024.0);
			program = false;
			mem_map.insert(pr);

			//fprintf(mem_time, "%lf \n", mem_map[p]);
			fflush(mem_time);
		}
	}
	program = true;

	return p;
}

void* operator new  [](std::size_t size) throw(std::bad_alloc) {
	number_of_allocs += size / 1024.0 / 1024.0;
	void *p = malloc(size);
	
	if (mem_time != NULL)
	{
		if (ismap)
		{
			map_size += size / 1024.0 / 1024.0;
		}
		if (ismap_print)
		{
			//fprintf(mem_time, "%lf \n", map_size);
			//fprintf(mem_time, "%lf \n", -map_size);
			fflush(mem_time);
			ismap = false;
			ismap_print = false;
		}
		if (program && !ismap)
		{
			program = false;
			std::pair<void *, double> pr(p, size / 1024.0 / 1024.0);
			program = false;
			mem_map.insert(pr);

			//fprintf(mem_time, "%lf \n", mem_map[p]);
			fflush(mem_time);
		}
	}
	program = true;

	if (!p) throw std::bad_alloc();
	return p;
}

void* operator new  [](std::size_t size, const std::nothrow_t&) throw() {
	number_of_allocs += size / 1024.0 / 1024.0;
	void *p = malloc(size);
	
	if (mem_time != NULL)
	{
		if (ismap)
		{
			map_size += size / 1024.0 / 1024.0;
		}
		if (ismap_print)
		{
			//fprintf(mem_time, "%lf \n", map_size);
			//fprintf(mem_time, "%lf \n", -map_size);
			fflush(mem_time);
			ismap = false;
		}
		if (program && !ismap)
		{
			program = false;
			std::pair<void *, double> pr(p, size / 1024.0 / 1024.0);
			program = false;
			mem_map.insert(pr);

			//fprintf(mem_time, "%lf \n", mem_map[p]);
			fflush(mem_time);
		}
	}
	program = true;

	return p;
}


void operator delete(void* ptr) throw() {  
	if (!isEnd) if (mem_map.count(ptr) > 0){
		//fprintf(mem_time, "%lf \n", -mem_map[ptr]);
		fflush(mem_time);
	}
free(ptr);
}
void operator delete (void* ptr, const std::nothrow_t&) throw() {  
	if (!isEnd)if (mem_map.count(ptr) > 0) {
		//fprintf(mem_time, "%lf \n", -mem_map[ptr]);
		fflush(mem_time);
	}
free(ptr);
}
void operator delete[](void* ptr) throw() { 
	if (!isEnd)if (mem_map.count(ptr) > 0) {
		//fprintf(mem_time, "%lf \n", -mem_map[ptr]);
		fflush(mem_time);
	}
free(ptr);
}
void operator delete[](void* ptr, const std::nothrow_t&) throw() { 
	if (!isEnd)if (mem_map.count(ptr) > 0) {
		//fprintf(mem_time, "%lf \n", -mem_map[ptr]);
		fflush(mem_time);
	}
free(ptr);
}
//*/

int main(int argc, char ** argv)
{
	mem_time = fopen("mem_time.txt", "w");
	//printf("current path (%s)\n\n", argv[0]);
	if (argc >1)
		omp_set_num_threads(atoi(argv[1]));

	double time, all_time;
	ConfigParam param;
	read_config(param, "config.txt");

	Model * model;
	//fprintf(mem_time, "model = createModel(param.N, param); \n");

	model = createModel(param.N, param);
	FILE * memlog = model->memlog;

	fprintf(memlog, "model = createModel(param.N, param); %lf\n", number_of_allocs);
	fflush(memlog);
	
	all_time = omp_get_wtime();

	number_of_allocs = 0;

	//fprintf(mem_time, "init_h_vector(model); \n");
	time = omp_get_wtime();
	init_h_vector_opt(model);
	time = omp_get_wtime() - time;
	printf("init_h_vector: %2.4lf\n", time);
	fflush(stdout);

	fprintf(memlog, "init_h_vector(model); %lf\n", number_of_allocs);
	fflush(memlog);
	number_of_allocs = 0;
	
	//fprintf(mem_time, "init_he_vector(model); \n");
	time = omp_get_wtime();
	init_he_vector_opt(model);
	time = omp_get_wtime() - time;
	printf("init_he_vector: %2.4lf\n", time);
	fflush(stdout);

	fprintf(memlog, "init_he_vector(model); %lf\n", number_of_allocs);
	fflush(memlog);
	number_of_allocs = 0;
	
	//fprintf(mem_time, "init_a1_a2(model); \n");
	time = omp_get_wtime();
	init_a1_a2_opt(model);
	time = omp_get_wtime() - time;
	printf("init_a1_a2: %2.4lf\n", time);
	fflush(stdout);

	fprintf(memlog, "init_a1_a2(model); %lf\n", number_of_allocs);
	fflush(memlog);
	number_of_allocs = 0;

	//fprintf(mem_time, "init_f_d_valentin(model); \n");
	time = omp_get_wtime();
	init_f_d_valentin(model);
	time = omp_get_wtime() - time;
	printf("init_f_d (valentin): %2.4lf\n", time);
	fflush(stdout);
	fprintf(memlog, "init_f_d_valentin(model); %lf\n", number_of_allocs);
	fflush(memlog);
	number_of_allocs = 0;

	//fprintf(mem_time, "calcQs(model); \n");
	time = omp_get_wtime();
	calcQs(model);
	time = omp_get_wtime() - time;
	printf("calcQs: %2.4lf\n", time);
	fflush(stdout);
	fprintf(memlog, "calcQs(model); %lf\n", number_of_allocs);
	fflush(memlog);
	number_of_allocs = 0;
	
	//fprintf(mem_time, "calcQEs(model); \n");
	time = omp_get_wtime();
	calcQEs(model);
	time = omp_get_wtime() - time;
	printf("calcQEs: %2.4lf\n", time);
	fflush(stdout);
	fprintf(memlog, "calcQEs(model); %lf\n", number_of_allocs);
	fflush(memlog);
	number_of_allocs = 0;
	
	
	//fprintf(mem_time, "calcKs(model); \n");
	time = omp_get_wtime();
	calcKs(model);
	time = omp_get_wtime() - time;
	printf("calcKs: %2.4lf\n", time);
	fflush(stdout);
	fprintf(memlog, "calcKs(model); %lf\n", number_of_allocs);
	fflush(memlog);
	number_of_allocs = 0;

	//fprintf(mem_time, "calcRs(model); \n");
	time = omp_get_wtime();
	calcRs(model);
	time = omp_get_wtime() - time;
	printf("calcRs: %2.4lf\n", time);
	fflush(stdout);
	fprintf(memlog, "calcRs(model); %lf\n", number_of_allocs);
	fflush(memlog);
	number_of_allocs = 0;

	//fprintf(mem_time, "calcGs(model); \n");
	time = omp_get_wtime();
	calcGs(model);
	time = omp_get_wtime() - time;
	printf("calcGs: %2.4lf\n", time);
	fflush(stdout);
	fprintf(memlog, "calcGs(model); %lf\n", number_of_allocs);
	fflush(memlog);
	number_of_allocs = 0;

	char file_name[512];

	//saveMatrix("Gs.txt", model->Gs);
	//save_complex_vector("Ks.txt", model->Ks, model->N_mat);

	//  saveMatrix_coor("Gs_p.txt", model->Gs);
	//  printMatrixVal(model->Gs);
	//  printVectorVal(model->Ks, model->N_mat);
	//  printMatrixVal(model->Rs);
	//  printMatrix(model->Rs);
	//  printMatrixVal(model->Gs);

	//fprintf(mem_time, "calcODE(model); \n");

	if (model->conf.hasDriving == 0)
	{
		time = omp_get_wtime();
		linSolv(model);
		time = omp_get_wtime() - time;
		printf("linSolv: %2.4lf\n", time);
		//saveMatrix("Gs.txt", model->Gs);
	}
	else
	{
		time = omp_get_wtime();
		initRhoODE(model);

		//save_complex_vector("rho_f_init.txt", model->RhoF, model->N_mat);

		time = omp_get_wtime() - time;
		printf("initRhoODE: %2.4lf\n", time);

		time = omp_get_wtime();

		complex_to_real(model->Gs->Value, model->Gs->NZ);
		complex_to_real(model->QEs->Value, model->QEs->NZ);
		complex_to_real(model->Ks, model->N_mat);
		complex_to_real(model->RhoF, model->N_mat);
		toOneBase(*(model->Gs));
		toOneBase(*(model->QEs));

		dcomplex diff_it;

		if (model->conf.deep_dump == 2)
		{
			string fn = "diag_rho_deep_evo.txt";
			ofstream ofs = ofstream(fn);

			if (ofs.is_open())
			{
				ofs.close();
			}

			model->conf.deep_dump = 20;
		}

		for (int itr = 0; itr < model->conf.N_T; itr++)
		{
			calcODE_real(model, model->conf.h,
				model->conf.NSTEP, itr * model->conf.T);

			diff_it = calcDiffIter(model);
			printf("diff on %d is %0.16le %0.16le \n", itr, diff_it.re, diff_it.im);

			if (model->conf.deep_dump == 1)
			{
				real_to_complex(model->RhoF, model->N_mat);

				sprintf(file_name, "rho_p%d.txt", itr);
				calcRho_fill(model);
				saveMatrix(file_name, model->Rho);

				complex_to_real(model->RhoF, model->N_mat);
			}
		}

		if (model->conf.deep_dump == 20)
		{
			model->conf.deep_dump = 2;

			for (int itr = model->conf.N_T; itr < model->conf.N_T + 10; itr++)
			{
				calcODE_real(model, model->conf.h,
					model->conf.NSTEP, itr * model->conf.T);
			}
		}

		toZeroBase(*(model->Gs));
		toZeroBase(*(model->QEs));
		real_to_complex(model->Gs->Value, model->Gs->NZ);
		real_to_complex(model->QEs->Value, model->QEs->NZ);
		real_to_complex(model->Ks, model->N_mat);
		real_to_complex(model->RhoF, model->N_mat);

		//save_complex_vector("rho_f_fin.txt", model->RhoF, model->N_mat);

		time = omp_get_wtime() - time;

		printf("calcODE: %2.4lf\n", time);
		fflush(stdout);	
	}

	fprintf(memlog, "calcODE_real; %lf\n", number_of_allocs);
	number_of_allocs = 0;

	linSolvCheck(model);

	time = omp_get_wtime();
	calcRho_fill(model);
	//check_rho_evals(model);
	time = omp_get_wtime() - time;
	printf("calcRho: %2.4lf\n", time);
	fprintf(memlog, "calcRho(model); %lf\n", number_of_allocs);
	fflush(memlog);
	number_of_allocs = 0;

	saveMatrix("rho.txt", model->Rho);
	save_mtx_diag("rho_diag.txt", model->Rho);

	if (model->conf.deep_dump == 1)
	{
		sprintf(file_name, "rho_p%d.txt", model->conf.N_T);
		calcRho_fill(model);
		saveMatrix(file_name, model->Rho);
	}

	all_time = omp_get_wtime() - all_time;

	printf("all time: %2.4lf\n", all_time);
	dcomplex tr = trace(*(model->Rho));
	printf("Trace: %lf\n", tr.re);

	//fprintf(mem_time, "calcTraseRO2(model); \n");
	time = omp_get_wtime();
	calcTraseRO2(model);
	time = omp_get_wtime() - time;
	printf("calcTraseRO2: %2.4lf\n", time);
	fprintf(memlog, "calcTraseRO2(model); %lf\n", number_of_allocs);
	fflush(memlog);
	number_of_allocs = 0;
	calc_negativity_final(model);

	if (model->conf.CalcEig == 1)
	{
		time = omp_get_wtime();
		calcEig(model);
		time = omp_get_wtime() - time;
		printf("calcEig: %2.4lf\n", time);
	}

	if (model->conf.hasDriving == 1)
	{
		if (model->conf.avg_purity == 1)
		{
			//  ================== Last period ===================== 

			int num_dumps = 100;
			int dump_iter_step = model->conf.NSTEP / num_dumps;

			double * mult_purity = new double[num_dumps];
			double * mult_negativity = new double[num_dumps];

			// regular dumps
			for (int dump_id = 0; dump_id < num_dumps; dump_id++)
			{
				complex_to_real(model->Gs->Value, model->Gs->NZ);
				complex_to_real(model->QEs->Value, model->QEs->NZ);
				complex_to_real(model->Ks, model->N_mat);
				complex_to_real(model->RhoF, model->N_mat);
				toOneBase(*(model->Gs));
				toOneBase(*(model->QEs));

				double int_time = model->conf.T * double(model->conf.N_T) + double(dump_id * dump_iter_step) * model->conf.h;
				
				//printf("int_time: %0.16le\n", int_time);

				calcODE_real(model, model->conf.h, dump_iter_step, int_time);

				toZeroBase(*(model->Gs));
				toZeroBase(*(model->QEs));
				real_to_complex(model->Gs->Value, model->Gs->NZ);
				real_to_complex(model->QEs->Value, model->QEs->NZ);
				real_to_complex(model->Ks, model->N_mat);
				real_to_complex(model->RhoF, model->N_mat);

				calcRho_fill(model);
				//check_rho_evals(model);

				double curr_trace = calc_purity(model);
				double curr_neg = calc_negativity(model);

				mult_purity[dump_id] = curr_trace;
				mult_negativity[dump_id] = curr_neg;
			}

			double avg_purity = 0.0;
			double avg_negativity = 0.0;
			for (int dump_id = 0; dump_id < num_dumps; dump_id++)
			{
				avg_purity += mult_purity[dump_id];
				avg_negativity += mult_negativity[dump_id];
			}
			avg_purity = avg_purity / double(num_dumps);
			avg_negativity = avg_negativity / double(num_dumps);

			FILE * f = fopen("purity_avg.txt", "w");
			if (f != NULL)
			{
				fprintf(f, "%0.16le\n", avg_purity);
			}
			fclose(f);

			f = fopen("negativity_avg.txt", "w");
			if (f != NULL)
			{
				fprintf(f, "%0.16le\n", avg_negativity);
			}
			fclose(f);

			delete[] mult_purity;
			delete[] mult_negativity;
		}
	}
	
	//fprintf(mem_time, "other \n");

	//saveVectorVal("RhoF.txt", model->RhoF, 1, model->N_mat);
	//saveMatrixVal("rho_f.txt", model->Rho);
	//saveAbsMatrixVal("absRho.txt", model->Rho);
	//saveAngleMatrixVal("angleRho.txt", model->Rho);

	printf("\n");

	if (model->conf.hasDriving == 1)
	{
		if (model->conf.multiplicators == 1)
		{
			int space_size = model->N_mat;
			int num_mults = model->N_mat;
			int full_size = model->N_mat + 1;

			dcomplex * monodromy_mtx = new dcomplex[full_size * full_size];
			for (int st_id_1 = 0; st_id_1 < full_size; st_id_1++)
			{
				for (int st_id_2 = 0; st_id_2 < full_size; st_id_2++)
				{
					monodromy_mtx[st_id_1 * full_size + st_id_2].re = 0.0;
					monodromy_mtx[st_id_1 * full_size + st_id_2].im = 0.0;
				}
			}
			monodromy_mtx[0].re = 1.0;

			dcomplex * mults = new dcomplex[full_size];
			for (int mult_id = 0; mult_id < full_size; mult_id++)
			{
				mults[mult_id].re = 0.0;
				mults[mult_id].im = 0.0;
			}

			char file_name[256];

			// ======================
			// Init multiplicator
			// ======================

			printf("first_multiplicator\n");

			init_first_multiplicator(model);
			//calcRho_fill(model);
			//check_rho_evals(model);

			//sprintf(file_name, "vec_first_init.txt");
			//save_complex_vector(file_name, model->RhoF, model->N_mat);
			//sprintf(file_name, "rho_first_init.txt");
			//saveMatrix(file_name, model->Rho);

			complex_to_real(model->Gs->Value, model->Gs->NZ);
			complex_to_real(model->QEs->Value, model->QEs->NZ);
			complex_to_real(model->Ks, model->N_mat);
			complex_to_real(model->RhoF, model->N_mat);
			toOneBase(*(model->Gs));
			toOneBase(*(model->QEs));

			calcODE_real(model, model->conf.h, model->conf.NSTEP, 0.0);

			toZeroBase(*(model->Gs));
			toZeroBase(*(model->QEs));
			real_to_complex(model->Gs->Value, model->Gs->NZ);
			real_to_complex(model->QEs->Value, model->QEs->NZ);
			real_to_complex(model->Ks, model->N_mat);
			real_to_complex(model->RhoF, model->N_mat);

			//sprintf(file_name, "vec_%d_fin.txt", mult_id);
			//save_complex_vector(file_name, model->RhoF, model->N_mat);
			//calcRho_fill(model);
			//check_rho_evals(model);
			//sprintf(file_name, "rho_first_fin.txt");
			//saveMatrix(file_name, model->Rho);

			set_monodromy_state(model, 0, monodromy_mtx);

			for (int mult_id = 0; mult_id < num_mults; mult_id++)
			{
				printf("mult_id: %d\n", mult_id);

				init_multiplicator(model, mult_id);

				//calcRho_fill(model);
				//check_rho_evals(model);
				//sprintf(file_name, "vec_%d_init.txt", mult_id);
				//save_complex_vector(file_name, model->RhoF, model->N_mat);
				//if (mult_id == -1)
				//{
				//	sprintf(file_name, "rho_%d_init.txt", mult_id);
				//	saveMatrix(file_name, model->Rho);
				//}

				complex_to_real(model->Gs->Value, model->Gs->NZ);
				complex_to_real(model->QEs->Value, model->QEs->NZ);
				complex_to_real(model->Ks, model->N_mat);
				complex_to_real(model->RhoF, model->N_mat);
				toOneBase(*(model->Gs));
				toOneBase(*(model->QEs));

				calcODE_real(model, model->conf.h, model->conf.NSTEP, 0.0);

				toZeroBase(*(model->Gs));
				toZeroBase(*(model->QEs));
				real_to_complex(model->Gs->Value, model->Gs->NZ);
				real_to_complex(model->QEs->Value, model->QEs->NZ);
				real_to_complex(model->Ks, model->N_mat);
				real_to_complex(model->RhoF, model->N_mat);

				set_monodromy_state(model, mult_id + 1, monodromy_mtx);

				//sprintf(file_name, "vec_%d_fin.txt", mult_id);
				//save_complex_vector(file_name, model->RhoF, model->N_mat);
				//calcRho_fill(model);
				//check_rho_evals(model);
				//if (mult_id == -1)
				//{
				//	sprintf(file_name, "rho_%d_fin.txt", mult_id);
				//	saveMatrix(file_name, model->Rho);
				//}
			}

			int info;
			info = LAPACKE_zgeev(LAPACK_ROW_MAJOR, 'N', 'N', full_size, (MKL_Complex16 *)monodromy_mtx, full_size,
				(MKL_Complex16 *)mults, NULL, full_size, NULL, full_size);
			/* Check for convergence */
			if (info > 0) {
				printf("The algorithm failed to compute eigenvalues.\n");
				exit(1);
			}

			FILE * file = fopen("multiplicators.txt", "w");

			for (int mult_id = 0; mult_id < full_size; mult_id++)
			{
				fprintf(file, "%1.16le %1.16le\n", mults[mult_id].re, mults[mult_id].im);
			}

			fclose(file);

			delete[] mults;
			delete[] monodromy_mtx;
		}
	}

	freeModel(model);
	isEnd = true;
	fclose(mem_time);
	return 0;
}
