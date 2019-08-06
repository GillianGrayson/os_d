#include "f_basis_init.h"
#include "initH.h"
#include "Init_a1_a2.h"
#include "CalcQs.h"
#include "CalcKs.h"
#include "CalcRs.h"
#include "CalcGs.h"
#include "calcODE.h"
#include "CalcGs.h"

void f_basis_init(Model* model, RunParam &rp, ConfigParam &cp, MainData &md)
{
	string fn;

	double time = omp_get_wtime();
	double init_time = time;

	init_h_base_vector(model, rp, cp, md);
	time = omp_get_wtime() - init_time;
	cout << "Time of init_h_base_vector: " << time << endl << endl;
	if (rp.debug == 1)
	{
		fn = "h_base_vec" + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, model->h_base, 16, false);
	}

	init_h_drv_vector(model, rp, cp, md);
	time = omp_get_wtime() - init_time;
	cout << "Time of init_h_drv_vector: " << time << endl << endl;
	if (rp.debug == 1)
	{
		fn = "h_drv_vector" + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, model->h_drv, 16, false);
	}

	if (cp.diss_type == 0)
	{
		init_diss(model, rp, cp, md);
	}
	else
	{
		init_diss(model, rp, cp, md);
	}

	time = omp_get_wtime() - init_time;
	cout << "time of init_diss: " << time << endl << endl;

	calc_Qs_base(model);
	time = omp_get_wtime() - init_time;
	cout << "time of calc_Qs_base: " << time << endl << endl;
	if (rp.debug == 1)
	{
		fn = "Qs_base" + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, model->Qs_base, 16, false);
	}

	calc_Qs_drv(model);
	time = omp_get_wtime() - init_time;
	cout << "time of calc_Qs_drv: " << time << endl << endl;
	if (rp.debug == 1)
	{
		fn = "Qs_drv" + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, model->Qs_drv, 16, false);
	}

	calcKs_dimer(model);
	time = omp_get_wtime() - init_time;
	cout << "time of calcKs: " << time << endl << endl;

	if (rp.debug == 1)
	{
		fn = "Ks" + file_name_suffix(cp, 4);
		save_complex_data(fn, (MKL_Complex16 *)model->Ks, model->N_mat, 15, false);
	}

	calcRs_dimer(model);
	time = omp_get_wtime() - init_time;
	cout << "time of calcRs: " << time << endl << endl;
	if (rp.debug == 1)
	{
		fn = "Rs" + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, model->Rs, 16, false);
	}

	calc_Gs(model, cp);
	time = omp_get_wtime() - init_time;
	cout << "time of calc_Gs: " << time << endl << endl;

	if (rp.debug == 1)
	{
		fn = "Gs" + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, model->Gs, 16, false);
	}

	initRhoODE(model, rp, cp, md);
	time = omp_get_wtime() - init_time;
	cout << "time of initRhoODE: " << time << endl << endl;

	if (rp.debug == 1)
	{
		fn = "init_rho_f" + file_name_suffix(cp, 4);
		save_complex_data(fn, (MKL_Complex16 *)model->RhoF, model->N_mat, 16, false);
	}
}
