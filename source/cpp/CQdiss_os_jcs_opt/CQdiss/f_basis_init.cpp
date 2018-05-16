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

	time = omp_get_wtime() - init_time;
	cout << "Time of createModel: " << time << endl << endl;

	init_h_0_vector(model, rp, cp, md);
	time = omp_get_wtime() - init_time;
	cout << "Time of init_h_0_vector: " << time << endl << endl;
	if (rp.debug == 1)
	{
		fn = rp.path + "h_0_vec" + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, model->h_0, 16, false);
	}

	init_h_1_vector(model, rp, cp, md);
	time = omp_get_wtime() - init_time;
	cout << "Time of init_h_1_vector: " << time << endl << endl;
	if (rp.debug == 1)
	{
		fn = rp.path + "h_1_vec" + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, model->h_1, 16, false);
	}

	init_H0(model, rp, cp, md);
	init_H1(model, rp, cp, md);

	if (rp.issmtx == 1)
	{
		fn = rp.path + "H0" + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, model->H0, 16, false);

		fn = rp.path + "H1" + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, model->H1, 16, false);
	}

	if (cp.dt == 1)
	{
		init_diss_1(model, rp, cp, md);
	}
	else
	{
		init_diss_1(model, rp, cp, md);
	}

	time = omp_get_wtime() - init_time;
	cout << "time of init_diss_" << to_string(cp.dt) << ": " << time << endl << endl;

	calc_Q_0(model);
	time = omp_get_wtime() - init_time;
	cout << "time of calc_Q_0: " << time << endl << endl;
	if (rp.debug == 1)
	{
		fn = rp.path + "Q_0" + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, model->Q_0, 16, false);
	}

	calc_Q_1(model);
	time = omp_get_wtime() - init_time;
	cout << "time of calc_Q_1: " << time << endl << endl;
	if (rp.debug == 1)
	{
		fn = rp.path + "Q_1" + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, model->Q_1, 16, false);
	}

	calcKs(model);
	time = omp_get_wtime() - init_time;
	cout << "time of calcKs: " << time << endl << endl;

	if (rp.debug == 1)
	{
		fn = rp.path + "Ks" + file_name_suffix(cp, 4);
		save_complex_data(fn, (MKL_Complex16 *)model->Ks, model->N_mat, 15, false);
	}

	calcRs(model);
	time = omp_get_wtime() - init_time;
	cout << "time of calcRs: " << time << endl << endl;
	if (rp.debug == 1)
	{
		fn = rp.path + "Rs" + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, model->Rs, 16, false);
	}

	calc_G_0_s(model, cp);
	time = omp_get_wtime() - init_time;
	cout << "time of calc_G_0_s: " << time << endl << endl;

	if (rp.debug == 1)
	{
		fn = rp.path + "G_0_s" + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, model->G_0_s, 16, false);
	}

	calc_G_1_s(model, cp);
	time = omp_get_wtime() - init_time;
	cout << "time of calc_G_1_s: " << time << endl << endl;

	if (rp.debug == 1)
	{
		fn = rp.path + "G_1_s" + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, model->G_1_s, 16, false);
	}

	initRhoODE(model, rp, cp, md);
	time = omp_get_wtime() - init_time;
	cout << "time of initRhoODE: " << time << endl << endl;

	if (rp.debug == 1)
	{
		fn = rp.path + "init_rho_f" + file_name_suffix(cp, 4);
		save_complex_data(fn, (MKL_Complex16 *)model->RhoF, model->N_mat, 16, false);
	}
}
