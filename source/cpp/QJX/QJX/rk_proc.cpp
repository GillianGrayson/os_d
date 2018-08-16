#include "rk_proc.h"
#include "experiment.h"

void init_rk_data(AllData * ad)
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int num_threads = rp->num_threads;
	int num_trajectories = cp->num_trajectories;
	int sys_size = md->sys_size;

	ed->k1 = new MKL_Complex16*[num_threads];
	ed->k2 = new MKL_Complex16*[num_threads];
	ed->k3 = new MKL_Complex16*[num_threads];
	ed->k4 = new MKL_Complex16*[num_threads];
	ed->args = new MKL_Complex16*[num_threads];
	ed->non_drv_tmp = new MKL_Complex16*[num_threads];
	ed->drv_tmp = new MKL_Complex16*[num_threads];
	for (int th_id = 0; th_id < num_threads; th_id++)
	{
		ed->k1[th_id] = new MKL_Complex16[sys_size];
		ed->k2[th_id] = new MKL_Complex16[sys_size];
		ed->k3[th_id] = new MKL_Complex16[sys_size];
		ed->k4[th_id] = new MKL_Complex16[sys_size];
		ed->args[th_id] = new MKL_Complex16[sys_size];
		ed->non_drv_tmp[th_id] = new MKL_Complex16[sys_size];
		ed->drv_tmp[th_id] = new MKL_Complex16[sys_size];
		for (int st_id = 0; st_id < sys_size; st_id++)
		{
			ed->k1[th_id][st_id].real = 0.0;
			ed->k1[th_id][st_id].imag = 0.0;
			ed->k2[th_id][st_id].real = 0.0;
			ed->k2[th_id][st_id].imag = 0.0;
			ed->k3[th_id][st_id].real = 0.0;
			ed->k3[th_id][st_id].imag = 0.0;
			ed->k4[th_id][st_id].real = 0.0;
			ed->k4[th_id][st_id].imag = 0.0;
			ed->args[th_id][st_id].real = 0.0;
			ed->args[th_id][st_id].imag = 0.0;
			ed->non_drv_tmp[th_id][st_id].real = 0.0;
			ed->non_drv_tmp[th_id][st_id].imag = 0.0;
			ed->drv_tmp[th_id][st_id].real = 0.0;
			ed->drv_tmp[th_id][st_id].imag = 0.0;
		}
	}
}

void init_rk(AllData * ad)
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	ed->rk_step = md->T / double(cp->rk_ns);

	init_rk_data(ad);
}

void init_rk_deep(AllData * ad)
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int num_branches = md->num_ham_qj;
	int num_sub_steps = num_branches * int(cp->params.find("deep_num_steps")->second);
	int total_steps = cp->rk_ns * num_sub_steps;
	ed->rk_step = md->T / (double(total_steps));

	init_rk_data(ad);
}

void free_rk_data(AllData * ad)
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int num_threads = rp->num_threads;
	int num_trajectories = cp->num_trajectories;
	int sys_size = md->sys_size;

	for (int th_id = 0; th_id < num_threads; th_id++)
	{
		delete[] ed->k1[th_id];
		delete[] ed->k2[th_id];
		delete[] ed->k3[th_id];
		delete[] ed->k4[th_id];
		delete[] ed->args[th_id];
		delete[] ed->non_drv_tmp[th_id];
		delete[] ed->drv_tmp[th_id];
	}
	delete[] ed->k1;
	delete[] ed->k2;
	delete[] ed->k3;
	delete[] ed->k4;
	delete[] ed->args;
	delete[] ed->non_drv_tmp;
	delete[] ed->drv_tmp;
}

void free_rk(AllData * ad)
{
	free_rk_data(ad);
}

void free_rk_deep(AllData * ad)
{
	free_rk_data(ad);
}

void set_init_args(AllData * ad, int tr_id, int th_id)
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int num_threads = rp->num_threads;
	int num_trajectories = cp->num_trajectories;
	int sys_size = md->sys_size;

	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		ed->args[th_id][st_id].real = ed->phi_all[tr_id * sys_size + st_id].real;
		ed->args[th_id][st_id].imag = ed->phi_all[tr_id * sys_size + st_id].imag;
	}
}

void arg_upd(AllData * ad, int sub_step, int tr_id, int th_id)
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;

	if (sub_step == 1)
	{
		for (int st_id = 0; st_id < sys_size; st_id++)
		{
			ed->args[th_id][st_id].real = ed->phi_all[tr_id * sys_size + st_id].real + 0.5 * ed->rk_step * ed->k1[th_id][st_id].real;
			ed->args[th_id][st_id].imag = ed->phi_all[tr_id * sys_size + st_id].imag + 0.5 * ed->rk_step * ed->k1[th_id][st_id].imag;
		}
	}
	else if (sub_step == 2)
	{
		for (int st_id = 0; st_id < sys_size; st_id++)
		{
			ed->args[th_id][st_id].real = ed->phi_all[tr_id * sys_size + st_id].real + 0.5 * ed->rk_step * ed->k2[th_id][st_id].real;
			ed->args[th_id][st_id].imag = ed->phi_all[tr_id * sys_size + st_id].imag + 0.5 * ed->rk_step * ed->k2[th_id][st_id].imag;
		}
	}
	else if (sub_step == 3)
	{
		for (int st_id = 0; st_id < sys_size; st_id++)
		{
			ed->args[th_id][st_id].real = ed->phi_all[tr_id * sys_size + st_id].real + 1.0 * ed->rk_step * ed->k3[th_id][st_id].real;
			ed->args[th_id][st_id].imag = ed->phi_all[tr_id * sys_size + st_id].imag + 1.0 * ed->rk_step * ed->k3[th_id][st_id].imag;
		}
	}
}

void rk_final(AllData * ad, int tr_id, int th_id, double step)
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;

	MKL_Complex16 * k1 = ed->k1[th_id];
	MKL_Complex16 * k2 = ed->k2[th_id];
	MKL_Complex16 * k3 = ed->k3[th_id];
	MKL_Complex16 * k4 = ed->k4[th_id];

	for (int st_id = 0; st_id < md->sys_size; st_id++)
	{
		ed->phi_all[tr_id * sys_size + st_id].real += (k1[st_id].real + 2.0 * k2[st_id].real + 2.0 * k3[st_id].real + k4[st_id].real) * step / 6.0;
		ed->phi_all[tr_id * sys_size + st_id].imag += (k1[st_id].imag + 2.0 * k2[st_id].imag + 2.0 * k3[st_id].imag + k4[st_id].imag) * step / 6.0;
	}
}

void rk_recovery(AllData * ad, int tr_id, int th_id)
{
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	MKL_Complex16 ZERO = { 0.0, 0.0 };
	MKL_Complex16 ONE = { 1.0, 0.0 };

	int sys_size = md->sys_size;

	VSLStreamStatePtr * stream = &(ed->streams[tr_id]);
	MKL_Complex16 * phi = &(ed->phi_all[tr_id * sys_size]);
	double * eta = &(ed->etas_all[tr_id]);
	int k = md->num_diss;

	MKL_Complex16 * res = new MKL_Complex16[sys_size];
	double norm = sqrt(norm_square(phi, sys_size));
	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		phi[st_id].real /= (norm);
		phi[st_id].imag /= (norm);
	}
	norm = 0.0;

	double * gnorms = new double[k];
	double tmp = 0.0;
	double ran = 0.0;
	vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, *stream, 1, &ran, 0.0, 1.0);

	for (int diss_id = 0; diss_id < k; diss_id++)
	{
		cblas_zgemv(
			CblasRowMajor,
			CblasNoTrans,
			sys_size,
			sys_size,
			&ONE,
			md->dissipators[diss_id],
			sys_size,
			phi,
			1,
			&ZERO,
			res,
			1
		);

		gnorms[diss_id] = norm_square(res, sys_size);
		gnorms[diss_id] *= 1.0;
		tmp += gnorms[diss_id];
	}

	ran *= tmp;

	int index = 0;
	while (ran - gnorms[index] > 0.0)
	{
		ran -= gnorms[index];
		index++;
		if (index == k - 1)
		{
			break;
		}
	}

	while (gnorms[index] == 0)
	{
		if (index == 0)
		{
			index++;
		}
		else
		{
			index--;
		}
	}

	memset(res, 0, sys_size * sizeof(MKL_Complex16));

	cblas_zgemv(
		CblasRowMajor,
		CblasNoTrans,
		sys_size,
		sys_size,
		&ONE,
		md->dissipators[index],
		sys_size,
		phi,
		1,
		&ZERO,
		res,
		1
	);

	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		phi[st_id].real = res[st_id].real / sqrt(gnorms[index] / 1.0);
		phi[st_id].imag = res[st_id].imag / sqrt(gnorms[index] / 1.0);
	}

	delete[] res;
	delete[] gnorms;
}

void save_phi_prev(AllData * ad, int tr_id, int th_id)
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;
	double step = ed->rk_step;

	for (int st_id = 0; st_id < md->sys_size; st_id++)
	{
		ed->phi_all_aux[tr_id * sys_size + st_id].real = ed->phi_all[tr_id * sys_size + st_id].real;
		ed->phi_all_aux[tr_id * sys_size + st_id].imag = ed->phi_all[tr_id * sys_size + st_id].imag;
	}
}

void restore_from_prev(AllData * ad, int tr_id, int th_id, double step)
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;

	for (int st_id = 0; st_id < md->sys_size; st_id++)
	{
		ed->phi_all[tr_id * sys_size + st_id].real = ed->phi_all_aux[tr_id * sys_size + st_id].real;
		ed->phi_all[tr_id * sys_size + st_id].imag = ed->phi_all_aux[tr_id * sys_size + st_id].imag;
	}

	ed->times_all[tr_id] -= step;
} 
