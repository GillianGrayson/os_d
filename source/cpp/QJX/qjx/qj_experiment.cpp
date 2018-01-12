#include "qj_experiment.h"

void LyapunovMCExperimentBehaviour::trans_process(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd) const
{
	int num_trajectories = cp->qj_num_trajectories;

#pragma omp parallel for
	for (int tr_id = 0; tr_id < 1; tr_id++)
	{
		int thread_id = omp_get_thread_num();

		tp_single_std(rp, cp, md, qjd, tr_id, thread_id);
	}

	resresh_times(rp, cp, md, qjd);
	copy_trajectories_lpn(rp, cp, md, qjd);
	var_trajectories_lpn(rp, cp, md, qjd);

#pragma omp parallel for
	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		calc_chars_start_std(rp, cp, md, qjd, tr_id);
	}

	

	// dump_adr_single(rp, cp, md, qjd, tr_id);
	// dump_adr_avg(rp, cp, md, qjd);
}

inline void QJ_step(MKL_Complex16 * phi, MKL_Complex16 * matrix, MKL_Complex16 * res, int sys_size)
{
	MKL_Complex16 ZERO = { 0.0, 0.0 };
	MKL_Complex16 ONE = { 1.0, 0.0 };
	cblas_zgemv(CblasRowMajor, CblasNoTrans, sys_size, sys_size, &ONE, matrix, sys_size, phi, 1, &ZERO, res, 1);
}

inline MKL_Complex16 mult_scalar_double(MKL_Complex16 a, double b)
{
	MKL_Complex16 res = { 0.0, 0.0 };
	res.real = a.real * b;
	res.imag = a.imag * b;
	return res;
}

inline MKL_Complex16 mult_scalar_complex(MKL_Complex16 * a, MKL_Complex16 * b, int N)
{
	MKL_Complex16 res = { 0.0, 0.0 };
	for (int i = 0; i < N; i++)
	{
		res.real += a[i].real * b[i].real + a[i].imag * b[i].imag;
		res.imag += b[i].real * a[i].imag - a[i].real * b[i].imag;
	}
	return res;
}

int is_norm_crossed(MKL_Complex16 * phi, double * eta, int sys_size)
{
	double norm = 0.0;
	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		norm += phi[st_id].real * phi[st_id].real + phi[st_id].imag * phi[st_id].imag;
	}

	if (norm > *(eta))
	{
		return 0;
	}
	else
	{
		return 1;
	}
}

double norm_square(MKL_Complex16 * phi, int sys_size)
{
	double norm = 0.0;
	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		norm += phi[st_id].real * phi[st_id].real + phi[st_id].imag * phi[st_id].imag;
	}
	return norm;
}

void recovery(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, Split * head, int tr_id)
{
	MKL_Complex16 ZERO = { 0.0, 0.0 };
	MKL_Complex16 ONE = { 1.0, 0.0 };

	int sys_size = md->sys_size;
	VSLStreamStatePtr * stream = &(qjd->streams[tr_id]);
	MKL_Complex16 * phi = &(qjd->phi_all_aux[tr_id * sys_size]);
	double * eta = &(qjd->etas_all[tr_id]);
	double * g = head->g;
	MKL_Complex16 * A = head->matrix;
	int k = head->steps;

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
	vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, &stream, 1, &ran, 0.0, 1.0);

	for (unsigned int i = 0; i < k; i++)
	{
		cblas_zgemv(
			CblasRowMajor,
			CblasNoTrans,
			sys_size,
			sys_size,
			&ONE,
			&A[i *  sys_size * sys_size + 0],
			sys_size,
			phi,
			1,
			&ZERO,
			res,
			1
		);

		gnorms[i] = norm_square(res, sys_size);
		gnorms[i] *= g[i];
		tmp += gnorms[i];
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
		&A[index *  sys_size * sys_size + 0],
		sys_size,
		phi,
		1,
		&ZERO,
		res,
		1
	);

	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		phi[st_id].real = res[st_id].real / sqrt(gnorms[index] / g[index]);
		phi[st_id].imag = res[st_id].imag / sqrt(gnorms[index] / g[index]);
	}

	delete[] res;
	delete[] gnorms;
}

void one_period_branch(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, Split * head, int tr_id, Split * branch)
{
	int sys_size = md->sys_size;

	VSLStreamStatePtr * stream = &(qjd->streams[tr_id]);
	MKL_Complex16 * phi = &(qjd->phi_all[tr_id * sys_size]);
	MKL_Complex16 * phi_aux = &(qjd->phi_all_aux[tr_id * sys_size]);
	double * eta = &(qjd->etas_all[tr_id]);
	double * g = head->g;
	MKL_Complex16 * A = head->matrix;
	int k = head->steps;

	if (branch->next == 0)
	{
		while (branch->counter != branch->steps)
		{
			QJ_step(phi, branch->matrix, phi_aux, branch->N);
			if (is_norm_crossed(phi_aux, eta, branch->N))
			{
				recovery(rp, cp, md, qjd, head, tr_id);
				vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, &stream, 1, eta, 0.0, 1.0);
				while (*eta == 0.0)
				{
					vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, &stream, 1, eta, 0.0, 1.0);
				}
			}

			memcpy(phi, phi_aux, sizeof(MKL_Complex16) * branch->N);
			branch->counter++;

			qjd->times_all[tr_id] += branch->dt;
		}
	}
	else
	{
		while (branch->counter != branch->steps)
		{
			QJ_step(phi, branch->matrix, phi_aux, branch->N);
			if (is_norm_crossed(phi_aux, eta, branch->N))
			{
				one_period_branch(rp, cp, md, qjd, head, tr_id, branch->next);
				qjd->times_all[tr_id] -= branch->dt;
			}
			else
			{
				memcpy(phi, phi_aux, sizeof(MKL_Complex16) * branch->N);
			}
			branch->counter++;
			qjd->times_all[tr_id] += branch->dt;
		}
	}

	branch->counter = 0;
}

void one_period(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int tr_id, int thread_id)
{
	int sys_size = md->sys_size;
	Split * head = &(md->splits)[thread_id];

	for (unsigned int b_id = 0; b_id < head->counter; b_id++)
	{
		one_period_branch(rp, cp, md, qjd, head, tr_id, &(head->next)[b_id]);
	}
}

double get_mean_std(double * adr, int sys_size)
{
	double mean = 0.0;
	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		mean += double(st_id) * adr[st_id];
	}

	return mean;
}

double get_dispersion_std(double mean_curr, double mean_start)
{
	double dispersion = (mean_curr - mean_start) * (mean_curr - mean_start);
	return dispersion;
}

double get_m2(double * adr, int sys_size, double mean)
{
	double m2 = 0.0;
	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		m2 += (double(st_id) - mean) * (double(st_id) - mean) * adr[st_id];
	}

	return m2;
}

double get_energy(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int tr_id)
{
	int sys_size = md->sys_size;

	double prm_E = double(cp->params.find("prm_E")->second);

	MKL_Complex16 * phi = &(qjd->phi_all[tr_id * sys_size]);
	double * hamiltonian = md->hamiltonian;
	double * hamiltonian_drv = md->hamiltonian_drv;

	double norm_2 = norm_square(phi, sys_size);

	MKL_Complex16 * sm = new MKL_Complex16[sys_size];

	for (int st_id_1 = 0; st_id_1 < sys_size; st_id_1++)
	{
		MKL_Complex16 tmp;
		tmp.real = 0;
		tmp.imag = 0;
		for (int st_id_2 = 0; st_id_2 < sys_size; st_id_2++)
		{
			int index = st_id_1 * sys_size + st_id_2;

			double ham_val = (hamiltonian[index] + prm_E * hamiltonian_drv[index]);

			tmp.real += (ham_val * phi[st_id_2].real / sqrt(norm_2) - (0.0) * phi[st_id_2].imag / sqrt(norm_2));
			tmp.imag += (ham_val * phi[st_id_2].imag / sqrt(norm_2) + (0.0) * phi[st_id_2].real / sqrt(norm_2));
		}

		sm[st_id_1].real = tmp.real;
		sm[st_id_1].imag = tmp.imag;
	}

	double energy = 0.0;
	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		energy += (phi[st_id].real / sqrt(norm_2) * sm[st_id].real + phi[st_id].imag / sqrt(norm_2) * sm[st_id].imag);
	}

	return energy;
}


void calc_chars_start_std(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int tr_id)
{
	int sys_size = md->sys_size;
	MKL_Complex16 * phi = &(qjd->phi_all[tr_id * sys_size]);
	double * adr = &(qjd->abs_diag_rho_all[tr_id * sys_size]);

	double norm_2 = norm_square(phi, sys_size);

	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		adr[st_id] = mult_scalar_double(mult_scalar_complex(&phi[st_id], &phi[st_id], 1), 1.0 / norm_2).real;
	}

	double mean = get_mean_std(adr, sys_size);
	double dispersion = get_dispersion_std(mean, mean);
	double m2 = get_m2(adr, sys_size, mean);

	qjd->mean_start[tr_id]		= mean;
	qjd->mean[tr_id]			= mean;
	qjd->dispersion[tr_id]		= dispersion;
	qjd->m2[tr_id]				= m2;
}

void calc_chars_start_lpn(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int tr_id)
{
	int sys_size = md->sys_size;
	MKL_Complex16 * phi = &(qjd->phi_all[tr_id * sys_size]);
	double * adr = &(qjd->abs_diag_rho_all[tr_id * sys_size]);

	double norm_2 = norm_square(phi, sys_size);

	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		adr[st_id] = mult_scalar_double(mult_scalar_complex(&phi[st_id], &phi[st_id], 1), 1.0 / norm_2).real;
	}

	double energy = get_energy(rp, cp, md, qjd, tr_id);
	double lambda = 0.0;
	double mean_lpn = qjd->mean[tr_id];
	double energy_lpn = energy;
	double lambda_lpn = lambda;

	qjd->energy[tr_id] = ;
	qjd->lambda[tr_id] = ;
	qjd->mean_lpn[tr_id] = ;
	qjd->energy_lpn[tr_id] = ;
	qjd->lambda_lpn[tr_id] = ;
}

void resresh_times(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd)
{
	int num_trajectories = cp->qj_num_trajectories;

	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		qjd->times_all[tr_id] = 0.0;
	}
}

void copy_trajectories_lpn(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd)
{
	int sys_size = md->sys_size;
	int num_trajectories = cp->qj_num_trajectories;

	for (int tr_id = 1; tr_id < num_trajectories; tr_id++)
	{
		VSLStreamStatePtr * streams = qjd->streams;
		MKL_Complex16 * phi = &(qjd->phi_all[tr_id * sys_size]);
		MKL_Complex16 * phi_original = &(qjd->phi_all[0]);
		double * etas = qjd->etas_all;

		vslCopyStream(&streams[tr_id], streams[0]);
		etas[tr_id] = etas[0];

		for (int st_id = 0; st_id < sys_size; st_id++)
		{
			phi[st_id].real = phi_original[st_id].real;
			phi[st_id].imag = phi_original[st_id].imag;
		}
	}
}

void var_trajectories_lpn(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd)
{
	int sys_size = md->sys_size;
	int num_trajectories = cp->qj_num_trajectories;

	double eps_lpn = double(cp->params.find("eps_lpn")->second);

	for (int tr_id = 1; tr_id < num_trajectories; tr_id++)
	{
		VSLStreamStatePtr * streams_var = qjd->streams_var;
		MKL_Complex16 * phi = &(qjd->phi_all[tr_id * sys_size]);

		MKL_Complex16 * phi_var = new MKL_Complex16[sys_size];

		double * phi_var_double = new double[2 * sys_size];
		vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, streams_var[tr_id], 2 * sys_size, phi_var_double, -1.0, 1.0);
		for (int st_id = 0; st_id < sys_size; st_id++)
		{
			phi_var[st_id].real = phi_var_double[0 * sys_size + st_id];
			phi_var[st_id].imag = phi_var_double[1 * sys_size + st_id];
		}
		delete[] phi_var_double;

		double norm_var_2 = norm_square(phi_var, sys_size);
		for (int st_id = 0; st_id < sys_size; st_id++)
		{
			phi_var[st_id].real = phi_var[st_id].real / sqrt(norm_var_2);
			phi_var[st_id].imag = phi_var[st_id].imag / sqrt(norm_var_2);
		}

		double norm_2 = norm_square(phi, sys_size);
		for (int st_id = 0; st_id < sys_size; st_id++)
		{
			phi[st_id].real += eps_lpn * phi_var[st_id].real;
			phi[st_id].imag += eps_lpn * phi_var[st_id].imag;
		}

		delete[] phi_var;

		double norm_mod_2 = norm_square(phi, sys_size);

		for (int st_id = 0; st_id < sys_size; st_id++)
		{
			phi[st_id].real /= sqrt(norm_mod_2);
			phi[st_id].imag /= sqrt(norm_mod_2);
		}

		for (int st_id = 0; st_id < sys_size; st_id++)
		{
			phi[st_id].real *= sqrt(norm_2);
			phi[st_id].imag *= sqrt(norm_2);
		}
	}
}

void tp_single_std(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int tr_id, int thread_id)
{
	int sys_size = md->sys_size;
	int num_tp_periods = cp->qj_num_tp_periods;

	VSLStreamStatePtr * stream = &(qjd->streams[tr_id]);
	MKL_Complex16 * phi = &(qjd->phi_all[tr_id * sys_size]);
	double * adr = &(qjd->abs_diag_rho_all[tr_id * sys_size]);
	double * eta = &(qjd->etas_all[tr_id]);

	*eta = 0.0;
	while (*eta == 0.0)
	{
		vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, &stream, 1, eta, 0.0, 1.0);
	}

	for (int period_id = 0; period_id < num_tp_periods; period_id++)
	{
		one_period(rp, cp, md, qjd, tr_id, thread_id);
	}
}