#include "qj_experiment.h"

void LyapunovMCExperimentBehaviour::trans_process(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd) const
{
	int num_trajectories = cp->qj_num_trajectories;

#pragma omp parallel for
	for (int tr_id = 0; tr_id < 1; tr_id++)
	{
		int thread_id = omp_get_thread_num();

		tp_single_std(rp, cp, md, qjd, tr_id, thread_id);
		calc_chars_start(rp, cp, md, qjd, tr_id);
	}

	// dump_adr_single(rp, cp, md, qjd, tr_id);
	//dump_adr_avg(rp, cp, md, qjd);
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

void calc_chars_start(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int tr_id)
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

	qjd->mean_start[tr_id] = mean;
	qjd->mean[tr_id] = mean;
	qjd->dispersion[tr_id] = 0.0;
	qjd->m2[tr_id] = 0.0;
}

void resresh_times(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd)
{
	int num_trajectories = cp->qj_num_trajectories;

	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		qjd->times_all[tr_id] = 0.0;
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