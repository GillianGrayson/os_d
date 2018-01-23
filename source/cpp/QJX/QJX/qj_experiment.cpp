#include "qj_experiment.h"

void LpnExperimentBehaviour::trans_process(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd) const
{
	int num_trajectories = cp->qj_num_trajectories;

	int is_evo_dump_sep = int(cp->params.find("is_evo_dump_sep")->second);
	int is_evo_dump_avg = int(cp->params.find("is_evo_dump_avg")->second);

#pragma omp parallel for
	for (int tr_id = 0; tr_id < 1; tr_id++)
	{
		int thread_id = omp_get_thread_num();
		trans_process_single_std(rp, cp, md, qjd, tr_id, thread_id);
	}

#pragma omp parallel for
	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		if (tr_id > 0)
		{
			copy_trajectory_lpn(rp, cp, md, qjd, tr_id);
			var_trajectory_lpn(rp, cp, md, qjd, tr_id);
		}

		resresh_times(rp, cp, md, qjd, tr_id);

		calc_chars_start_std(rp, cp, md, qjd, tr_id);
		calc_chars_start_lpn(rp, cp, md, qjd, tr_id);

		evo_chars_std(rp, cp, md, qjd, tr_id, 0);
		evo_chars_lpn(rp, cp, md, qjd, tr_id, 0);

		if (is_evo_dump_sep == 1)
		{
			dump_adr_single(rp, cp, md, qjd, tr_id, false);
		}
	}

	if (is_evo_dump_avg == 1)
	{
		dump_adr_avg(rp, cp, md, qjd, false);
	}
}

void StdExperimentBehaviour::trans_process(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd) const
{
	int num_trajectories = cp->qj_num_trajectories;

	int is_evo_dump_sep = int(cp->params.find("is_evo_dump_sep")->second);
	int is_evo_dump_avg = int(cp->params.find("is_evo_dump_avg")->second);

#pragma omp parallel for
	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		int thread_id = omp_get_thread_num();
		trans_process_single_std(rp, cp, md, qjd, tr_id, thread_id);
	}

#pragma omp parallel for
	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		resresh_times(rp, cp, md, qjd, tr_id);

		calc_chars_start_std(rp, cp, md, qjd, tr_id);

		evo_chars_std(rp, cp, md, qjd, tr_id, 0);

		if (is_evo_dump_sep == 1)
		{
			dump_adr_single(rp, cp, md, qjd, tr_id, false);
		}
	}

	if (is_evo_dump_avg == 1)
	{
		dump_adr_avg(rp, cp, md, qjd, false);
	}
}

void CorrDimExperimentBehaviour::trans_process(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd) const
{
	int num_trajectories = cp->qj_num_trajectories;

	int is_evo_dump_sep = int(cp->params.find("is_evo_dump_sep")->second);
	int is_evo_dump_avg = int(cp->params.find("is_evo_dump_avg")->second);

#pragma omp parallel for
	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		int thread_id = omp_get_thread_num();
		trans_process_single_cd(rp, cp, md, qjd, tr_id, thread_id);
	}

#pragma omp parallel for
	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		resresh_times(rp, cp, md, qjd, tr_id);

		calc_chars_start_std(rp, cp, md, qjd, tr_id);

		evo_chars_std(rp, cp, md, qjd, tr_id, 0);

		if (is_evo_dump_sep == 1)
		{
			dump_adr_single(rp, cp, md, qjd, tr_id, false);
		}
	}

	if (is_evo_dump_avg == 1)
	{
		dump_adr_avg(rp, cp, md, qjd, false);
	}
}

void LpnExperimentBehaviour::obs_process(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd) const
{
	int num_trajectories = cp->qj_num_trajectories;

	int num_dumps_total = qjd->num_dumps_total;
	int * dump_periods = qjd->dump_periods;

	int is_evo_dump_sep = int(cp->params.find("is_evo_dump_sep")->second);
	int is_evo_dump_avg = int(cp->params.find("is_evo_dump_avg")->second);

	int begin_period_id = 0;
	int end_period_id = 0;
	for (int dump_id = 1; dump_id < num_dumps_total; dump_id++)
	{
		if (rp->is_pp == 1)
		{
			cout << "dump_id: " << dump_id << endl;
		}

		begin_period_id = dump_periods[dump_id - 1];
		end_period_id = dump_periods[dump_id];

		for (int period_id = begin_period_id; period_id < end_period_id; period_id++)
		{
#pragma omp parallel for
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				int thread_id = omp_get_thread_num();
				one_period(rp, cp, md, qjd, tr_id, thread_id);
				calc_chars_std(rp, cp, md, qjd, tr_id);

				qjd->energy[tr_id] = get_energy(rp, cp, md, qjd, tr_id);
			}

			calc_chars_lpn(rp, cp, md, qjd, 0);

			qjd->period_id = (period_id + 1);

#pragma omp parallel for
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				if (tr_id > 0)
				{
					lambda_lpn(rp, cp, md, qjd, tr_id);
				}
			}
		}

#pragma omp parallel for
		for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
		{
			evo_chars_std(rp, cp, md, qjd, tr_id, dump_id);
			evo_chars_lpn(rp, cp, md, qjd, tr_id, dump_id);

			if (is_evo_dump_sep == 1)
			{
				dump_adr_single(rp, cp, md, qjd, tr_id, true);
			}
		}

		if (is_evo_dump_avg == 1)
		{
			dump_adr_avg(rp, cp, md, qjd, true);
		}
	}


	dump_std(rp, cp, md, qjd);
	dump_lpn(rp, cp, md, qjd);

	if (is_evo_dump_sep == 1)
	{
		dump_evo_std(rp, cp, md, qjd);
		dump_evo_lpn(rp, cp, md, qjd);
	}

	if (is_evo_dump_avg == 0)
	{
		dump_adr_avg(rp, cp, md, qjd, true);
	}

}

void StdExperimentBehaviour::obs_process(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd) const
{
	int num_trajectories = cp->qj_num_trajectories;

	int num_dumps_total = qjd->num_dumps_total;
	int * dump_periods = qjd->dump_periods;

	int is_evo_dump_sep = int(cp->params.find("is_evo_dump_sep")->second);
	int is_evo_dump_avg = int(cp->params.find("is_evo_dump_avg")->second);

	int begin_period_id = 0;
	int end_period_id = 0;
	for (int dump_id = 1; dump_id < num_dumps_total; dump_id++)
	{
		if (rp->is_pp == 1)
		{
			cout << "dump_id: " << dump_id << endl;
		}
		
		begin_period_id = dump_periods[dump_id - 1];
		end_period_id = dump_periods[dump_id];

		for (int period_id = begin_period_id; period_id < end_period_id; period_id++)
		{
#pragma omp parallel for
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				int thread_id = omp_get_thread_num();
				one_period(rp, cp, md, qjd, tr_id, thread_id);			
			}
		}

#pragma omp parallel for
		for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
		{
			calc_chars_std(rp, cp, md, qjd, tr_id);

			evo_chars_std(rp, cp, md, qjd, tr_id, dump_id);

			if (is_evo_dump_sep == 1)
			{
				dump_adr_single(rp, cp, md, qjd, tr_id, true);
			}
		}

		if (is_evo_dump_avg == 1)
		{
			dump_adr_avg(rp, cp, md, qjd, true);
		}
	}

	dump_std(rp, cp, md, qjd);

	if (is_evo_dump_sep == 1)
	{
		dump_evo_std(rp, cp, md, qjd);
	}

	if (is_evo_dump_avg == 0)
	{
		dump_adr_avg(rp, cp, md, qjd, false);
	}
}

void CorrDimExperimentBehaviour::obs_process(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd) const
{
	int cd_dump_deep = int(cp->params.find("cd_dump_deep")->second);

	int num_trajectories = cp->qj_num_trajectories;

	int num_dumps_total = qjd->num_dumps_total;
	if (cd_dump_deep == 1)
	{
		num_dumps_total = cp->qj_num_obs_periods + 1;
	}

	int * dump_periods = qjd->dump_periods;

	int is_evo_dump_sep = int(cp->params.find("is_evo_dump_sep")->second);
	int is_evo_dump_avg = int(cp->params.find("is_evo_dump_avg")->second);

	int begin_period_id = 0;
	int end_period_id = 0;
	for (int dump_id = 1; dump_id < num_dumps_total; dump_id++)
	{
		if (rp->is_pp == 1)
		{
			cout << "dump_id: " << dump_id << endl;
		}

		if (cd_dump_deep == 1)
		{
			begin_period_id = dump_id - 1;
			end_period_id = dump_id;
		}
		else
		{
			begin_period_id = dump_periods[dump_id - 1];
			end_period_id = dump_periods[dump_id];
		}

		for (int period_id = begin_period_id; period_id < end_period_id; period_id++)
		{
#pragma omp parallel for
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				int thread_id = omp_get_thread_num();
				one_period_cd_obs(rp, cp, md, qjd, tr_id, thread_id, period_id);
			}
		}

		if (cd_dump_deep == 0)
		{
#pragma omp parallel for
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				evo_chars_std(rp, cp, md, qjd, tr_id, dump_id);

				if (is_evo_dump_sep == 1)
				{
					dump_adr_single(rp, cp, md, qjd, tr_id, true);
				}
			}
		}

		if (is_evo_dump_avg == 1)
		{
			dump_adr_avg(rp, cp, md, qjd, true);
		}
	}

#pragma omp parallel for
	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		calc_ci(rp, cp, md, qjd, tr_id);
	}

	dump_std(rp, cp, md, qjd);

	dump_cd(rp, cp, md, qjd);

	if (is_evo_dump_sep == 1)
	{
		dump_evo_std(rp, cp, md, qjd);
	}

	if (is_evo_dump_avg == 0)
	{
		dump_adr_avg(rp, cp, md, qjd, false);
	}
}

void prop_step(MKL_Complex16 * phi, MKL_Complex16 * matrix, MKL_Complex16 * res, int sys_size)
{
	MKL_Complex16 ZERO = { 0.0, 0.0 };
	MKL_Complex16 ONE = { 1.0, 0.0 };
	cblas_zgemv(CblasRowMajor, CblasNoTrans, sys_size, sys_size, &ONE, matrix, sys_size, phi, 1, &ZERO, res, 1);
}

MKL_Complex16 mult_scalar_double(MKL_Complex16 a, double b)
{
	MKL_Complex16 res = { 0.0, 0.0 };
	res.real = a.real * b;
	res.imag = a.imag * b;
	return res;
}

MKL_Complex16 mult_scalar_complex(MKL_Complex16 * a, MKL_Complex16 * b, int N)
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
	vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, *stream, 1, &ran, 0.0, 1.0);

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

	delete(res);
	delete(gnorms);
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
			prop_step(phi, branch->matrix, phi_aux, branch->N);
			if (is_norm_crossed(phi_aux, eta, branch->N))
			{
				recovery(rp, cp, md, qjd, head, tr_id);
				vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, *stream, 1, eta, 0.0, 1.0);
				while (*eta == 0.0)
				{
					vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, *stream, 1, eta, 0.0, 1.0);
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
			prop_step(phi, branch->matrix, phi_aux, branch->N);
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

void one_sub_period_cd(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int tr_id, int part_id, int thread_id)
{
	int num_threads = rp->num_threads;

	int split_id = num_threads * part_id + thread_id;

	Split * head = &(md->splits)[split_id];

	for (unsigned int b_id = 0; b_id < head->counter; b_id++)
	{
		one_period_branch(rp, cp, md, qjd, head, tr_id, &(head->next)[b_id]);
	}
}

void one_period_cd_tp(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int tr_id, int thread_id)
{
	int num_branches = md->num_ham_qj;
	int num_sub_steps = int(cp->params.find("cd_num_sub_steps")->second);

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		for (int sub_step_id = 0; sub_step_id < num_sub_steps; sub_step_id++)
		{
			one_sub_period_cd(rp, cp, md, qjd, tr_id, part_id, thread_id);
		}
	}
}

void one_period_cd_obs(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int tr_id, int thread_id, int period_id)
{
	int cd_dump_deep = int(cp->params.find("cd_dump_deep")->second);
	int is_evo_dump_sep = int(cp->params.find("is_evo_dump_sep")->second);

	int num_branches = md->num_ham_qj;
	int num_sub_steps_per_part = int(cp->params.find("cd_num_sub_steps")->second);
	int num_sub_steps = num_branches * int(cp->params.find("cd_num_sub_steps")->second);

	int dump_point_id = 0;
	int curr_point_id = 0;
	int global_point_id = 0;

	int step_id = 0;

	for (int part_id = 0; part_id < num_branches; part_id++)
	{
		for (int sub_step_id = 0; sub_step_id < num_sub_steps_per_part; sub_step_id++)
		{
			step_id = part_id * num_sub_steps_per_part + sub_step_id;

			global_point_id = period_id * num_sub_steps + dump_point_id;

			dump_point_id++;

			for (int cd_st_id = 0; cd_st_id < qjd->cd_dim; cd_st_id++)
			{
				curr_point_id = global_point_id - cd_st_id;

				if (curr_point_id >= 0 && curr_point_id < qjd->cd_num_points)
				{
					qjd->cd_rec_data[tr_id][curr_point_id][cd_st_id] = qjd->mean[tr_id];
				}
			}

			one_sub_period_cd(rp, cp, md, qjd, tr_id, part_id, thread_id);
			calc_chars_std(rp, cp, md, qjd, tr_id);

			if (cd_dump_deep == 1)
			{
				int dump_id = global_point_id + 1;

				evo_chars_std(rp, cp, md, qjd, tr_id, dump_id);

				if (is_evo_dump_sep == 1)
				{
					dump_adr_single(rp, cp, md, qjd, tr_id, true);
				}
			}
		}
	}
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

double get_norm_cd(double * vec, int size)
{
	double sum = 0.0;
	for (int i = 0; i < size; i++)
	{
		sum += vec[i] * vec[i];
	}

	return sqrt(sum);
}

double get_mean_simple(double * adr, int sys_size)
{
	double mean = 0.0;
	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		mean += double(st_id) * adr[st_id];
	}

	return mean;
}

double get_dispersion_simple(double mean_curr, double mean_start)
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

	delete[] sm;

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

	double mean = get_mean_simple(adr, sys_size);
	double dispersion = get_dispersion_simple(mean, mean);
	double m2 = get_m2(adr, sys_size, mean);

	qjd->mean_start[tr_id]		= mean;
	qjd->mean[tr_id]			= mean;
	qjd->dispersion[tr_id]		= dispersion;
	qjd->m2[tr_id]				= m2;
}

void calc_chars_std(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int tr_id)
{
	int sys_size = md->sys_size;
	MKL_Complex16 * phi = &(qjd->phi_all[tr_id * sys_size]);
	double * adr = &(qjd->abs_diag_rho_all[tr_id * sys_size]);

	double norm_2 = norm_square(phi, sys_size);

	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		adr[st_id] = mult_scalar_double(mult_scalar_complex(&phi[st_id], &phi[st_id], 1), 1.0 / norm_2).real;
	}

	double mean = get_mean_simple(adr, sys_size);
	double dispersion = get_dispersion_simple(mean, qjd->mean_start[tr_id]);
	double m2 = get_m2(adr, sys_size, mean);

	qjd->mean[tr_id] = mean;
	qjd->dispersion[tr_id] = dispersion;
	qjd->m2[tr_id] = m2;
}

void evo_chars_std(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int tr_id, int dump_id)
{
	int num_trajectories = cp->qj_num_trajectories;
	int num_dumps_total = qjd->num_dumps_total;

	int index = tr_id * num_dumps_total + dump_id;

	qjd->mean_evo[index] = qjd->mean[tr_id];
	qjd->dispersion_evo[index] = qjd->dispersion[tr_id];
	qjd->m2_evo[index] = qjd->m2[tr_id];
}

void calc_chars_start_lpn(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int tr_id)
{
	int sys_size = md->sys_size;

	int type_lpn = double(cp->params.find("type_lpn")->second);

	MKL_Complex16 * phi = &(qjd->phi_all[tr_id * sys_size]);
	double * adr = &(qjd->abs_diag_rho_all[tr_id * sys_size]);

	double norm_2 = norm_square(phi, sys_size);

	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		adr[st_id] = mult_scalar_double(mult_scalar_complex(&phi[st_id], &phi[st_id], 1), 1.0 / norm_2).real;
	}

	double energy = get_energy(rp, cp, md, qjd, tr_id);
	double lambda = 0.0;
	double lambda_now = 0.0;
	double mean_lpn = qjd->mean[tr_id];
	double energy_lpn = energy;

	double delta_s = 0.0;
	if (type_lpn == 0)
	{
		delta_s = fabs(qjd->mean[tr_id] - qjd->mean[0]) / double(sys_size);
	}
	else if (type_lpn == 0)
	{
		delta_s = fabs(qjd->energy[tr_id] - qjd->energy[0]) / qjd->max_energy;
	}
	else 
	{
		delta_s = fabs(qjd->mean[tr_id] - qjd->mean[0]) / double(sys_size);
	}

	qjd->energy[tr_id] = energy;
	qjd->lambda[tr_id] = lambda;
	qjd->lambda_now[tr_id] = lambda_now;
	qjd->mean_lpn[tr_id] = mean_lpn;
	qjd->energy_lpn[tr_id] = energy_lpn;
	qjd->delta_s[tr_id] = delta_s;
}

void calc_chars_lpn(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int tr_id)
{
	int sys_size = md->sys_size;

	int type_lpn = double(cp->params.find("type_lpn")->second);

	MKL_Complex16 * phi = &(qjd->phi_all[tr_id * sys_size]);
	double * adr = &(qjd->abs_diag_rho_all[tr_id * sys_size]);

	double norm_2 = norm_square(phi, sys_size);

	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		adr[st_id] = mult_scalar_double(mult_scalar_complex(&phi[st_id], &phi[st_id], 1), 1.0 / norm_2).real;
	}

	double mean_lpn = get_mean_simple(adr, sys_size);
	double energy_lpn = get_energy(rp, cp, md, qjd, tr_id);

	qjd->mean_lpn[tr_id] = mean_lpn;
	qjd->energy_lpn[tr_id] = energy_lpn;
}

void calc_ci(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int tr_id)
{
	double eps = double(cp->params.find("cd_eps")->second);

	double * curr_diff = new double[qjd->cd_dim];
	double curr_norm = 0.0;
	double integral = 0.0;

	for (int cd_p_id_1 = 0; cd_p_id_1 < qjd->cd_num_points; cd_p_id_1++)
	{
		for (int cd_p_id_2 = 0; cd_p_id_2 < qjd->cd_num_points; cd_p_id_2++)
		{
			if (cd_p_id_1 != cd_p_id_2)
			{
				for (int cd_st_id = 0; cd_st_id < qjd->cd_dim; cd_st_id++)
				{
					curr_diff[cd_st_id] = qjd->cd_rec_data[tr_id][cd_p_id_1][cd_st_id] - qjd->cd_rec_data[tr_id][cd_p_id_2][cd_st_id];
				}

				curr_norm = get_norm_cd(curr_diff, qjd->cd_dim);

				if (curr_norm < eps)
				{
					integral += 1.0;
				}
			}
		}
	}


	integral /= (double(qjd->cd_num_points) * double(qjd->cd_num_points - 1));

	qjd->cd_i[tr_id] = integral;

	delete curr_diff;
}

void evo_chars_lpn(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int tr_id, int dump_id)
{
	int num_trajectories = cp->qj_num_trajectories;
	int num_dumps_total = qjd->num_dumps_total;

	int index = tr_id * num_dumps_total + dump_id;

	qjd->energy_evo[index] = qjd->energy[tr_id];
	qjd->lambda_evo[index] = qjd->lambda_now[tr_id];
	qjd->mean_lpn_evo[index] = qjd->mean_lpn[tr_id];
	qjd->energy_lpn_evo[index] = qjd->energy_lpn[tr_id];
}

void resresh_times(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int tr_id)
{
		qjd->times_all[tr_id] = 0.0;
}

void copy_trajectory_lpn(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int tr_id)
{
	int sys_size = md->sys_size;

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

void var_trajectory_lpn(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int tr_id)
{
	int sys_size = md->sys_size;

	double eps_lpn = double(cp->params.find("eps_lpn")->second);

	VSLStreamStatePtr * streams_var = qjd->streams_var;
	MKL_Complex16 * phi_original = &(qjd->phi_all[0]);
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

	double norm_2 = norm_square(phi_original, sys_size);
	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		phi[st_id].real = phi_original[st_id].real + eps_lpn * phi_var[st_id].real;
		phi[st_id].imag = phi_original[st_id].imag + eps_lpn * phi_var[st_id].imag;
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

void lambda_lpn(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int tr_id)
{
	int sys_size = md->sys_size;

	int type_lpn = double(cp->params.find("type_lpn")->second);
	double delta_up_lpn = double(cp->params.find("delta_up_lpn")->second);
	double delta_down_lpn = double(cp->params.find("delta_down_lpn")->second);

	double delta_f = 0.0;
	if (type_lpn == 0)
	{
		delta_f = fabs(qjd->mean[tr_id] - qjd->mean[0]) / double(sys_size);
	}
	else if (type_lpn == 0)
	{
		delta_f = fabs(qjd->energy[tr_id] - qjd->energy[0]) / qjd->max_energy;
	}
	else
	{
		delta_f = fabs(qjd->mean[tr_id] - qjd->mean[0]) / double(sys_size);
	}

	if ((delta_f > delta_up_lpn) || (delta_f < delta_down_lpn))
	{
		qjd->lambda[tr_id] += log(delta_f / qjd->delta_s[tr_id] + 1.0e-16);

		var_trajectory_lpn(rp, cp, md, qjd, tr_id);

		calc_chars_lpn(rp, cp, md, qjd, tr_id);

		qjd->delta_s[tr_id] = fabs(qjd->mean_lpn[tr_id] - qjd->mean_lpn[0]) / double(sys_size);

		qjd->lambda_now[tr_id] = qjd->lambda[tr_id] / (double(qjd->period_id) * md->T);
	}
	else
	{
		calc_chars_lpn(rp, cp, md, qjd, tr_id);
		qjd->lambda_now[tr_id] = (qjd->lambda[tr_id] + log(delta_f / qjd->delta_s[tr_id] + 1.0e-16)) / (double(qjd->period_id) * md->T);
	}
}

void trans_process_single_std(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int tr_id, int thread_id)
{
	int sys_size = md->sys_size;
	int num_tp_periods = cp->qj_num_tp_periods;

	MKL_Complex16 * phi = &(qjd->phi_all[tr_id * sys_size]);
	double * eta = &(qjd->etas_all[tr_id]);

	*eta = 0.0;
	while (*eta == 0.0)
	{
		vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, qjd->streams[tr_id], 1, eta, 0.0, 1.0);
	}

	for (int period_id = 0; period_id < num_tp_periods; period_id++)
	{
		one_period(rp, cp, md, qjd, tr_id, thread_id);
	}
}

void trans_process_single_cd(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int tr_id, int thread_id)
{
	int sys_size = md->sys_size;
	int num_tp_periods = cp->qj_num_tp_periods;

	MKL_Complex16 * phi = &(qjd->phi_all[tr_id * sys_size]);
	double * eta = &(qjd->etas_all[tr_id]);

	*eta = 0.0;
	while (*eta == 0.0)
	{
		vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, qjd->streams[tr_id], 1, eta, 0.0, 1.0);
	}

	for (int period_id = 0; period_id < num_tp_periods; period_id++)
	{
		one_period_cd_tp(rp, cp, md, qjd, tr_id, thread_id);
	}
}
