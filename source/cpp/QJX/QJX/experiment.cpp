#include "experiment.h"

void LpnExperimentBehaviour::trans_process(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const
{
	ConfigParam * cp = ad->cp;

	int num_trajectories = cp->num_trajectories;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);
	int dump_evo_avg = int(cp->params.find("dump_evo_avg")->second);

#pragma omp parallel for
	for (int tr_id = 0; tr_id < 1; tr_id++)
	{
		int thread_id = omp_get_thread_num();
		trans_process_single(ad, pb, cb, tr_id, thread_id);
	}

#pragma omp parallel for
	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		if (tr_id > 0)
		{
			copy_trajectory_lpn(ad, tr_id);
			var_trajectory_lpn(ad, tr_id);
		}

		resresh_times(ad, tr_id);

		cb->calc_chars_std_start(ad, tr_id);
		cb->calc_chars_lpn_start(ad, tr_id);

		cb->evo_chars_std(ad, tr_id, 0);
		cb->evo_chars_lpn(ad, tr_id, 0);

		if (dump_evo_sep == 1)
		{
			dump_adr_single(ad, tr_id, false);
		}
	}

	if (dump_evo_avg == 1)
	{
		dump_adr_avg(ad, false);
	}
}

void LpnExperimentBehaviour::obser_process(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const
{	
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	ExpData * ed = ad->ed;

	ed->is_obs = 1;

	int num_trajectories = cp->num_trajectories;

	int dump_num_total = ed->dump_num_total;
	int * dump_periods = ed->dump_periods;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);
	int dump_evo_avg = int(cp->params.find("dump_evo_avg")->second);

	int begin_period_id = 0;
	int end_period_id = 0;
	for (int dump_id = 1; dump_id < dump_num_total; dump_id++)
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
				pb->one_period(ad, cb, tr_id, thread_id, period_id);
				cb->calc_chars_std(ad, tr_id);
			}

			cb->calc_chars_lpn(ad, 0);

			ed->period_id = (period_id + 1);

#pragma omp parallel for
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				if (tr_id > 0)
				{
					lambda_lpn(ad, cb, tr_id);
				}
			}
		}

#pragma omp parallel for
		for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
		{
			cb->evo_chars_std(ad, tr_id, dump_id);
			cb->evo_chars_lpn(ad, tr_id, dump_id);

			if (dump_evo_sep == 1)
			{
				dump_adr_single(ad, tr_id, true);
			}
		}

		if (dump_evo_avg == 1)
		{
			dump_adr_avg(ad, true);
		}
	}

	cb->dump_std(ad);
	cb->dump_lpn(ad);

	if (dump_evo_sep == 1)
	{
		cb->dump_std_evo(ad);
		cb->dump_lpn_evo(ad);
	}

	if (dump_evo_avg == 0)
	{
		dump_adr_avg(ad, true);
	}
}

void StdExperimentBehaviour::trans_process(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const
{
	ConfigParam * cp = ad->cp;

	int num_trajectories = cp->num_trajectories;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);
	int dump_evo_avg = int(cp->params.find("dump_evo_avg")->second);

#pragma omp parallel for
	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		int thread_id = omp_get_thread_num();
		trans_process_single(ad, pb, cb, tr_id, thread_id);
	}

#pragma omp parallel for
	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		resresh_times(ad, tr_id);

		cb->calc_chars_std_start(ad, tr_id);

		cb->evo_chars_std(ad, tr_id, 0);

		if (dump_evo_sep == 1)
		{
			dump_adr_single(ad, tr_id, false);
		}
	}

	if (dump_evo_avg == 1)
	{
		dump_adr_avg(ad, false);
	}
}

void StdExperimentBehaviour::obser_process(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	ExpData * ed = ad->ed;

	ed->is_obs = 1;

	int num_trajectories = cp->num_trajectories;

	int dump_num_total = ed->dump_num_total;
	int * dump_periods = ed->dump_periods;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);
	int dump_evo_avg = int(cp->params.find("dump_evo_avg")->second);

	int begin_period_id = 0;
	int end_period_id = 0;
	for (int dump_id = 1; dump_id < dump_num_total; dump_id++)
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
				pb->one_period(ad, cb, tr_id, thread_id, period_id);
			}
		}

#pragma omp parallel for
		for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
		{
			cb->calc_chars_std(ad, tr_id);

			cb->evo_chars_std(ad, tr_id, dump_id);

			if (dump_evo_sep == 1)
			{
				dump_adr_single(ad, tr_id, true);
			}
		}

		if (dump_evo_avg == 1)
		{
			dump_adr_avg(ad, true);
		}
	}

	cb->dump_std(ad);

	if (dump_evo_sep == 1)
	{
		cb->dump_std_evo(ad);
	}

	if (dump_evo_avg == 0)
	{
		dump_adr_avg(ad, false);
	}
}

void CorrDimExperimentBehaviour::trans_process(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const
{
	ConfigParam * cp = ad->cp;

	int num_trajectories = cp->num_trajectories;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);
	int dump_evo_avg = int(cp->params.find("dump_evo_avg")->second);

#pragma omp parallel for
	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		int thread_id = omp_get_thread_num();
		trans_process_single_deep(ad, pb, cb, tr_id, thread_id);
	}

#pragma omp parallel for
	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		resresh_times(ad, tr_id);

		cb->calc_chars_std_start(ad, tr_id);

		cb->evo_chars_std(ad, tr_id, 0);

		if (dump_evo_sep == 1)
		{
			dump_adr_single(ad, tr_id, false);
		}
	}

	if (dump_evo_avg == 1)
	{
		dump_adr_avg(ad, false);
	}
}

void CorrDimExperimentBehaviour::obser_process(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	ExpData * ed = ad->ed;

	ed->is_obs = 1;

	int deep_dump = int(cp->params.find("deep_dump")->second);

	int num_trajectories = cp->num_trajectories;

	int dump_num_total = ed->dump_num_total;
	if (deep_dump == 1)
	{
		dump_num_total = cp->num_obs_periods + 1;
	}

	int * dump_periods = ed->dump_periods;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);
	int dump_evo_avg = int(cp->params.find("dump_evo_avg")->second);

	int begin_period_id = 0;
	int end_period_id = 0;
	for (int dump_id = 1; dump_id < dump_num_total; dump_id++)
	{
		if (rp->is_pp == 1)
		{
			cout << "dump_id: " << dump_id << endl;
		}

		if (deep_dump == 1)
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
				pb->one_period_obs_deep_cd(ad, cb, tr_id, thread_id, period_id);
			}
		}

		if (deep_dump == 0)
		{
#pragma omp parallel for
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				cb->evo_chars_std(ad, tr_id, dump_id);

				if (dump_evo_sep == 1)
				{
					dump_adr_single(ad, tr_id, true);
				}
			}
		}

		if (dump_evo_avg == 1)
		{
			dump_adr_avg(ad, true);
		}
	}

#pragma omp parallel for
	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		cb->calc_ci(ad, tr_id);
	}

	cb->dump_std(ad);

	dump_cd(ad);

	if (dump_evo_sep == 1)
	{
		cb->dump_std_evo(ad);
	}

	if (dump_evo_avg == 0)
	{
		dump_adr_avg(ad, false);
	}
}

void SigmaExperimentBehaviour::trans_process(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const
{
	ConfigParam * cp = ad->cp;

	int num_trajectories = cp->num_trajectories;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);
	int dump_evo_avg = int(cp->params.find("dump_evo_avg")->second);

#pragma omp parallel for
	for (int tr_id = 0; tr_id < 1; tr_id++)
	{
		int thread_id = omp_get_thread_num();
		trans_process_single_deep(ad, pb, cb, tr_id, thread_id);
	}

#pragma omp parallel for
	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		if (tr_id > 0)
		{
			copy_trajectory_data(ad, tr_id);
		}

		resresh_times(ad, tr_id);

		cb->calc_chars_std_start(ad, tr_id);

		cb->evo_chars_std(ad, tr_id, 0);

		if (dump_evo_sep == 1)
		{
			dump_adr_single(ad, tr_id, false);
		}
	}

	if (dump_evo_avg == 1)
	{
		dump_adr_avg(ad, false);
	}
}

void SigmaExperimentBehaviour::obser_process(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	ExpData * ed = ad->ed;

	ed->is_obs = 1;

	int deep_dump = int(cp->params.find("deep_dump")->second);

	int num_trajectories = cp->num_trajectories;

	int dump_num_total = ed->dump_num_total;
	if (deep_dump == 1)
	{
		dump_num_total = cp->num_obs_periods + 1;
	}

	int * dump_periods = ed->dump_periods;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);
	int dump_evo_avg = int(cp->params.find("dump_evo_avg")->second);

	int begin_period_id = 0;
	int end_period_id = 0;
	for (int dump_id = 1; dump_id < dump_num_total; dump_id++)
	{
		if (rp->is_pp == 1)
		{
			cout << "dump_id: " << dump_id << endl;
		}

		if (deep_dump == 1)
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
				pb->one_period_obs_deep_sigma(ad, cb, tr_id, thread_id, period_id);
			}
		}

		if (deep_dump == 0)
		{
#pragma omp parallel for
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				cb->evo_chars_std(ad, tr_id, dump_id);

				if (dump_evo_sep == 1)
				{
					dump_adr_single(ad, tr_id, true);
				}
			}
		}

		if (dump_evo_avg == 1)
		{
			dump_adr_avg(ad, true);
		}
	}

	cb->dump_std(ad);

	if (dump_evo_sep == 1)
	{
		cb->dump_std_evo(ad);
	}

	if (dump_evo_avg == 0)
	{
		dump_adr_avg(ad, false);
	}
}

void StdDeepExperimentBehaviour::trans_process(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const
{
	ConfigParam * cp = ad->cp;

	int num_trajectories = cp->num_trajectories;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);
	int dump_evo_avg = int(cp->params.find("dump_evo_avg")->second);

#pragma omp parallel for
	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		int thread_id = omp_get_thread_num();
		trans_process_single_deep(ad, pb, cb, tr_id, thread_id);
	}

#pragma omp parallel for
	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		resresh_times(ad, tr_id);

		cb->calc_chars_std_start(ad, tr_id);

		cb->evo_chars_std(ad, tr_id, 0);

		if (dump_evo_sep == 1)
		{
			dump_adr_single(ad, tr_id, false);
		}
	}

	if (dump_evo_avg == 1)
	{
		dump_adr_avg(ad, false);
	}
}

void StdDeepExperimentBehaviour::obser_process(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	ExpData * ed = ad->ed;

	ed->is_obs = 1;

	int num_trajectories = cp->num_trajectories;

	int dump_num_total = cp->num_obs_periods + 1;

	int * dump_periods = ed->dump_periods;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);
	int dump_evo_avg = int(cp->params.find("dump_evo_avg")->second);

	int begin_period_id = 0;
	int end_period_id = 0;
	for (int dump_id = 1; dump_id < dump_num_total; dump_id++)
	{
		if (rp->is_pp == 1)
		{
			cout << "dump_id: " << dump_id << endl;
		}

		begin_period_id = dump_id - 1;
		end_period_id = dump_id;

		for (int period_id = begin_period_id; period_id < end_period_id; period_id++)
		{
#pragma omp parallel for
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				int thread_id = omp_get_thread_num();
				pb->one_period_obs_deep(ad, cb, tr_id, thread_id, period_id);
			}
		}

		if (dump_evo_avg == 1)
		{
			dump_adr_avg(ad, true);
		}
	}

	cb->dump_std(ad);

	if (dump_evo_sep == 1)
	{
		cb->dump_std_evo(ad);
	}

	if (dump_evo_avg == 0)
	{
		dump_adr_avg(ad, false);
	}
}

void LpnDeepExperimentBehaviour::trans_process(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const
{
	ConfigParam * cp = ad->cp;

	int num_trajectories = cp->num_trajectories;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);
	int dump_evo_avg = int(cp->params.find("dump_evo_avg")->second);

#pragma omp parallel for
	for (int tr_id = 0; tr_id < 1; tr_id++)
	{
		int thread_id = omp_get_thread_num();
		trans_process_single_deep(ad, pb, cb, tr_id, thread_id);
	}

#pragma omp parallel for
	for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
	{
		if (tr_id > 0)
		{
			copy_trajectory_lpn(ad, tr_id);
			var_trajectory_lpn(ad, tr_id);
		}

		resresh_times(ad, tr_id);

		cb->calc_chars_std_start(ad, tr_id);
		cb->calc_chars_lpn_start(ad, tr_id);

		cb->evo_chars_std(ad, tr_id, 0);
		cb->evo_chars_lpn(ad, tr_id, 0);

		if (dump_evo_sep == 1)
		{
			dump_adr_single(ad, tr_id, false);
		}
	}

	if (dump_evo_avg == 1)
	{
		dump_adr_avg(ad, false);
	}
}

void LpnDeepExperimentBehaviour::obser_process(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb) const
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	ExpData * ed = ad->ed;

	ed->is_obs = 1;

	int num_trajectories = cp->num_trajectories;

	int dump_num_total = cp->num_obs_periods + 1;

	int * dump_periods = ed->dump_periods;

	int dump_evo_sep = int(cp->params.find("dump_evo_sep")->second);
	int dump_evo_avg = int(cp->params.find("dump_evo_avg")->second);

	int begin_period_id = 0;
	int end_period_id = 0;
	for (int dump_id = 1; dump_id < dump_num_total; dump_id++)
	{
		if (rp->is_pp == 1)
		{
			cout << "dump_id: " << dump_id << endl;
		}

		begin_period_id = dump_id - 1;
		end_period_id = dump_id;

		for (int period_id = begin_period_id; period_id < end_period_id; period_id++)
		{
#pragma omp parallel for
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				int thread_id = omp_get_thread_num();
				pb->one_period_obs_deep_lpn(ad, cb, tr_id, thread_id, period_id);
			}

			ed->period_id = (period_id + 1);

#pragma omp parallel for
			for (int tr_id = 0; tr_id < num_trajectories; tr_id++)
			{
				if (tr_id > 0)
				{
					lambda_lpn(ad, cb, tr_id);
				}
			}
		}

		if (dump_evo_avg == 1)
		{
			dump_adr_avg(ad, true);
		}
	}

	cb->dump_std(ad);
	cb->dump_lpn(ad);

	if (dump_evo_sep == 1)
	{
		cb->dump_std_evo(ad);
		cb->dump_lpn_evo(ad);
	}

	if (dump_evo_avg == 0)
	{
		dump_adr_avg(ad, true);
	}
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

void recovery(AllData * ad, Split * head, int tr_id)
{
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	MKL_Complex16 ZERO = { 0.0, 0.0 };
	MKL_Complex16 ONE = { 1.0, 0.0 };

	int sys_size = md->sys_size;
	VSLStreamStatePtr * stream = &(ed->streams[tr_id]);
	MKL_Complex16 * phi = &(ed->phi_all_aux[tr_id * sys_size]);
	double * eta = &(ed->etas_all[tr_id]);
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

double get_energy(AllData * ad, int tr_id)
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;

	double dimer_prm_E = double(cp->params.find("dimer_prm_E")->second);

	MKL_Complex16 * phi = &(ed->phi_all[tr_id * sys_size]);
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

			double ham_val = (hamiltonian[index] + dimer_prm_E * hamiltonian_drv[index]);

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

MKL_Complex16 get_spec(AllData * ad, int tr_id)
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;

	MKL_Complex16 * phi = &(ed->phi_all[tr_id * sys_size]);
	MKL_Complex16 * phi_normed = new MKL_Complex16[sys_size];
	MKL_Complex16 * phi_normed_conj = new MKL_Complex16[sys_size];
	double norm = sqrt(norm_square(phi, sys_size));
	MKL_Complex16 * mult_tmp = new MKL_Complex16[sys_size];
	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		mult_tmp[st_id].real = 0.0;
		mult_tmp[st_id].imag = 0.0;

		phi_normed[st_id].real = phi[st_id].real / norm;
		phi_normed[st_id].imag = phi[st_id].imag / norm;

		phi_normed_conj[st_id].real = phi_normed[st_id].real;
		phi_normed_conj[st_id].imag = -phi_normed[st_id].imag;
	}

	MKL_Complex16 ZERO = { 0.0, 0.0 };
	MKL_Complex16 ONE = { 1.0, 0.0 };
	cblas_zgemv(CblasRowMajor, CblasNoTrans, sys_size, sys_size, &ONE, md->special, sys_size, phi_normed, 1, &ZERO, mult_tmp, 1);

	MKL_Complex16 result = { 0.0, 0.0 };
	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		result.real += (phi_normed_conj[st_id].real * mult_tmp[st_id].real - phi_normed_conj[st_id].imag * mult_tmp[st_id].imag);
		result.imag += (phi_normed_conj[st_id].imag * mult_tmp[st_id].real + phi_normed_conj[st_id].real * mult_tmp[st_id].imag);
	}

	delete[] mult_tmp;
	delete[] phi_normed;
	delete[] phi_normed_conj;

	return result;
}

void resresh_times(AllData * ad, int tr_id)
{
	ExpData * ed = ad->ed;

	ed->times_all[tr_id] = 0.0;
}

void copy_trajectory_lpn(AllData * ad, int tr_id)
{
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;

	VSLStreamStatePtr * streams = ed->streams;
	MKL_Complex16 * phi = &(ed->phi_all[tr_id * sys_size]);
	MKL_Complex16 * phi_original = &(ed->phi_all[0]);
	double * etas = ed->etas_all;

	vslCopyStream(&streams[tr_id], streams[0]);

	etas[tr_id] = etas[0];

	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		phi[st_id].real = phi_original[st_id].real;
		phi[st_id].imag = phi_original[st_id].imag;
	}
}

void copy_trajectory_data(AllData * ad, int tr_id)
{
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;

	MKL_Complex16 * phi = &(ed->phi_all[tr_id * sys_size]);
	MKL_Complex16 * phi_original = &(ed->phi_all[0]);

	double * etas = ed->etas_all;

	etas[tr_id] = etas[0];

	for (int st_id = 0; st_id < sys_size; st_id++)
	{
		phi[st_id].real = phi_original[st_id].real;
		phi[st_id].imag = phi_original[st_id].imag;
	}
}

void var_trajectory_lpn(AllData * ad, int tr_id)
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;

	double lpn_eps = double(cp->params.find("lpn_eps")->second);

	VSLStreamStatePtr * streams_var = ed->streams_var;
	MKL_Complex16 * phi_original = &(ed->phi_all[0]);
	MKL_Complex16 * phi = &(ed->phi_all[tr_id * sys_size]);

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
		phi[st_id].real = phi_original[st_id].real + lpn_eps * phi_var[st_id].real;
		phi[st_id].imag = phi_original[st_id].imag + lpn_eps * phi_var[st_id].imag;
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

void lambda_lpn(AllData * ad, CoreBehavior *cb, int tr_id)
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;

	int lpn_type = double(cp->params.find("lpn_type")->second);
	double lpn_delta_up = double(cp->params.find("lpn_delta_up")->second);
	double lpn_delta_down = double(cp->params.find("lpn_delta_down")->second);

	double delta_f = cb->calc_delta(ad, tr_id);

	if ((delta_f > lpn_delta_up) || (delta_f < lpn_delta_down))
	{
		ed->lambda[tr_id] += log(delta_f / ed->delta_s[tr_id] + 1.0e-16);

		var_trajectory_lpn(ad, tr_id);

		cb->calc_chars_lpn(ad, tr_id);

		ed->delta_s[tr_id] = fabs(ed->mean_lpn[tr_id] - ed->mean_lpn[0]) / double(sys_size);

		ed->lambda_now[tr_id] = ed->lambda[tr_id] / (double(ed->period_id) * md->T);
	}
	else
	{
		cb->calc_chars_lpn(ad, tr_id);
		ed->lambda_now[tr_id] = (ed->lambda[tr_id] + log(delta_f / ed->delta_s[tr_id] + 1.0e-16)) / (double(ed->period_id) * md->T);
	}
}

void trans_process_single(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb, int tr_id, int thread_id)
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;
	int num_tp_periods = cp->num_tp_periods;

	MKL_Complex16 * phi = &(ed->phi_all[tr_id * sys_size]);
	double * eta = &(ed->etas_all[tr_id]);

	*eta = 0.0;
	while (*eta == 0.0)
	{
		vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, ed->streams[tr_id], 1, eta, 0.0, 1.0);
	}

	for (int period_id = 0; period_id < num_tp_periods; period_id++)
	{
		pb->one_period(ad, cb, tr_id, thread_id, period_id);
	}
}

void trans_process_single_deep(AllData * ad, PropagateBehavior * pb, CoreBehavior * cb, int tr_id, int thread_id)
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;
	ExpData * ed = ad->ed;

	int sys_size = md->sys_size;
	int num_tp_periods = cp->num_tp_periods;

	MKL_Complex16 * phi = &(ed->phi_all[tr_id * sys_size]);
	double * eta = &(ed->etas_all[tr_id]);

	*eta = 0.0;
	while (*eta == 0.0)
	{
		vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, ed->streams[tr_id], 1, eta, 0.0, 1.0);
	}

	for (int period_id = 0; period_id < num_tp_periods; period_id++)
	{
		pb->one_period_trp_deep(ad, cb, tr_id, thread_id, period_id);
	}
}
