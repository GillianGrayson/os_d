#include "qj_proc.h"
#include "experiment.h"

void init_split_branches(Split * branch, int branch_id, AllData * ad)
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;

	Split * node = branch;

	int deep = cp->ex_deep;
	int	sys_size = (node->prev)->N;

	double T = (node->prev)->dt;
	double dt = T * 0.5;

	int mtx_syze = sys_size * sys_size;

	for (int i = 0; i < deep - 1; i++)
	{
		node->steps = 2;
		dt *= 0.5;
		node->dt = dt;
		node->counter = 0;
		node->g = 0;
		node->N = sys_size;
		node->type = true;

		Eigen::MatrixXcd ham_qj_eigen(sys_size, sys_size);
		for (int st_id_1 = 0; st_id_1 < sys_size; st_id_1++)
		{
			for (int st_id_2 = 0; st_id_2 < sys_size; st_id_2++)
			{
				int index = st_id_1 * sys_size + st_id_2;

				std::complex<double> val;
				val.real(md->hamiltonians_qj[branch_id][index].real * (node->dt));
				val.imag(md->hamiltonians_qj[branch_id][index].imag * (node->dt));
				ham_qj_eigen(st_id_1, st_id_2) = val;
			}
		}

		Eigen::MatrixXcd matrix_eigen = ham_qj_eigen.exp();
		cout << "exp() branch_id: " << branch_id << " i:" << i << endl;

		node->matrix = new MKL_Complex16[mtx_syze];
		for (int st_id_1 = 0; st_id_1 < sys_size; st_id_1++)
		{
			for (int st_id_2 = 0; st_id_2 < sys_size; st_id_2++)
			{
				int index = st_id_1 * sys_size + st_id_2;

				std::complex<double> val = matrix_eigen(st_id_1, st_id_2);
				node->matrix[index].real = val.real();
				node->matrix[index].imag = val.imag();
			}
		}

		node->next = new Split;
		node->next->prev = node;
		node = node->next;
	}

	node->steps = 2;
	dt *= 0.5;
	node->dt = dt;
	node->counter = 0;
	node->g = 0;
	node->N = sys_size;
	node->type = true;

	Eigen::MatrixXcd ham_qj_eigen(sys_size, sys_size);
	for (int st_id_1 = 0; st_id_1 < sys_size; st_id_1++)
	{
		for (int st_id_2 = 0; st_id_2 < sys_size; st_id_2++)
		{
			int index = st_id_1 * sys_size + st_id_2;

			std::complex<double> val;
			val.real(md->hamiltonians_qj[branch_id][index].real * (node->dt));
			val.imag(md->hamiltonians_qj[branch_id][index].imag * (node->dt));
			ham_qj_eigen(st_id_1, st_id_2) = val;
		}
	}

	Eigen::MatrixXcd matrix_eigen = ham_qj_eigen.exp();

	node->matrix = new MKL_Complex16[mtx_syze];
	for (int st_id_1 = 0; st_id_1 < sys_size; st_id_1++)
	{
		for (int st_id_2 = 0; st_id_2 < sys_size; st_id_2++)
		{
			int index = st_id_1 * sys_size + st_id_2;

			std::complex<double> val = matrix_eigen(st_id_1, st_id_2);
			node->matrix[index].real = val.real();
			node->matrix[index].imag = val.imag();
		}
	}

	node->next = 0;

	if (rp->is_debug)
	{
		string fn = rp->path + "exp_mtx_" + to_string(branch_id) + cp->fn_suffix;
		save_complex_data(fn, node->matrix, md->sys_size * md->sys_size, 16, 0);
	}
}

void copy_branch_not_member(Split * src, Split * dst)
{
	Split * node1 = src;
	Split * node2 = dst;
	while (node1->next)
	{
		node2->dt = node1->dt;
		node2->steps = node1->steps;
		node2->counter = node1->counter;
		node2->g = node1->g;
		node2->N = node1->N;
		node2->type = node1->type;
		node2->matrix = node1->matrix;
		node2->next = new Split;
		node2->next->prev = node2;
		node2 = node2->next;
		node1 = node1->next;
	}
	node2->dt = node1->dt;
	node2->steps = node1->steps;
	node2->counter = node1->counter;
	node2->g = node1->g;
	node2->N = node1->N;
	node2->type = node1->type;
	node2->matrix = node1->matrix;
	node2->next = 0;
}

void copy_struct_not_member(Split * src, Split * dst)
{
	dst->prev = src->prev;
	dst->type = src->type;
	dst->dt = src->dt;
	dst->counter = src->counter;
	dst->N = src->N;
	dst->next = new Split[dst->counter];
	for (int i = 0; i < dst->counter; i++)
	{
		(dst->next)[i].prev = dst;
		copy_branch_not_member(&((src->next)[i]), &((dst->next)[i]));
	}

	dst->steps = src->steps;
	dst->matrix = src->matrix;
	dst->g = src->g;
}

void delete_branch(Split * branch)
{
	delete (branch->matrix);
	branch->prev = 0;
	branch->matrix = 0;
	branch->steps = 0;
	branch->type = false;
	branch->dt = 0;
	branch->counter = 0;
	branch->N = 0;

	branch = branch->next;
	while (branch->next->next == 0)
	{
		delete (branch->matrix);
		branch->prev = 0;
		branch->matrix = 0;
		branch->steps = 0;
		branch->type = false;
		branch->dt = 0;
		branch->counter = 0;
		branch->N = 0;

		branch = branch->next;
		delete (branch->prev);
	}
	delete (branch->matrix);
	branch->prev = 0;
	branch->matrix = 0;
	branch->steps = 0;
	branch->type = false;
	branch->dt = 0;
	branch->counter = 0;
	branch->N = 0;

	delete (branch);
}

void delete_split_struct(Split * head)
{
	for (int i = 0; i < head->counter; i++)
	{
		delete_branch(&(head->next)[i]);
	}
	delete (head->next);
	delete (head->matrix);
	delete (head->g);
}

void delete_branch_not_member(Split * branch)
{
	branch->prev = 0;
	branch->matrix = 0;
	branch->steps = 0;
	branch->type = false;
	branch->dt = 0;
	branch->counter = 0;
	branch->N = 0;

	branch = branch->next;
	while (branch->next->next == 0)
	{
		branch->prev = 0;
		branch->matrix = 0;
		branch->steps = 0;
		branch->type = false;
		branch->dt = 0;
		branch->counter = 0;
		branch->N = 0;

		branch = branch->next;
		delete (branch->prev);
	}
	branch->prev = 0;
	branch->matrix = 0;
	branch->steps = 0;
	branch->type = false;
	branch->dt = 0;
	branch->counter = 0;
	branch->N = 0;

	delete (branch);
}

void delete_split_struct_not_member(Split * head)
{
	for (int i = 0; i < head->counter; i++)
	{
		delete_branch_not_member(&(head->next)[i]);
	}
	delete (head->next);
	head->matrix = 0;
	head->g = 0;
}

void prop_step(MKL_Complex16 * phi, MKL_Complex16 * matrix, MKL_Complex16 * res, int sys_size)
{
	MKL_Complex16 ZERO = { 0.0, 0.0 };
	MKL_Complex16 ONE = { 1.0, 0.0 };
	cblas_zgemv(CblasRowMajor, CblasNoTrans, sys_size, sys_size, &ONE, matrix, sys_size, phi, 1, &ZERO, res, 1);
}

void one_period_branch(AllData * ad, Split * head, int tr_id, Split * branch)
{
	MainData * md = ad->md;
	ExpData * ed = ad->ed;
	ConfigParam * cp = ad->cp;

	int sys_size = md->sys_size;

	int jump = int(cp->params.find("jump")->second);

	VSLStreamStatePtr * stream = &(ed->streams[tr_id]);
	MKL_Complex16 * phi = &(ed->phi_all[tr_id * sys_size]);
	MKL_Complex16 * phi_aux = &(ed->phi_all_aux[tr_id * sys_size]);
	double * eta = &(ed->etas_all[tr_id]);
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
				if (jump > 0 && ed->is_obs == 1)
				{
					double jump_time = ed->times_all[tr_id] + branch->dt;
					double jump_norm = norm_square(phi_aux, branch->N);
					double jump_eta = *eta;

					ed->jump_times[tr_id].push_back(jump_time);
					ed->jump_norms[tr_id].push_back(jump_norm);
					ed->jump_etas[tr_id].push_back(jump_eta);

					ed->jumps_counts[tr_id]++;
				}

				recovery(ad, head, tr_id);
				vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, *stream, 1, eta, 0.0, 1.0);
				while (*eta == 0.0)
				{
					vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, *stream, 1, eta, 0.0, 1.0);
				}
			}

			memcpy(phi, phi_aux, sizeof(MKL_Complex16) * branch->N);
			branch->counter++;

			ed->times_all[tr_id] += branch->dt;
		}
	}
	else
	{
		while (branch->counter != branch->steps)
		{
			prop_step(phi, branch->matrix, phi_aux, branch->N);
			if (is_norm_crossed(phi_aux, eta, branch->N))
			{
				one_period_branch(ad, head, tr_id, branch->next);
				ed->times_all[tr_id] -= branch->dt;
			}
			else
			{
				memcpy(phi, phi_aux, sizeof(MKL_Complex16) * branch->N);
			}
			branch->counter++;
			ed->times_all[tr_id] += branch->dt;
		}
	}

	branch->counter = 0;
}

void one_sub_period_deep(AllData * ad, int tr_id, int part_id, int thread_id)
{
	RunParam * rp = ad->rp;
	MainData * md = ad->md;

	int num_threads = rp->num_threads;

	int split_id = num_threads * part_id + thread_id;

	Split * head = &(md->splits)[split_id];

	for (int b_id = 0; b_id < head->counter; b_id++)
	{
		one_period_branch(ad, head, tr_id, &(head->next)[b_id]);
	}
}