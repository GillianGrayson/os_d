#include "split_proc.h"

Split * init_split_structure_cd(RunParam * rp, ConfigParam * cp, MainData * md)
{
	int cd_num_sub_steps = int(cp->params.find("cd_num_sub_steps")->second);

	double T = md->T / double(cd_num_sub_steps);
	int N = md->sys_size;
	int num_branches = md->num_ham_qj;

	Split * head = new Split[num_branches];

	for (int br_id = 0; br_id < num_branches; br_id++)
	{
		Split * branch = (&head[br_id]);
		branch->prev = 0;
		branch->type = false;
		branch->dt = T;
		branch->counter = 2;
		branch->N = N;
		branch->next = new Split[2];

		for (unsigned int i = 0; i < 2; i++)
		{
			(branch->next)[i].prev = branch;
			init_split_branches(&((branch->next)[i]), br_id, rp, cp, md);
		}

		branch->steps = md->num_diss;

		branch->matrix = new MKL_Complex16[branch->steps * N * N];
		branch->g = new double[branch->steps];

		for (int diss_id = 0; diss_id < branch->steps; diss_id++)
		{
			branch->g[diss_id] = 1.0;

			for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
			{
				for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
				{
					int index_xtd = diss_id * (md->sys_size * md->sys_size) + st_id_1 * md->sys_size + st_id_2;
					int index = st_id_1 * md->sys_size + st_id_2;

					branch->matrix[index_xtd].real = md->dissipators[diss_id][index].real;
					branch->matrix[index_xtd].imag = md->dissipators[diss_id][index].imag;
				}
			}
		}
	}

	return head;
}

Split * init_split_structure(RunParam * rp, ConfigParam * cp, MainData * md)
{
	Split * head = new Split[1];

	double T = md->T;
	int N = md->sys_size;
	int num_branches = md->num_ham_qj;

	head->prev = 0;
	head->type = false;
	head->dt = T;
	head->counter = num_branches;
	head->N = N;
	head->next = new Split[num_branches];

	for (unsigned int i = 0; i < num_branches; i++)
	{
		(head->next)[i].prev = head;
		init_split_branches(&((head->next)[i]), i, rp, cp, md);
	}

	head->steps = md->num_diss;

	head->matrix = new MKL_Complex16[head->steps * N * N];
	head->g = new double[head->steps];

	for (int diss_id = 0; diss_id < head->steps; diss_id++)
	{
		head->g[diss_id] = 1.0;

		for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
		{
			for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
			{
				int index_xtd = diss_id * (md->sys_size * md->sys_size) + st_id_1 * md->sys_size + st_id_2;
				int index = st_id_1 * md->sys_size + st_id_2;

				head->matrix[index_xtd].real = md->dissipators[diss_id][index].real;
				head->matrix[index_xtd].imag = md->dissipators[diss_id][index].imag;
			}
		}
	}

	return head;
}

void init_split_branches(Split * branch, int branch_id, RunParam * rp, ConfigParam * cp, MainData * md)
{
	Split * node = branch;

	int deep = cp->qj_deep;
	int	sys_size = (node->prev)->N;

	double T = (node->prev)->dt;
	double dt = T * 0.5;

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

		node->matrix = new MKL_Complex16[sys_size * sys_size];
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

	node->matrix = new MKL_Complex16[sys_size * sys_size];
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
	for (unsigned int i = 0; i < dst->counter; i++)
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
	for (unsigned int i = 0; i < head->counter; i++)
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
	for (unsigned int i = 0; i < head->counter; i++)
	{
		delete_branch_not_member(&(head->next)[i]);
	}
	delete (head->next);
	head->matrix = 0;
	head->g = 0;
}