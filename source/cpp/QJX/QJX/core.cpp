#include "core.h"
#include "qj_proc.h"

void DimerCoreBehaviour::init_splits(AllData * ad) const
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;

	int num_threads = rp->num_threads;

	md->structure = init_split_structure_dimer(ad);
	md->splits = new Split[num_threads];
	for (int th_id = 0; th_id < num_threads; th_id++)
	{
		copy_struct_not_member(md->structure, &(md->splits)[th_id]);
	}

	cout << "Split initialization complete" << endl;
}

void DimerCoreBehaviour::free_splits(AllData * ad) const
{
	RunParam * rp = ad->rp;
	MainData * md = ad->md;

	int num_threads = rp->num_threads;

	for (int i = 0; i < num_threads; i++)
	{
		delete_split_struct_not_member(&(md->splits[i]));
	}

	delete(md->splits);
	delete_split_struct(md->structure);
}

void DimerCoreBehaviour::init_splits_deep(AllData * ad) const
{
	RunParam * rp = ad->rp;
	MainData * md = ad->md;

	int num_branches = md->num_ham_qj;
	int num_threads = rp->num_threads;

	int num_total = num_threads * num_branches;

	md->structure = init_split_structure_dimer_deep(ad);
	md->splits = new Split[num_total];

	for (int b_id = 0; b_id < num_branches; b_id++)
	{
		for (int th_id = 0; th_id < num_threads; th_id++)
		{
			int index = b_id * num_threads + th_id;
			copy_struct_not_member(&(md->structure)[b_id], &(md->splits)[index]);
		}
	}
}

void DimerCoreBehaviour::free_splits_deep(AllData * ad) const
{
	RunParam * rp = ad->rp;
	MainData * md = ad->md;

	int num_branches = md->num_ham_qj;
	int num_threads = rp->num_threads;

	for (int b_id = 0; b_id < num_branches; b_id++)
	{
		for (int th_id = 0; th_id < num_threads; th_id++)
		{
			int index = b_id * num_threads + th_id;
			delete_split_struct_not_member(&(md->splits[index]));
		}
	}

	delete(md->splits);
	delete_split_struct(md->structure);
}

void JCSCoreBehaviour::init_splits(AllData * ad) const
{
	RunParam * rp = ad->rp;
	MainData * md = ad->md;

	int num_branches = md->num_ham_qj;
	int num_threads = rp->num_threads;

	int num_total = num_threads * num_branches;

	md->structure = init_split_structure_jcs(ad);
	md->splits = new Split[num_total];

	for (int b_id = 0; b_id < num_branches; b_id++)
	{
		for (int th_id = 0; th_id < num_threads; th_id++)
		{
			int index = b_id * num_threads + th_id;
			copy_struct_not_member(&(md->structure)[b_id], &(md->splits)[index]);
		}
	}
}

void JCSCoreBehaviour::free_splits(AllData * ad) const
{
	RunParam * rp = ad->rp;
	MainData * md = ad->md;

	int num_branches = md->num_ham_qj;
	int num_threads = rp->num_threads;

	for (int b_id = 0; b_id < num_branches; b_id++)
	{
		for (int th_id = 0; th_id < num_threads; th_id++)
		{
			int index = b_id * num_threads + th_id;
			delete_split_struct_not_member(&(md->splits[index]));
		}
	}

	delete(md->splits);
	delete_split_struct(md->structure);
}

void JCSCoreBehaviour::init_splits_deep(AllData * ad) const
{
	stringstream msg;
	msg << "Error: need to add functionality" << endl;
	Error(msg.str());
}

void JCSCoreBehaviour::free_splits_deep(AllData * ad) const
{
	stringstream msg;
	msg << "Error: need to add functionality" << endl;
	Error(msg.str());
}

Split * init_split_structure_dimer(AllData * ad)
{
	MainData * md = ad->md;

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
		init_split_branches(&((head->next)[i]), i, ad);
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

Split * init_split_structure_dimer_deep(AllData * ad)
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;

	int deep_num_steps = int(cp->params.find("deep_num_steps")->second);

	double T = 0.5 * md->T / double(deep_num_steps);
	int N = md->sys_size;
	int num_branches = md->num_ham_qj;

	Split * head = new Split[num_branches];

	for (int br_id = 0; br_id < num_branches; br_id++)
	{
		Split * branch = &head[br_id];
		branch->prev = 0;
		branch->type = false;
		branch->dt = T;
		branch->counter = 2;
		branch->N = N;
		branch->next = new Split[2];

		for (unsigned int i = 0; i < 2; i++)
		{
			(branch->next)[i].prev = branch;
			init_split_branches(&((branch->next)[i]), br_id, ad);
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

Split * init_split_structure_jcs(AllData * ad)
{
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;

	int deep_num_steps = int(cp->params.find("deep_num_steps")->second);

	double jcs_drv_part_1 = double(cp->params.find("jcs_drv_part_1")->second);
	double jcs_drv_part_2 = double(cp->params.find("jcs_drv_part_2")->second);

	double T_1 = jcs_drv_part_1 * md->T;
	double T_2 = jcs_drv_part_2 * md->T;

	int N = md->sys_size;
	int num_branches = md->num_ham_qj;

	Split * head = new Split[num_branches];

	for (int br_id = 0; br_id < num_branches; br_id++)
	{
		Split * branch = &head[br_id];
		branch->prev = 0;
		branch->type = false;

		if (br_id == 0)
		{
			branch->dt = T_1;
		}
		else if (br_id == 1)
		{
			branch->dt = T_2;
		}
		else
		{
			stringstream msg;
			msg << "Error: Wrong br_id" << endl;
			Error(msg.str());
		}

		branch->counter = 2;
		branch->N = N;
		branch->next = new Split[2];

		for (unsigned int i = 0; i < 2; i++)
		{
			(branch->next)[i].prev = branch;
			init_split_branches(&((branch->next)[i]), br_id, ad);
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

Split * init_split_structure_jcs_deep(AllData * ad)
{
	stringstream msg;
	msg << "Error: need to add functionality" << endl;
	Error(msg.str());
	return nullptr;
}
