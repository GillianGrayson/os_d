#include "debugger.h"

void DimerDebugBehaviour::save(AllData * ad) const
{
	save_hamiltonian_and_dissipation(ad);
}

void JCSDebugBehaviour::save(AllData * ad) const
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;

	save_hamiltonian_and_dissipation(ad);

	if (rp->is_debug)
	{
		string fn = rp->path + "spec_mtx" + cp->fn_suffix;
		save_complex_data(fn, md->special, md->sys_size * md->sys_size, 16, 0);
	}
}

void PSDebugBehaviour::save(AllData * ad) const
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;

	save_hamiltonian_and_dissipation(ad);

	if (rp->is_debug)
	{
		string fn = rp->path + "spec_mtx" + cp->fn_suffix;
		save_complex_data(fn, md->special, md->sys_size * md->sys_size, 16, 0);

		fn = rp->path + "spec_2_mtx" + cp->fn_suffix;
		save_complex_data(fn, md->special_2, md->sys_size * md->sys_size, 16, 0);

		fn = rp->path + "spec_3_mtx" + cp->fn_suffix;
		save_complex_data(fn, md->special_3, md->sys_size * md->sys_size, 16, 0);

		fn = rp->path + "spec_4_mtx" + cp->fn_suffix;
		save_complex_data(fn, md->special_4, md->sys_size * md->sys_size, 16, 0);
	}
}

void MBLDebugBehaviour::save(AllData * ad) const
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;

	save_hamiltonian_and_dissipation(ad);

	if (rp->is_debug)
	{
		string fn = rp->path + "spec_mtx" + cp->fn_suffix;
		save_complex_data(fn, md->special, md->sys_size * md->sys_size, 16, 0);
	}
}

void LndHamDebugBehaviour::save(AllData* ad) const
{
	RunParam* rp = ad->rp;
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;

	save_hamiltonian_and_dissipation(ad, false);

	if (rp->is_debug)
	{
		string fn = rp->path + "spec_mtx" + cp->fn_suffix;
		save_complex_data(fn, md->special, md->sys_size * md->sys_size, 16, 0);
	}
}

void IntegrableDebugBehaviour::save(AllData* ad) const
{
	RunParam* rp = ad->rp;
	ConfigParam* cp = ad->cp;
	MainData* md = ad->md;

	save_hamiltonian_and_dissipation(ad, true);

	if (rp->is_debug)
	{
		string fn = rp->path + "spec_mtx" + cp->fn_suffix;
		save_complex_data(fn, md->special, md->sys_size * md->sys_size, 16, 0);
	}
}

void save_hamiltonian_and_dissipation(AllData * ad, bool save_diss)
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;

	if (rp->is_debug)
	{
		string hamiltonian_fn = rp->path + "hamiltonian" + cp->fn_suffix;
		save_complex_data(hamiltonian_fn, md->hamiltonian, md->sys_size * md->sys_size, 16, 0);

		if (save_diss)
		{
			for (int diss_id = 0; diss_id < md->num_diss; diss_id++)
			{
				string dissipator_fn = rp->path + "dissipator_" + to_string(diss_id) + cp->fn_suffix;
				save_complex_data(dissipator_fn, md->dissipators[diss_id], md->sys_size * md->sys_size, 16, 0);
			}
		}

		for (int qj_ham_id = 0; qj_ham_id < md->num_ham_qj; qj_ham_id++)
		{
			string dissipator_fn = rp->path + "hamiltonians_qj_" + to_string(qj_ham_id) + cp->fn_suffix;
			save_complex_data(dissipator_fn, md->hamiltonians_qj[qj_ham_id], md->sys_size * md->sys_size, 16, 0);
		}
	}
}