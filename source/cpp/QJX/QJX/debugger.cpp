#include "debugger.h"

void DimerDebugBehaviour::save(AllData * ad) const
{
	save_hamiltonian_and_dissipation(ad);
}

void JCSBehaviour::save(AllData * ad) const
{
	save_hamiltonian_and_dissipation(ad);
}

void save_hamiltonian_and_dissipation(AllData * ad)
{
	RunParam * rp = ad->rp;
	ConfigParam * cp = ad->cp;
	MainData * md = ad->md;

	if (rp->is_debug)
	{
		string hamiltonian_fn = rp->path + "hamiltonian" + cp->fn_suffix;
		save_double_data(hamiltonian_fn, md->hamiltonian, md->sys_size * md->sys_size, 16, 0);

		string hamiltonian_drv_fn = rp->path + "hamiltonian_drv" + cp->fn_suffix;
		save_double_data(hamiltonian_drv_fn, md->hamiltonian_drv, md->sys_size * md->sys_size, 16, 0);

		string dissipator_fn = rp->path + "dissipator_drv" + cp->fn_suffix;
		save_complex_data(dissipator_fn, md->dissipators[0], md->sys_size * md->sys_size, 16, 0);

		for (int qj_ham_id = 0; qj_ham_id < md->num_ham_qj; qj_ham_id++)
		{
			string dissipator_fn = rp->path + "hamiltonians_qj_" + to_string(qj_ham_id) + cp->fn_suffix;
			save_complex_data(dissipator_fn, md->hamiltonians_qj[qj_ham_id], md->sys_size * md->sys_size, 16, 0);
		}
	}
}