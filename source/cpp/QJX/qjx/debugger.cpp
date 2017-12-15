#include "debugger.h"

void DimerDebugBehaviour::save(RunParam * rp, ConfigParam * cp, MainData * md) const
{
	string hamiltonian_fn = rp->path + "hamiltonian" + cp->fn_suffix;
	save_double_data(hamiltonian_fn, md->hamiltonian, md->sys_size * md->sys_size, 16, 0);

	string hamiltonian_drv_fn = rp->path + "hamiltonian_drv" + cp->fn_suffix;
	save_double_data(hamiltonian_drv_fn, md->hamiltonian_drv, md->sys_size * md->sys_size, 16, 0);

	string dissipator_fn = rp->path + "dissipator_drv" + cp->fn_suffix;
	save_complex_data(dissipator_fn, md->dissipators[0], md->sys_size * md->sys_size, 16, 0);
}