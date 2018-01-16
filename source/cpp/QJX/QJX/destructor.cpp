#include "destructor.h"

void DimerFreeBehaviour::free_hamiltonians(RunParam * rp, ConfigParam * cp, MainData * md) const
{
	delete[] md->hamiltonian;
	delete[] md->hamiltonian_drv;
}

void DimerFreeBehaviour::free_dissipators(RunParam * rp, ConfigParam * cp, MainData * md) const
{
	for (int diss_id = 0; diss_id < md->num_diss; diss_id++)
	{
		delete[] md->dissipators[diss_id];
	}
	delete[] md->dissipators;
}

void DimerFreeBehaviour::free_hamiltonians_qj(RunParam * rp, ConfigParam * cp, MainData * md) const
{
	for (int qj_ham_id = 0; qj_ham_id < md->num_ham_qj; qj_ham_id++)
	{
		delete[] md->hamiltonians_qj[qj_ham_id];
	}
	delete[] md->hamiltonians_qj;
}
