#include "initiator.h"

void DimerInitBehaviour::init_sizes(RunParam * rp, ConfigParam * cp, MainData * md) const
{
	int N = int(cp->params.find("N")->second);
	
	md->sys_size = N + 1;
	md->num_diss = 1;
}

void DimerInitBehaviour::init_hamiltonian(RunParam * rp, ConfigParam * cp, MainData * md) const
{
	int N = int(cp->params.find("N")->second);

	double prm_E = double(cp->params.find("prm_E")->second);
	double prm_U = double(cp->params.find("prm_U")->second);
	prm_U = prm_U / double(md->sys_size - 1);
	double prm_J = double(cp->params.find("prm_J")->second);

	md->hamiltonian = new double[md->sys_size * md->sys_size];
	md->hamiltonian_drv = new double[md->sys_size * md->sys_size];

	for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
	{
		for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
		{
			int index = st_id_1 * md->sys_size + st_id_2;

			md->hamiltonian[index] = 0.0;
			md->hamiltonian_drv[index] = 0.0;
		}
	}

	for (int st_id = 0; st_id < md->sys_size; st_id++)
	{
		int index = st_id * md->sys_size + st_id;

		md->hamiltonian[index] = 2.0 * prm_U * double(st_id * (st_id - 1) + (N - (st_id + 1)) * (N - (st_id + 1) - 1));
		md->hamiltonian_drv[index] = double((N - (st_id + 1)) - st_id);
	}

	for (int st_id = 0; st_id < (md->sys_size - 1); st_id++)
	{
		int down_left = (st_id + 1) * md->sys_size + st_id;
		int up_right = st_id * md->sys_size + (st_id + 1);

		md->hamiltonian[down_left] -= prm_J * sqrt(double((md->sys_size - (st_id + 1)) * (st_id + 1)));
		md->hamiltonian[up_right] -= prm_J * sqrt(double((st_id + 1) * (md->sys_size - (st_id + 1))));
	}
}

void DimerInitBehaviour::init_dissipator(RunParam * rp, ConfigParam * cp, MainData * md) const
{
	int N = int(cp->params.find("N")->second);

	md->dissipators = new MKL_Complex16 *[md->num_diss];
	for (int diss_id = 0; diss_id < md->num_diss; diss_id++)
	{
		md->dissipators[diss_id] = new MKL_Complex16[md->sys_size * md->sys_size];
	}

	for (int st_id_1 = 0; st_id_1 < md->sys_size; st_id_1++)
	{
		for (int st_id_2 = 0; st_id_2 < md->sys_size; st_id_2++)
		{
			int index = st_id_1 * md->sys_size + st_id_2;

			for (int diss_id = 0; diss_id < md->num_diss; diss_id++)
			{
				md->dissipators[diss_id][index].real = 0.0;
				md->dissipators[diss_id][index].imag = 0.0;
			}
		}
	}

	for (int st_id = 0; st_id < md->sys_size; st_id++)
	{
		int index = st_id * md->sys_size + st_id;

		md->dissipators[0][index].real = double(st_id - (N - (st_id + 1)));
		md->dissipators[0][index].imag = 0.0;
	}

	for (int st_id = 0; st_id < (md->sys_size - 1); st_id++)
	{
		int down_left = (st_id + 1) * md->sys_size + st_id;
		int up_right = st_id * md->sys_size + (st_id + 1);

		md->dissipators[0][down_left].real -= sqrt(double((md->sys_size - (st_id + 1)) * (st_id + 1)));
		md->dissipators[0][up_right].real += sqrt(double((st_id + 1) * (md->sys_size - (st_id + 1))));
	}
}