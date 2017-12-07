#include "data.h"

void init_main_data_0(RunParam &rp, ConfigParam &cp, MainData &md)
{
	if (rp.sys_id == 0)
	{
		int N = int(cp.params.find("N")->second);
		md.sys_size = N + 1;
		md.hamiltonian = new double[md.sys_size * md.sys_size];
		md.hamiltonian_drv = new double[md.sys_size * md.sys_size];

		double prm_E = double(cp.params.find("prm_E")->second);
		double prm_U = double(cp.params.find("prm_U")->second);
		double prm_J = double(cp.params.find("prm_J")->second);

		prm_U = prm_U / double(N);

		for (int st_id_1 = 0; st_id_1 < md.sys_size; st_id_1++)
		{
			for (int st_id_2 = 0; st_id_2 < md.sys_size; st_id_2++)
			{
				int index = st_id_1 * md.sys_size + st_id_2;
				md.hamiltonian[index] = 0.0;
				md.hamiltonian_drv[index] = 0.0;
			}
		}

		for (int st_id = 0; st_id < md.sys_size; st_id++)
		{
			int index = st_id * md.sys_size + st_id;
			md.hamiltonian[index] = 2.0 * prm_U * double(st_id * (st_id - 1) + (N - st_id)*(N - st_id - 1));
			md.hamiltonian_drv[index] = double((N - st_id) - st_id);
		}

		for (int st_id = 0; st_id < md.sys_size; st_id++)
		{
			int down_left = (st_id + 1) * md.sys_size + st_id;
			int up_right = st_id * md.sys_size + (st_id + 1);
			md.hamiltonian[down_left]

		}

		for i = 1:N - 1
			% if (i + 1 <= Nn)
			% H1(i, i + 1) = H1(i, i + 1) - J*sqrt((Nn - i)*(i + 1)); %b1'*b2
			%     end
			%     if (i - 1>0)
			% H1(i, i - 1) = H1(i, i - 1) - J*sqrt((i)*(Nn - i + 1)); %b2'*b1
			%     end

			H1(i + 1, i) = H1(i + 1, i) - J*sqrt((N - i)*(i));
		H1(i, i + 1) = H1(i, i + 1) - J*sqrt((i)*(N - i));
		end

	}
}