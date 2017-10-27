#include <omp.h>
#include "config.h"
#include "experiments.h"

int main(int argc, char ** argv)
{
	cout << "current path: " << argv[0] << endl << endl;
	if (argc > 1)
	{
		omp_set_num_threads(atoi(argv[1]));
	}

	ConfigParam cp;
	RunParam rp;
	init_params(rp, cp, "config.txt");

	double time = omp_get_wtime();

	if (rp.task == BASIC_EXP_ID)
	{
		basic_exp(rp, cp);
	}
	else if (rp.task == LPN_FIN_EXP_ID)
	{
		lpn_fin_exp(rp, cp);
	}
	else if (rp.task == CD_EXP_ID)
	{
		cd_exp(rp, cp);
	}
	else if (rp.task == BASIC_AND_LPN_FIN_EXP_ID)
	{
		basic_and_lpn_fin_exp(rp, cp);
	}
	else if (rp.task == CD_D_EXP_ID)
	{
		cd_d_exp(rp, cp);
	}
	else
	{
		stringstream msg;
		msg << "wrong task value: " << rp.task << endl;
		Error(msg.str());
	}

	time = omp_get_wtime() - time;
	cout << "total time: " << time << endl << endl;

	return 0;
}