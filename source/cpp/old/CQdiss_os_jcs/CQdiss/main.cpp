#include <omp.h>

#include "read_config.h"
#include "Model.h"
#include "InitFs.h"
#include "initH.h"
#include "Init_a1_a2.h"
#include "Init_f_d.h"
#include "TranspFs.h"

#include "CalcQs.h"
#include "CalcKs.h"
#include "CalcRs.h"
#include "CalcGs.h"

#include "calcODE.h"

#include "CalcRho.h"

int main(int argc, char ** argv)
{
	cout << "Current path: " << argv[0] << endl << endl;
	if (argc > 1)
	{
		omp_set_num_threads(atoi(argv[1]));
	}

	RunParam rp;
	ConfigParam cp;
	MainData md;
	init_params(rp, cp, "config.txt");

	double time = omp_get_wtime();

	if (rp.task == TASK_ID_STD)
	{
		f_basis_prop_std(rp, cp, md);
	}
	else if (rp.task == TASK_ID_DEEP)
	{
		f_basis_prop_deep(rp, cp, md);
	}
	else
	{
		stringstream msg;
		msg << "Wrong task value: " << rp.task << endl;
		Error(msg.str());
	}

	time = omp_get_wtime() - time;
	cout << "Total time: " << time << endl << endl;

	return 0;
}
