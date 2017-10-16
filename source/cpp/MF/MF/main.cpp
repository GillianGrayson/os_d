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

	if (rp.task == 0)
	{
		only_data_exp(rp, cp);
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