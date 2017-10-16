#include "data.h"

void init_main_data(ConfigParam &cp, Data &dt)
{
	dt.size = 2;

	dt.step = (2.0 * PI) / (cp.omega * cp.num_steps);

	dt.time = 0;

	dt.data = new double[dt.size];
	dt.args = new double[dt.size];
	dt.k1s = new  double[dt.size];
	dt.k2s = new  double[dt.size];
	dt.k3s = new  double[dt.size];
	dt.k4s = new  double[dt.size];

	for (int v_id = 0; v_id < dt.size; v_id++)
	{
		dt.data[v_id] = 0.0;
		dt.args[v_id] = 0.0;
		dt.k1s[v_id] = 0.0;
		dt.k2s[v_id] = 0.0;
		dt.k3s[v_id] = 0.0;
		dt.k4s[v_id] = 0.0;
	}
}

void delete_main_data(Data &dt)
{
	delete_data(dt.data);
	delete_data(dt.args);
	delete_data(dt.k1s);
	delete_data(dt.k2s);
	delete_data(dt.k3s);
	delete_data(dt.k4s);
}

void init_cond(ConfigParam &cp, Data &dt)
{
	double left_border = -PI;
	double right_border = PI;
	int max_num_seeds = 1000000;

	VSLStreamStatePtr stream;
	vslNewStream(&stream, VSL_BRNG_MCG31, 77778888);
	vslLeapfrogStream(stream, cp.seed, max_num_seeds);
	vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, dt.size, dt.data, left_border, right_border);
}
