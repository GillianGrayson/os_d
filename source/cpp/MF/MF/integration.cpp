#include "integration.h"

void right_part(ConfigParam &cp, double * ks, double * x, double time)
{
	double driving = sin(cp.omega * time + cp.phase);

	if (cp.mt == 1)
	{
		driving *= cp.A;
	}
	else
	{
		driving = cp.A * (2.0 * (driving > 0.0) - 1.0);
	}

	ks[0] = 2.0 * cp.J * sin(x[1]) + 4.0 * cp.gamma * cos(x[1]) * cos(x[0]);
	ks[1] = 2.0 * cp.J * cos(x[0]) * cos(x[1]) / sin(x[0]) - 2.0 * cp.E - 2.0 * driving + cp.U * cos(x[0]) - 4.0 * cp.gamma * sin(x[1]) / sin(x[0]);
}

void upd_arg(int size, double * x_arg, double * x, double * ks, double coeff)
{
	for (int st_id = 0; st_id < size; st_id++)
	{
		x_arg[st_id] = x[st_id] + coeff * ks[st_id];
	}
}

void rk_final(int size, double * x, double * k1s, double * k2s, double * k3s, double * k4s, double step)
{
	for (int st_id = 0; st_id < size; st_id++)
	{
		x[st_id] += (k1s[st_id] + 2.0 * k2s[st_id] + 2.0 * k3s[st_id] + k4s[st_id]) * step / 6.0;
	}
}

void rk_step(ConfigParam &cp, Data &dt)
{

	right_part(cp, dt.k1s, dt.data, dt.time);
	upd_arg(dt.size, dt.args, dt.data, dt.k1s, dt.step * 0.5);

	dt.time += dt.step * 0.5;

	right_part(cp, dt.k2s, dt.args, dt.time);
	upd_arg(dt.size, dt.args, dt.data, dt.k2s, dt.step * 0.5);

	right_part(cp, dt.k3s, dt.args, dt.time);
	upd_arg(dt.size, dt.args, dt.data, dt.k3s, dt.step * 1.0);

	dt.time += dt.step * 0.5;

	right_part(cp, dt.k4s, dt.args, dt.time);

	rk_final(dt.size, dt.data, dt.k1s, dt.k2s, dt.k3s, dt.k4s, dt.step);
}

void int_period(ConfigParam &cp, Data &dt, int per_id)
{
	for (int step_id = 0; step_id < cp.num_steps; step_id++)
	{
		dt.time = per_id * cp.T + step_id * dt.step;
		rk_step(cp, dt);
	}
}

void int_trans_proc(ConfigParam &cp, Data &dt)
{
	for (int per_id = 0; per_id < cp.npt; per_id++)
	{
		int_period(cp, dt, per_id);
	}
}