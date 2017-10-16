#pragma once
#include "config.h"
#include "utils.h"
#include "data.h"

void right_part(ConfigParam &cp, double * ks, double * x, double time);

void upd_arg(int size, double * x_arg, double * x, double * ks, double coeff);

void rk_final(int size, double * x, double * k1s, double * k2s, double * k3s, double * k4s, double step);

void rk_step(ConfigParam &cp, Data &dt);

void int_period(ConfigParam &cp, Data &dt, int per_id);

void int_trans_proc(ConfigParam &cp, Data &dt);