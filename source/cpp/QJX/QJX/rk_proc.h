#pragma once
#include "data.h"

void init_rk_data(AllData * ad);

void init_rk(AllData * ad);

void init_rk_cd(AllData * ad);

void free_rk(AllData * ad);

void free_rk_cd(AllData * ad);

void set_init_args(AllData * ad, int tr_id, int th_id);

void arg_upd(AllData * ad, int sub_step, int tr_id, int th_id);

void right_part(AllData * ad, int sub_step, int tr_id, int th_id);

void rk_final(AllData * ad, int tr_id, int th_id, double step);

void rk_recovery(AllData * ad, int tr_id, int th_id);

void save_phi_prev(AllData * ad, int tr_id, int th_id);

void restore_from_prev(AllData * ad, int tr_id, int th_id, double step);

void rk_int(AllData * ad, int tr_id, int th_id, double step);

void rk_step(AllData * ad, int tr_id, int th_id, double step);

void rk_period(AllData * ad, int tr_id, int th_id, int period_id);

void rk_period_cd(AllData * ad, int tr_id, int th_id, double start_time);