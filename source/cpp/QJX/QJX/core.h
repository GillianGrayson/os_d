#pragma once
#include "config.h"
#include "data.h"

class CoreBehavior
{
public:
	virtual void init_splits(AllData * ad) const = 0;
	virtual void free_splits(AllData * ad) const = 0;

	virtual void init_splits_deep(AllData * ad) const = 0;
	virtual void free_splits_deep(AllData * ad) const = 0;

	virtual void ex_period(AllData * ad, int tr_id, int th_id, int period_id) const = 0;
	virtual void ex_period_trp_deep(AllData * ad, int tr_id, int th_id, int period_id) const = 0;
	virtual void ex_period_obs_deep(AllData * ad, int tr_id, int th_id, int period_id) const = 0;
	virtual void ex_period_obs_deep_lpn(AllData * ad, int period_id) const = 0;
	virtual void ex_period_obs_deep_mult_lpn(AllData * ad, int period_id) const = 0;
	virtual void ex_period_obs_deep_lpn_per_period(AllData * ad, int period_id, int num_periods) const = 0;
	virtual void ex_period_obs_deep_cd(AllData * ad, int tr_id, int th_id, int period_id) const = 0;
	virtual void ex_period_obs_deep_sigma(AllData * ad, int tr_id, int th_id, int period_id) const = 0;

	virtual void rk_period(AllData * ad, int tr_id, int th_id, int period_id) const = 0;
	virtual void rk_period_trp_deep(AllData * ad, int tr_id, int th_id, int period_id) const = 0;
	virtual void rk_period_obs_deep(AllData * ad, int tr_id, int th_id, int period_id) const = 0;
	virtual void rk_period_obs_deep_lpn(AllData * ad, int period_id) const = 0;
	virtual void rk_period_obs_deep_cd(AllData * ad, int tr_id, int th_id, int period_id) const = 0;
	virtual void rk_period_obs_deep_sigma(AllData * ad, int tr_id, int th_id, int period_id) const = 0;

	virtual void calc_chars_std_start(AllData * ad, int tr_id) const = 0;
	virtual void calc_chars_std(AllData * ad, int tr_id) const = 0;
	virtual void calc_chars_lpn_start(AllData * ad, int tr_id, int base_tr_id) const = 0;
	virtual void calc_chars_lpn(AllData * ad, int tr_id, int base_tr_id) const = 0;

	virtual double calc_T(AllData * ad) const = 0;

	virtual void evo_chars_std(AllData * ad, int tr_id, int dump_id) const = 0;
	virtual void evo_chars_lpn(AllData * ad, int tr_id, int dump_id) const = 0;

	virtual double calc_delta_s(AllData * ad, int tr_id, int base_tr_id) const = 0;
	virtual double calc_delta_f(AllData * ad, int tr_id, int base_tr_id) const = 0;

	virtual void calc_ci(AllData * ad, int tr_id) const = 0;

	virtual void dump_std(AllData * ad) const = 0;
	virtual void dump_lpn(AllData * ad) const = 0;
	virtual void dump_std_evo(AllData * ad) const = 0;
	virtual void dump_lpn_evo(AllData * ad) const = 0;
};

class DimerCoreBehaviour : public CoreBehavior
{
public:
	virtual void init_splits(AllData * ad) const;
	virtual void free_splits(AllData * ad) const;

	virtual void init_splits_deep(AllData * ad) const;
	virtual void free_splits_deep(AllData * ad) const;

	virtual void ex_period(AllData * ad, int tr_id, int th_id, int period_id) const;
	virtual void ex_period_trp_deep(AllData * ad, int tr_id, int th_id, int period_id) const;
	virtual void ex_period_obs_deep(AllData * ad, int tr_id, int th_id, int period_id) const;
	virtual void ex_period_obs_deep_lpn(AllData * ad, int period_id) const;
	virtual void ex_period_obs_deep_mult_lpn(AllData * ad, int period_id) const;
	virtual void ex_period_obs_deep_lpn_per_period(AllData* ad, int period_id, int num_periods) const;
	virtual void ex_period_obs_deep_cd(AllData * ad, int tr_id, int th_id, int period_id) const;
	virtual void ex_period_obs_deep_sigma(AllData * ad, int tr_id, int th_id, int period_id) const;

	virtual void rk_period(AllData * ad, int tr_id, int th_id, int period_id) const;
	virtual void rk_period_trp_deep(AllData * ad, int tr_id, int th_id, int period_id) const;
	virtual void rk_period_obs_deep(AllData * ad, int tr_id, int th_id, int period_id) const;
	virtual void rk_period_obs_deep_lpn(AllData * ad, int period_id) const;
	virtual void rk_period_obs_deep_cd(AllData * ad, int tr_id, int th_id, int period_id) const;
	virtual void rk_period_obs_deep_sigma(AllData * ad, int tr_id, int th_id, int period_id) const;

	virtual void calc_chars_std_start(AllData * ad, int tr_id) const;
	virtual void calc_chars_std(AllData * ad, int tr_id) const;
	virtual void calc_chars_lpn_start(AllData * ad, int tr_id, int base_tr_id) const;
	virtual void calc_chars_lpn(AllData * ad, int tr_id, int base_tr_id) const;

	virtual double calc_T(AllData * ad) const;

	virtual void evo_chars_std(AllData * ad, int tr_id, int dump_id) const;
	virtual void evo_chars_lpn(AllData * ad, int tr_id, int dump_id) const;

	virtual double calc_delta_s(AllData * ad, int tr_id, int base_tr_id) const;
	virtual double calc_delta_f(AllData * ad, int tr_id, int base_tr_id) const;

	virtual void calc_ci(AllData * ad, int tr_id) const;

	virtual void dump_std(AllData * ad) const;
	virtual void dump_lpn(AllData * ad) const;
	virtual void dump_std_evo(AllData * ad) const;
	virtual void dump_lpn_evo(AllData * ad) const;
};


class DimerSyncCoreBehaviour : public CoreBehavior
{
public:
	virtual void init_splits(AllData* ad) const;
	virtual void free_splits(AllData* ad) const;

	virtual void init_splits_deep(AllData* ad) const;
	virtual void free_splits_deep(AllData* ad) const;

	virtual void ex_period(AllData* ad, int tr_id, int th_id, int period_id) const;
	virtual void ex_period_trp_deep(AllData* ad, int tr_id, int th_id, int period_id) const;
	virtual void ex_period_obs_deep(AllData* ad, int tr_id, int th_id, int period_id) const;
	virtual void ex_period_obs_deep_lpn(AllData* ad, int period_id) const;
	virtual void ex_period_obs_deep_mult_lpn(AllData* ad, int period_id) const;
	virtual void ex_period_obs_deep_lpn_per_period(AllData* ad, int period_id, int num_periods) const;
	virtual void ex_period_obs_deep_cd(AllData* ad, int tr_id, int th_id, int period_id) const;
	virtual void ex_period_obs_deep_sigma(AllData* ad, int tr_id, int th_id, int period_id) const;

	virtual void rk_period(AllData* ad, int tr_id, int th_id, int period_id) const;
	virtual void rk_period_trp_deep(AllData* ad, int tr_id, int th_id, int period_id) const;
	virtual void rk_period_obs_deep(AllData* ad, int tr_id, int th_id, int period_id) const;
	virtual void rk_period_obs_deep_lpn(AllData* ad, int period_id) const;
	virtual void rk_period_obs_deep_cd(AllData* ad, int tr_id, int th_id, int period_id) const;
	virtual void rk_period_obs_deep_sigma(AllData* ad, int tr_id, int th_id, int period_id) const;

	virtual void calc_chars_std_start(AllData* ad, int tr_id) const;
	virtual void calc_chars_std(AllData* ad, int tr_id) const;
	virtual void calc_chars_lpn_start(AllData* ad, int tr_id, int base_tr_id) const;
	virtual void calc_chars_lpn(AllData* ad, int tr_id, int base_tr_id) const;

	virtual double calc_T(AllData* ad) const;

	virtual void evo_chars_std(AllData* ad, int tr_id, int dump_id) const;
	virtual void evo_chars_lpn(AllData* ad, int tr_id, int dump_id) const;

	virtual double calc_delta_s(AllData* ad, int tr_id, int base_tr_id) const;
	virtual double calc_delta_f(AllData* ad, int tr_id, int base_tr_id) const;

	virtual void calc_ci(AllData* ad, int tr_id) const;

	virtual void dump_std(AllData* ad) const;
	virtual void dump_lpn(AllData* ad) const;
	virtual void dump_std_evo(AllData* ad) const;
	virtual void dump_lpn_evo(AllData* ad) const;
};


class JCSCoreBehaviour : public CoreBehavior
{
public:
	virtual void init_splits(AllData * ad) const;
	virtual void free_splits(AllData * ad) const;

	virtual void init_splits_deep(AllData * ad) const;
	virtual void free_splits_deep(AllData * ad) const;

	virtual void ex_period(AllData * ad, int tr_id, int th_id, int period_id) const;
	virtual void ex_period_trp_deep(AllData * ad, int tr_id, int th_id, int period_id) const;
	virtual void ex_period_obs_deep(AllData * ad, int tr_id, int th_id, int period_id) const;
	virtual void ex_period_obs_deep_lpn(AllData * ad, int period_id) const;
	virtual void ex_period_obs_deep_mult_lpn(AllData * ad, int period_id) const;
	virtual void ex_period_obs_deep_lpn_per_period(AllData* ad, int period_id, int num_periods) const;
	virtual void ex_period_obs_deep_cd(AllData * ad, int tr_id, int th_id, int period_id) const;
	virtual void ex_period_obs_deep_sigma(AllData * ad, int tr_id, int th_id, int period_id) const;

	virtual void rk_period(AllData * ad, int tr_id, int th_id, int period_id) const;
	virtual void rk_period_trp_deep(AllData * ad, int tr_id, int th_id, int period_id) const;
	virtual void rk_period_obs_deep(AllData * ad, int tr_id, int th_id, int period_id) const;
	virtual void rk_period_obs_deep_lpn(AllData * ad, int period_id) const;
	virtual void rk_period_obs_deep_cd(AllData * ad, int tr_id, int th_id, int period_id) const;
	virtual void rk_period_obs_deep_sigma(AllData * ad, int tr_id, int th_id, int period_id) const;

	virtual void calc_chars_std_start(AllData * ad, int tr_id) const;
	virtual void calc_chars_std(AllData * ad, int tr_id) const;
	virtual void calc_chars_lpn_start(AllData * ad, int tr_id, int base_tr_id) const;
	virtual void calc_chars_lpn(AllData * ad, int tr_id, int base_tr_id) const;

	virtual double calc_T(AllData * ad) const;

	virtual void evo_chars_std(AllData * ad, int tr_id, int dump_id) const;
	virtual void evo_chars_lpn(AllData * ad, int tr_id, int dump_id) const;

	virtual double calc_delta_s(AllData * ad, int tr_id, int base_tr_id) const;
	virtual double calc_delta_f(AllData * ad, int tr_id, int base_tr_id) const;

	virtual void calc_ci(AllData * ad, int tr_id) const;

	virtual void dump_std(AllData * ad) const;
	virtual void dump_lpn(AllData * ad) const;
	virtual void dump_std_evo(AllData * ad) const;
	virtual void dump_lpn_evo(AllData * ad) const;
};

class PSCoreBehaviour : public CoreBehavior
{
public:
	virtual void init_splits(AllData * ad) const;
	virtual void free_splits(AllData * ad) const;

	virtual void init_splits_deep(AllData * ad) const;
	virtual void free_splits_deep(AllData * ad) const;

	virtual void ex_period(AllData * ad, int tr_id, int th_id, int period_id) const;
	virtual void ex_period_trp_deep(AllData * ad, int tr_id, int th_id, int period_id) const;
	virtual void ex_period_obs_deep(AllData * ad, int tr_id, int th_id, int period_id) const;
	virtual void ex_period_obs_deep_lpn(AllData * ad, int period_id) const;
	virtual void ex_period_obs_deep_mult_lpn(AllData * ad, int period_id) const;
	virtual void ex_period_obs_deep_lpn_per_period(AllData* ad, int period_id, int num_periods) const;
	virtual void ex_period_obs_deep_cd(AllData * ad, int tr_id, int th_id, int period_id) const;
	virtual void ex_period_obs_deep_sigma(AllData * ad, int tr_id, int th_id, int period_id) const;

	virtual void rk_period(AllData * ad, int tr_id, int th_id, int period_id) const;
	virtual void rk_period_trp_deep(AllData * ad, int tr_id, int th_id, int period_id) const;
	virtual void rk_period_obs_deep(AllData * ad, int tr_id, int th_id, int period_id) const;
	virtual void rk_period_obs_deep_lpn(AllData * ad, int period_id) const;
	virtual void rk_period_obs_deep_cd(AllData * ad, int tr_id, int th_id, int period_id) const;
	virtual void rk_period_obs_deep_sigma(AllData * ad, int tr_id, int th_id, int period_id) const;

	virtual void calc_chars_std_start(AllData * ad, int tr_id) const;
	virtual void calc_chars_std(AllData * ad, int tr_id) const;
	virtual void calc_chars_lpn_start(AllData * ad, int tr_id, int base_tr_id) const;
	virtual void calc_chars_lpn(AllData * ad, int tr_id, int base_tr_id) const;

	virtual double calc_T(AllData * ad) const;

	virtual void evo_chars_std(AllData * ad, int tr_id, int dump_id) const;
	virtual void evo_chars_lpn(AllData * ad, int tr_id, int dump_id) const;

	virtual double calc_delta_s(AllData * ad, int tr_id, int base_tr_id) const;
	virtual double calc_delta_f(AllData * ad, int tr_id, int base_tr_id) const;

	virtual void calc_ci(AllData * ad, int tr_id) const;

	virtual void dump_std(AllData * ad) const;
	virtual void dump_lpn(AllData * ad) const;
	virtual void dump_std_evo(AllData * ad) const;
	virtual void dump_lpn_evo(AllData * ad) const;
};

class MBLCoreBehaviour : public CoreBehavior
{
public:
	virtual void init_splits(AllData * ad) const;
	virtual void free_splits(AllData * ad) const;

	virtual void init_splits_deep(AllData * ad) const;
	virtual void free_splits_deep(AllData * ad) const;

	virtual void ex_period(AllData * ad, int tr_id, int th_id, int period_id) const;
	virtual void ex_period_trp_deep(AllData * ad, int tr_id, int th_id, int period_id) const;
	virtual void ex_period_obs_deep(AllData * ad, int tr_id, int th_id, int period_id) const;
	virtual void ex_period_obs_deep_lpn(AllData * ad, int period_id) const;
	virtual void ex_period_obs_deep_mult_lpn(AllData * ad, int period_id) const;
	virtual void ex_period_obs_deep_lpn_per_period(AllData* ad, int period_id, int num_periods) const;
	virtual void ex_period_obs_deep_cd(AllData * ad, int tr_id, int th_id, int period_id) const;
	virtual void ex_period_obs_deep_sigma(AllData * ad, int tr_id, int th_id, int period_id) const;

	virtual void rk_period(AllData * ad, int tr_id, int th_id, int period_id) const;
	virtual void rk_period_trp_deep(AllData * ad, int tr_id, int th_id, int period_id) const;
	virtual void rk_period_obs_deep(AllData * ad, int tr_id, int th_id, int period_id) const;
	virtual void rk_period_obs_deep_lpn(AllData * ad, int period_id) const;
	virtual void rk_period_obs_deep_cd(AllData * ad, int tr_id, int th_id, int period_id) const;
	virtual void rk_period_obs_deep_sigma(AllData * ad, int tr_id, int th_id, int period_id) const;

	virtual void calc_chars_std_start(AllData * ad, int tr_id) const;
	virtual void calc_chars_std(AllData * ad, int tr_id) const;
	virtual void calc_chars_lpn_start(AllData * ad, int tr_id, int base_tr_id) const;
	virtual void calc_chars_lpn(AllData * ad, int tr_id, int base_tr_id) const;

	virtual double calc_T(AllData * ad) const;

	virtual void evo_chars_std(AllData * ad, int tr_id, int dump_id) const;
	virtual void evo_chars_lpn(AllData * ad, int tr_id, int dump_id) const;

	virtual double calc_delta_s(AllData * ad, int tr_id, int base_tr_id) const;
	virtual double calc_delta_f(AllData * ad, int tr_id, int base_tr_id) const;

	virtual void calc_ci(AllData * ad, int tr_id) const;

	virtual void dump_std(AllData * ad) const;
	virtual void dump_lpn(AllData * ad) const;
	virtual void dump_std_evo(AllData * ad) const;
	virtual void dump_lpn_evo(AllData * ad) const;
};

class LndHamCoreBehaviour : public CoreBehavior
{
public:
	virtual void init_splits(AllData* ad) const;
	virtual void free_splits(AllData* ad) const;

	virtual void init_splits_deep(AllData* ad) const;
	virtual void free_splits_deep(AllData* ad) const;

	virtual void ex_period(AllData* ad, int tr_id, int th_id, int period_id) const;
	virtual void ex_period_trp_deep(AllData* ad, int tr_id, int th_id, int period_id) const;
	virtual void ex_period_obs_deep(AllData* ad, int tr_id, int th_id, int period_id) const;
	virtual void ex_period_obs_deep_lpn(AllData* ad, int period_id) const;
	virtual void ex_period_obs_deep_mult_lpn(AllData* ad, int period_id) const;
	virtual void ex_period_obs_deep_lpn_per_period(AllData* ad, int period_id, int num_periods) const;
	virtual void ex_period_obs_deep_cd(AllData* ad, int tr_id, int th_id, int period_id) const;
	virtual void ex_period_obs_deep_sigma(AllData* ad, int tr_id, int th_id, int period_id) const;

	virtual void rk_period(AllData* ad, int tr_id, int th_id, int period_id) const;
	virtual void rk_period_trp_deep(AllData* ad, int tr_id, int th_id, int period_id) const;
	virtual void rk_period_obs_deep(AllData* ad, int tr_id, int th_id, int period_id) const;
	virtual void rk_period_obs_deep_lpn(AllData* ad, int period_id) const;
	virtual void rk_period_obs_deep_cd(AllData* ad, int tr_id, int th_id, int period_id) const;
	virtual void rk_period_obs_deep_sigma(AllData* ad, int tr_id, int th_id, int period_id) const;

	virtual void calc_chars_std_start(AllData* ad, int tr_id) const;
	virtual void calc_chars_std(AllData* ad, int tr_id) const;
	virtual void calc_chars_lpn_start(AllData* ad, int tr_id, int base_tr_id) const;
	virtual void calc_chars_lpn(AllData* ad, int tr_id, int base_tr_id) const;

	virtual double calc_T(AllData* ad) const;

	virtual void evo_chars_std(AllData* ad, int tr_id, int dump_id) const;
	virtual void evo_chars_lpn(AllData* ad, int tr_id, int dump_id) const;

	virtual double calc_delta_s(AllData* ad, int tr_id, int base_tr_id) const;
	virtual double calc_delta_f(AllData* ad, int tr_id, int base_tr_id) const;

	virtual void calc_ci(AllData* ad, int tr_id) const;

	virtual void dump_std(AllData* ad) const;
	virtual void dump_lpn(AllData* ad) const;
	virtual void dump_std_evo(AllData* ad) const;
	virtual void dump_lpn_evo(AllData* ad) const;
};

class IntegrableCoreBehaviour : public CoreBehavior
{
public:
	virtual void init_splits(AllData* ad) const;
	virtual void free_splits(AllData* ad) const;

	virtual void init_splits_deep(AllData* ad) const;
	virtual void free_splits_deep(AllData* ad) const;

	virtual void ex_period(AllData* ad, int tr_id, int th_id, int period_id) const;
	virtual void ex_period_trp_deep(AllData* ad, int tr_id, int th_id, int period_id) const;
	virtual void ex_period_obs_deep(AllData* ad, int tr_id, int th_id, int period_id) const;
	virtual void ex_period_obs_deep_lpn(AllData* ad, int period_id) const;
	virtual void ex_period_obs_deep_mult_lpn(AllData* ad, int period_id) const;
	virtual void ex_period_obs_deep_lpn_per_period(AllData* ad, int period_id, int num_periods) const;
	virtual void ex_period_obs_deep_cd(AllData* ad, int tr_id, int th_id, int period_id) const;
	virtual void ex_period_obs_deep_sigma(AllData* ad, int tr_id, int th_id, int period_id) const;

	virtual void rk_period(AllData* ad, int tr_id, int th_id, int period_id) const;
	virtual void rk_period_trp_deep(AllData* ad, int tr_id, int th_id, int period_id) const;
	virtual void rk_period_obs_deep(AllData* ad, int tr_id, int th_id, int period_id) const;
	virtual void rk_period_obs_deep_lpn(AllData* ad, int period_id) const;
	virtual void rk_period_obs_deep_cd(AllData* ad, int tr_id, int th_id, int period_id) const;
	virtual void rk_period_obs_deep_sigma(AllData* ad, int tr_id, int th_id, int period_id) const;

	virtual void calc_chars_std_start(AllData* ad, int tr_id) const;
	virtual void calc_chars_std(AllData* ad, int tr_id) const;
	virtual void calc_chars_lpn_start(AllData* ad, int tr_id, int base_tr_id) const;
	virtual void calc_chars_lpn(AllData* ad, int tr_id, int base_tr_id) const;

	virtual double calc_T(AllData* ad) const;

	virtual void evo_chars_std(AllData* ad, int tr_id, int dump_id) const;
	virtual void evo_chars_lpn(AllData* ad, int tr_id, int dump_id) const;

	virtual double calc_delta_s(AllData* ad, int tr_id, int base_tr_id) const;
	virtual double calc_delta_f(AllData* ad, int tr_id, int base_tr_id) const;

	virtual void calc_ci(AllData* ad, int tr_id) const;

	virtual void dump_std(AllData* ad) const;
	virtual void dump_lpn(AllData* ad) const;
	virtual void dump_std_evo(AllData* ad) const;
	virtual void dump_lpn_evo(AllData* ad) const;
};

Split * init_split_structure_dimer(AllData * ad);
Split * init_split_structure_dimer_deep(AllData * ad);
Split* init_split_structure_dimersync(AllData* ad);
Split* init_split_structure_dimersync_deep(AllData* ad);
Split * init_split_structure_jcs(AllData * ad);
Split * init_split_structure_jcs_deep(AllData * ad);
Split * init_split_structure_ps(AllData * ad);
Split * init_split_structure_ps_deep(AllData * ad);
Split * init_split_structure_mbl(AllData * ad);
Split * init_split_structure_mbl_deep(AllData * ad);
Split* init_split_structure_lndham(AllData* ad);
Split* init_split_structure_lndham_deep(AllData* ad);
Split* init_split_structure_integrable(AllData* ad);
Split* init_split_structure_integrable_deep(AllData* ad);

void rk_right_part_dimer(AllData * ad, int sub_step, int tr_id, int th_id);
void rk_right_part_dimersync(AllData* ad, int sub_step, int tr_id, int th_id);
void rk_right_part_jcs(AllData * ad, int sub_step, int tr_id, int th_id);
void rk_right_part_ps(AllData * ad, int sub_step, int tr_id, int th_id);
void rk_right_part_mbl(AllData * ad, int sub_step, int tr_id, int th_id);
void rk_right_part_lndham(AllData* ad, int sub_step, int tr_id, int th_id);
void rk_right_part_integrable(AllData* ad, int sub_step, int tr_id, int th_id);

void rk_int_dimer(AllData * ad, int tr_id, int th_id, double step);
void rk_int_dimersync(AllData* ad, int tr_id, int th_id, double step);
void rk_int_jcs(AllData * ad, int tr_id, int th_id, double step);
void rk_int_ps(AllData * ad, int tr_id, int th_id, double step);
void rk_int_mbl(AllData * ad, int tr_id, int th_id, double step);
void rk_int_lndham(AllData* ad, int tr_id, int th_id, double step);
void rk_int_integrable(AllData* ad, int tr_id, int th_id, double step);

void rk_step_dimer(AllData * ad, int tr_id, int th_id, double step);
void rk_step_dimersync(AllData* ad, int tr_id, int th_id, double step);
void rk_step_jcs(AllData * ad, int tr_id, int th_id, double step);
void rk_step_ps(AllData * ad, int tr_id, int th_id, double step);
void rk_step_mbl(AllData * ad, int tr_id, int th_id, double step);
void rk_step_lndham(AllData* ad, int tr_id, int th_id, double step);
void rk_step_integrable(AllData* ad, int tr_id, int th_id, double step);

void calc_ci_double(AllData * ad, int tr_id);

void dump_phi(AllData * ad);

void dump_phi_evo(AllData * ad, bool append);


void calc_random_obs(AllData * ad, int tr_id);

void calc_random_obs_lpn_start(AllData * ad, int tr_id);

void calc_random_obs_lpn(AllData * ad, int tr_id);