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
	virtual void ex_period_obs_deep_lpn(AllData * ad, int tr_id, int th_id, int period_id) const = 0;
	virtual void ex_period_obs_deep_cd(AllData * ad, int tr_id, int th_id, int period_id) const = 0;
	virtual void ex_period_obs_deep_sigma(AllData * ad, int tr_id, int th_id, int period_id) const = 0;

	virtual void rk_period(AllData * ad, int tr_id, int th_id, int period_id) const = 0;
	virtual void rk_period_trp_deep(AllData * ad, int tr_id, int th_id, int period_id) const = 0;
	virtual void rk_period_obs_deep(AllData * ad, int tr_id, int th_id, int period_id) const = 0;
	virtual void rk_period_obs_deep_lpn(AllData * ad, int tr_id, int th_id, int period_id) const = 0;
	virtual void rk_period_obs_deep_cd(AllData * ad, int tr_id, int th_id, int period_id) const = 0;
	virtual void rk_period_obs_deep_sigma(AllData * ad, int tr_id, int th_id, int period_id) const = 0;

	virtual void calc_chars_std_start(AllData * ad, int tr_id) const = 0;
	virtual void calc_chars_std(AllData * ad, int tr_id) const = 0;
	virtual void calc_chars_lpn_start(AllData * ad, int tr_id) const = 0;
	virtual void calc_chars_lpn(AllData * ad, int tr_id) const = 0;

	virtual void evo_chars_std(AllData * ad, int tr_id, int dump_id) const = 0;
	virtual void evo_chars_lpn(AllData * ad, int tr_id, int dump_id) const = 0;

	virtual double calc_delta(AllData * ad, int tr_id) const = 0;

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
	virtual void ex_period_obs_deep_lpn(AllData * ad, int tr_id, int th_id, int period_id) const;
	virtual void ex_period_obs_deep_cd(AllData * ad, int tr_id, int th_id, int period_id) const;
	virtual void ex_period_obs_deep_sigma(AllData * ad, int tr_id, int th_id, int period_id) const;

	virtual void rk_period(AllData * ad, int tr_id, int th_id, int period_id) const;
	virtual void rk_period_trp_deep(AllData * ad, int tr_id, int th_id, int period_id) const;
	virtual void rk_period_obs_deep(AllData * ad, int tr_id, int th_id, int period_id) const;
	virtual void rk_period_obs_deep_lpn(AllData * ad, int tr_id, int th_id, int period_id) const;
	virtual void rk_period_obs_deep_cd(AllData * ad, int tr_id, int th_id, int period_id) const;
	virtual void rk_period_obs_deep_sigma(AllData * ad, int tr_id, int th_id, int period_id) const;

	virtual void calc_chars_std_start(AllData * ad, int tr_id) const;
	virtual void calc_chars_std(AllData * ad, int tr_id) const;
	virtual void calc_chars_lpn_start(AllData * ad, int tr_id) const;
	virtual void calc_chars_lpn(AllData * ad, int tr_id) const;

	virtual void evo_chars_std(AllData * ad, int tr_id, int dump_id) const;
	virtual void evo_chars_lpn(AllData * ad, int tr_id, int dump_id) const;

	virtual double calc_delta(AllData * ad, int tr_id) const;

	virtual void calc_ci(AllData * ad, int tr_id) const;

	virtual void dump_std(AllData * ad) const;
	virtual void dump_lpn(AllData * ad) const;
	virtual void dump_std_evo(AllData * ad) const;
	virtual void dump_lpn_evo(AllData * ad) const;
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
	virtual void ex_period_obs_deep_lpn(AllData * ad, int tr_id, int th_id, int period_id) const;
	virtual void ex_period_obs_deep_cd(AllData * ad, int tr_id, int th_id, int period_id) const;
	virtual void ex_period_obs_deep_sigma(AllData * ad, int tr_id, int th_id, int period_id) const;

	virtual void rk_period(AllData * ad, int tr_id, int th_id, int period_id) const;
	virtual void rk_period_trp_deep(AllData * ad, int tr_id, int th_id, int period_id) const;
	virtual void rk_period_obs_deep(AllData * ad, int tr_id, int th_id, int period_id) const;
	virtual void rk_period_obs_deep_lpn(AllData * ad, int tr_id, int th_id, int period_id) const;
	virtual void rk_period_obs_deep_cd(AllData * ad, int tr_id, int th_id, int period_id) const;
	virtual void rk_period_obs_deep_sigma(AllData * ad, int tr_id, int th_id, int period_id) const;

	virtual void calc_chars_std_start(AllData * ad, int tr_id) const;
	virtual void calc_chars_std(AllData * ad, int tr_id) const;
	virtual void calc_chars_lpn_start(AllData * ad, int tr_id) const;
	virtual void calc_chars_lpn(AllData * ad, int tr_id) const;

	virtual void evo_chars_std(AllData * ad, int tr_id, int dump_id) const;
	virtual void evo_chars_lpn(AllData * ad, int tr_id, int dump_id) const;

	virtual double calc_delta(AllData * ad, int tr_id) const;

	virtual void calc_ci(AllData * ad, int tr_id) const;

	virtual void dump_std(AllData * ad) const;
	virtual void dump_lpn(AllData * ad) const;
	virtual void dump_std_evo(AllData * ad) const;
	virtual void dump_lpn_evo(AllData * ad) const;
};

Split * init_split_structure_dimer(AllData * ad);
Split * init_split_structure_dimer_deep(AllData * ad);
Split * init_split_structure_jcs(AllData * ad);
Split * init_split_structure_jcs_deep(AllData * ad);

void rk_right_part_dimer(AllData * ad, int sub_step, int tr_id, int th_id);
void rk_right_part_jcs(AllData * ad, int sub_step, int tr_id, int th_id);
void rk_int_dimer(AllData * ad, int tr_id, int th_id, double step);
void rk_int_jcs(AllData * ad, int tr_id, int th_id, double step);
void rk_step_dimer(AllData * ad, int tr_id, int th_id, double step);
void rk_step_jcs(AllData * ad, int tr_id, int th_id, double step);