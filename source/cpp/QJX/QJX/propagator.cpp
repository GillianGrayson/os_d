#include "propagator.h"
#include "experiment.h"
#include "qj_proc.h"
#include "rk_proc.h"
#include "core.h"

void QJPropagateBehavior::init_prop_data(AllData * ad, CoreBehavior * cb) const
{
	cb->init_splits(ad);
}

void QJPropagateBehavior::free_prop_data(AllData * ad, CoreBehavior * cb) const
{
	cb->free_splits(ad);
}

void QJPropagateBehavior::init_prop_data_deep(AllData * ad, CoreBehavior * cb) const
{
	cb->init_splits_deep(ad);
}

void QJPropagateBehavior::free_prop_data_deep(AllData * ad, CoreBehavior * cb) const
{
	cb->free_splits_deep(ad);
}

void QJPropagateBehavior::one_period(AllData * ad, CoreBehavior * cb, int tr_id, int th_id, int period_id) const
{
	cb->ex_period(ad, tr_id, th_id, period_id);
}

void QJPropagateBehavior::one_period_trp_deep(AllData * ad, CoreBehavior * cb, int tr_id, int th_id, int period_id) const
{
	cb->ex_period_trp_deep(ad, tr_id, th_id, period_id);
}

void QJPropagateBehavior::one_period_obs_deep(AllData * ad, CoreBehavior * cb, int tr_id, int th_id, int period_id) const
{
	cb->ex_period_trp_deep(ad, tr_id, th_id, period_id);
}

void QJPropagateBehavior::one_period_obs_deep_lpn(AllData * ad, CoreBehavior * cb, int tr_id, int th_id, int period_id) const
{
	cb->ex_period_obs_deep_lpn(ad, tr_id, th_id, period_id);
}

void QJPropagateBehavior::one_period_obs_deep_cd(AllData * ad, CoreBehavior * cb, int tr_id, int th_id, int period_id) const
{
	cb->ex_period_obs_deep_cd(ad, tr_id, th_id, period_id);
}

void QJPropagateBehavior::one_period_obs_deep_sigma(AllData * ad, CoreBehavior * cb, int tr_id, int th_id, int period_id) const
{
	cb->ex_period_obs_deep_sigma(ad, tr_id, th_id, period_id);
}

void RKPropagateBehavior::init_prop_data(AllData * ad, CoreBehavior * cb) const
{
	init_rk(ad);
}

void RKPropagateBehavior::free_prop_data(AllData * ad, CoreBehavior * cb) const
{
	free_rk(ad);
}

void RKPropagateBehavior::init_prop_data_deep(AllData * ad, CoreBehavior * cb) const
{
	init_rk_deep(ad);
}

void RKPropagateBehavior::free_prop_data_deep(AllData * ad, CoreBehavior * cb) const
{
	free_rk_deep(ad);
}

void RKPropagateBehavior::one_period(AllData * ad, CoreBehavior * cb, int tr_id, int th_id, int period_id) const
{
	cb->rk_period(ad, tr_id, th_id, period_id);
}

void RKPropagateBehavior::one_period_trp_deep(AllData * ad, CoreBehavior * cb, int tr_id, int th_id, int period_id) const
{
	cb->rk_period_trp_deep(ad, tr_id, th_id, period_id);
}

void RKPropagateBehavior::one_period_obs_deep(AllData * ad, CoreBehavior * cb, int tr_id, int th_id, int period_id) const
{
	cb->rk_period_obs_deep(ad, tr_id, th_id, period_id);
}

void RKPropagateBehavior::one_period_obs_deep_lpn(AllData * ad, CoreBehavior * cb, int tr_id, int th_id, int period_id) const
{
	cb->rk_period_obs_deep_lpn(ad, tr_id, th_id, period_id);
}

void RKPropagateBehavior::one_period_obs_deep_cd(AllData * ad, CoreBehavior * cb, int tr_id, int th_id, int period_id) const
{
	cb->rk_period_obs_deep_cd(ad, tr_id, th_id, period_id);
}

void RKPropagateBehavior::one_period_obs_deep_sigma(AllData * ad, CoreBehavior * cb, int tr_id, int th_id, int period_id) const
{
	cb->rk_period_obs_deep_sigma(ad, tr_id, th_id, period_id);
}
