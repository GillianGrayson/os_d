#include "outputter.h"

void DimerOutputBehavior::suffix_param(RunParam * rp, ConfigParam * cp, int precision) const
{
	string setup = suffix_setup(rp);
	string qj = suffix_qj(rp, cp, precision);

	stringstream params;

	params << "_N(" << int(cp->params.find("N")->second) << ")";

	params << "_diss("
		<< int(cp->params.find("diss_type")->second) << "_"
		<< setprecision(precision) << fixed << double(cp->params.find("diss_gamma")->second) << "_"
		<< setprecision(precision) << fixed << double(cp->params.find("diss_phase")->second) << ")";

	params << "_drv("
		<< int(cp->params.find("dimer_drv_type")->second) << "_"
		<< setprecision(precision) << fixed << double(cp->params.find("dimer_drv_ampl")->second) << "_"
		<< setprecision(precision) << fixed << double(cp->params.find("dimer_drv_freq")->second) << "_"
		<< setprecision(precision) << fixed << double(cp->params.find("dimer_drv_phase")->second) << ")";

	params << "_prm("
		<< setprecision(precision) << fixed << double(cp->params.find("dimer_prm_E")->second) << "_"
		<< setprecision(precision) << fixed << double(cp->params.find("dimer_prm_U")->second) << "_"
		<< setprecision(precision) << fixed << double(cp->params.find("dimer_prm_J")->second) << ")";

	params << "_start("
		<< int(cp->params.find("start_type")->second) << "_"
		<< int(cp->params.find("start_state")->second) << ")";

	string suffix = setup + qj + params.str();

	if (rp->task_id == LPN_MULT_TASK_ID || rp->task_id == LPN_MULT_DEEP_TASK_ID)
	{
		string lpn = suffix_lpn(rp, cp);
		suffix += lpn;
	}

	string ext = extension();
	suffix += ext;

	cp->fn_suffix = suffix;
}

void DimerSyncOutputBehavior::suffix_param(RunParam* rp, ConfigParam* cp, int precision) const
{
	string setup = suffix_setup(rp);
	string qj = suffix_qj(rp, cp, precision);

	stringstream params;

	params << "_N(" << int(cp->params.find("dimersync_N")->second) << ")";

	params << "_diss("
		<< int(cp->params.find("diss_type")->second) << "_"
		<< setprecision(precision) << fixed << double(cp->params.find("diss_gamma")->second) << "_"
		<< setprecision(precision) << fixed << double(cp->params.find("diss_phase")->second) << ")";

	params << "_drv1("
		<< setprecision(precision) << fixed << double(cp->params.find("dimersync_drv_ampl_1")->second) << "_"
		<< setprecision(precision) << fixed << double(cp->params.find("dimersync_drv_freq_1")->second) << "_"
		<< setprecision(precision) << fixed << double(cp->params.find("dimersync_drv_phase_1")->second) << ")";

	params << "_drv2("
		<< setprecision(precision) << fixed << double(cp->params.find("dimersync_drv_ampl_2")->second) << "_"
		<< setprecision(precision) << fixed << double(cp->params.find("dimersync_drv_freq_2")->second) << "_"
		<< setprecision(precision) << fixed << double(cp->params.find("dimersync_drv_phase_2")->second) << ")";

	params << "_T(" << int(cp->params.find("dimersync_main_period")->second) << ")";

	params << "_prm("
		<< setprecision(precision) << fixed << double(cp->params.find("dimersync_prm_E")->second) << "_"
		<< setprecision(precision) << fixed << double(cp->params.find("dimersync_prm_U")->second) << "_"
		<< setprecision(precision) << fixed << double(cp->params.find("dimersync_prm_J")->second) << ")";

	params << "_start("
		<< int(cp->params.find("start_type")->second) << "_"
		<< int(cp->params.find("start_state")->second) << ")";

	string suffix = setup + qj + params.str();

	if (rp->task_id == LPN_MULT_TASK_ID || rp->task_id == LPN_MULT_DEEP_TASK_ID)
	{
		string lpn = suffix_lpn(rp, cp);
		suffix += lpn;
	}

	string ext = extension();
	suffix += ext;

	cp->fn_suffix = suffix;
}

void JCSOutputBehavior::suffix_param(RunParam * rp, ConfigParam * cp, int precision) const
{
	string setup = suffix_setup(rp);
	string qj = suffix_qj(rp, cp, precision);

	stringstream params;

	params << "_N(" << int(cp->params.find("N")->second) << ")";

	params << "_diss("
		<< int(cp->params.find("diss_type")->second) << "_"
		<< setprecision(precision) << fixed << double(cp->params.find("diss_gamma")->second) << "_"
		<< setprecision(precision) << fixed << double(cp->params.find("diss_phase")->second) << ")";

	params << "_drv("
		<< setprecision(precision) << fixed << double(cp->params.find("jcs_drv_part_1")->second) << "_"
		<< setprecision(precision) << fixed << double(cp->params.find("jcs_drv_part_2")->second) << "_"
		<< setprecision(precision) << fixed << double(cp->params.find("jcs_drv_ampl")->second) << ")";

	params << "_prm("
		<< setprecision(precision) << fixed << double(cp->params.find("jcs_prm_alpha")->second) << ")";

	params << "_start("
		<< int(cp->params.find("start_type")->second) << "_"
		<< int(cp->params.find("start_state")->second) << ")";

	string suffix = setup + qj + params.str();

	if (rp->task_id == LPN_MULT_TASK_ID || rp->task_id == LPN_MULT_DEEP_TASK_ID)
	{
		string lpn = suffix_lpn(rp, cp);
		suffix += lpn;
	}

	string ext = extension();
	suffix += ext;

	cp->fn_suffix = suffix;
}

void PSOutputBehavior::suffix_param(RunParam * rp, ConfigParam * cp, int precision) const
{
	string setup = suffix_setup(rp);
	string qj = suffix_qj(rp, cp, precision);

	stringstream params;

	params << "_s(" << int(cp->params.find("ps_num_spins")->second) << ")";
	params << "_nps(" << int(cp->params.find("ps_num_photons_states")->second) << ")";

	params << "_diss("
		<< int(cp->params.find("diss_type")->second) << "_"
		<< setprecision(precision) << fixed << double(cp->params.find("ps_diss_w")->second) << ")";

	params << "_drv("
		<< setprecision(precision) << fixed << double(cp->params.find("ps_drv_part_1")->second) << "_"
		<< setprecision(precision) << fixed << double(cp->params.find("ps_drv_part_2")->second) << "_"
		<< setprecision(precision) << fixed << double(cp->params.find("ps_drv_ampl")->second) << ")";

	params << "_prm("
		<< setprecision(precision) << fixed << double(cp->params.find("ps_prm_alpha")->second) << "_"
		<< setprecision(precision) << fixed << double(cp->params.find("ps_prm_d")->second) << "_"
		<< setprecision(precision) << fixed << double(cp->params.find("ps_prm_g")->second) << ")";

	params << "_start("
		<< int(cp->params.find("start_type")->second) << "_"
		<< int(cp->params.find("start_state")->second) << ")";

	string suffix = setup + qj + params.str();

	if (rp->task_id == LPN_MULT_TASK_ID || rp->task_id == LPN_MULT_DEEP_TASK_ID)
	{
		string lpn = suffix_lpn(rp, cp);
		suffix += lpn;
	}

	string ext = extension();
	suffix += ext;

	cp->fn_suffix = suffix;
}

void MBLOutputBehavior::suffix_param(RunParam * rp, ConfigParam * cp, int precision) const
{
	string setup = suffix_setup(rp);
	string qj = suffix_qj(rp, cp, precision);

	stringstream params;

	params << "_Nc(" << int(cp->params.find("mbl_Nc")->second) << ")";

	params << "_rnd("
		<< int(cp->params.find("mbl_seed")->second) << "_"
		<< int(cp->params.find("mbl_mns")->second) << ")";

	params << "_diss("
		<< int(cp->params.find("diss_type")->second) << "_"
		<< setprecision(precision) << fixed << double(cp->params.find("diss_phase")->second) << "_"
		<< setprecision(precision) << fixed << double(cp->params.find("diss_gamma")->second) << ")";

	params << "_prm("
		<< setprecision(precision) << fixed << double(cp->params.find("mbl_prm_W")->second) << "_"
		<< setprecision(precision) << fixed << double(cp->params.find("mbl_prm_U")->second) << "_"
		<< setprecision(precision) << fixed << double(cp->params.find("mbl_prm_J")->second) << ")";

	params << "_start("
		<< int(cp->params.find("start_type")->second) << "_"
		<< int(cp->params.find("start_state")->second) << ")";

	string suffix = setup + qj + params.str();

	if (rp->task_id == LPN_MULT_TASK_ID || rp->task_id == LPN_MULT_DEEP_TASK_ID)
	{
		string lpn = suffix_lpn(rp, cp);
		suffix += lpn;
	}

	string ext = extension();
	suffix += ext;

	cp->fn_suffix = suffix;
}

void LndHamOutputBehavior::suffix_param(RunParam* rp, ConfigParam* cp, int precision) const
{
	string setup = suffix_setup(rp);
	string qj = suffix_qj(rp, cp, precision);

	stringstream params;

	params << "_N(" << int(cp->params.find("lndham_N")->second) << ")";
	params << "_rnd(" << int(cp->params.find("lndham_seed")->second) << ")";
	params << "_alpha(" << setprecision(precision) << fixed << double(cp->params.find("lndham_alpha")->second) << ")";
	params << "_T(" << setprecision(precision) << fixed << double(cp->params.find("lndham_T")->second) << ")";

	string suffix = setup + qj + params.str();

	if (rp->task_id == LPN_MULT_TASK_ID || rp->task_id == LPN_MULT_DEEP_TASK_ID)
	{
		string lpn = suffix_lpn(rp, cp);
		suffix += lpn;
	}

	string ext = extension();
	suffix += ext;

	cp->fn_suffix = suffix;
}

void IntegrableOutputBehavior::suffix_param(RunParam* rp, ConfigParam* cp, int precision) const
{
	string setup = suffix_setup(rp);
	string qj = suffix_qj(rp, cp, precision);

	stringstream params;

	params << "_N(" << int(cp->params.find("integrable_N")->second) << ")";
	params << "_seed(" << int(cp->params.find("integrable_seed")->second) << ")";
	params << "_tau(" << int(cp->params.find("integrable_tau")->second) << ")";
	params << "_k(" << int(cp->params.find("integrable_k")->second) << ")";
	params << "_T(" << setprecision(precision) << fixed << double(cp->params.find("integrable_T")->second) << ")";

	string suffix = setup + qj + params.str();

	if (rp->task_id == LPN_MULT_TASK_ID || rp->task_id == LPN_MULT_DEEP_TASK_ID)
	{
		string lpn = suffix_lpn(rp, cp);
		suffix += lpn;
	}

	string ext = extension();
	suffix += ext;

	cp->fn_suffix = suffix;
}

void Floq2SpinsOutputBehavior::suffix_param(RunParam* rp, ConfigParam* cp, int precision) const
{
	string setup = suffix_setup(rp);
	string qj = suffix_qj(rp, cp, precision);

	stringstream params;

	params << "_ampl(" << setprecision(precision) << fixed << double(cp->params.find("flq2sp_ampl")->second) << ")";
	params << "_freq(" << setprecision(precision) << fixed << double(cp->params.find("flq2sp_freq")->second) << ")";
	params << "_D1(" << setprecision(precision) << fixed << double(cp->params.find("flq2sp_Delta_1")->second) << ")";
	params << "_D2(" << setprecision(precision) << fixed << double(cp->params.find("flq2sp_Delta_2")->second) << ")";
	params << "_J(" << setprecision(precision) << fixed << double(cp->params.find("flq2sp_J")->second) << ")";
	params << "_gamma(" << setprecision(precision) << fixed << double(cp->params.find("flq2sp_gamma")->second) << ")";

	string suffix = setup + qj + params.str();

	if (rp->task_id == LPN_MULT_TASK_ID || rp->task_id == LPN_MULT_DEEP_TASK_ID)
	{
		string lpn = suffix_lpn(rp, cp);
		suffix += lpn;
	}

	string ext = extension();
	suffix += ext;

	cp->fn_suffix = suffix;
}

void FloqManySpinsOutputBehavior::suffix_param(RunParam* rp, ConfigParam* cp, int precision) const
{
	string setup = suffix_setup(rp);
	string qj = suffix_qj(rp, cp, precision);

	stringstream params;
	params << "_n(" << int(cp->params.find("flqnsp_n")->second) << ")";
	params << "_ampl(" << setprecision(precision) << fixed << double(cp->params.find("flqnsp_ampl")->second) << ")";
	params << "_freq(" << setprecision(precision) << fixed << double(cp->params.find("flqnsp_freq")->second) << ")";
	params << "_D(" << setprecision(precision) << fixed << double(cp->params.find("flqnsp_Delta")->second) << ")";
	params << "_J(" << setprecision(precision) << fixed << double(cp->params.find("flqnsp_J")->second) << ")";
	params << "_gamma(" << setprecision(precision) << fixed << double(cp->params.find("flqnsp_gamma")->second) << ")";

	string suffix = setup + qj + params.str();

	if (rp->task_id == LPN_MULT_TASK_ID || rp->task_id == LPN_MULT_DEEP_TASK_ID)
	{
		string lpn = suffix_lpn(rp, cp);
		suffix += lpn;
	}

	string ext = extension();
	suffix += ext;

	cp->fn_suffix = suffix;
}

string suffix_qj(RunParam * rp, ConfigParam * cp, int precision)
{
	stringstream suffix;

	suffix << "_rnd(" << cp->seed << "_" << cp->mns << ")";

	return suffix.str();
}

string suffix_setup(RunParam * rp)
{
	stringstream suffix;

	suffix << "_setup(" << rp->sys_id << "_" << rp->task_id << "_" << rp->prop_id << ")";

	return suffix.str();
}

string suffix_lpn(RunParam * rp, ConfigParam * cp)
{
	stringstream suffix;

	suffix << "_lpn(" << double(cp->params.find("lpn_type")->second) << "_"
		<< setprecision(4) << fixed << log10(double(cp->params.find("lpn_delta_s")->second)) << "_"
		<< setprecision(4) << fixed << log10(double(cp->params.find("lpn_delta_f_high")->second)) << "_"
		<< setprecision(4) << fixed << log10(double(cp->params.find("lpn_delta_f_low")->second)) << ")";

	return suffix.str();
}

string extension()
{
	stringstream suffix;
	suffix << ".txt";
	return suffix.str();
}
