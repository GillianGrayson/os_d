#include "outputter.h"

void DimerOutputBehavior::suffix_param(RunParam * rp, ConfigParam * cp, int precision) const
{
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

	string ext = extension();

	string suffix = qj + params.str() + ext;
	//string suffix = qj + ext;

	cp->fn_suffix = suffix;
}

void JCSOutputBehavior::suffix_param(RunParam * rp, ConfigParam * cp, int precision) const
{
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

	string ext = extension();

	string suffix = qj + params.str() + ext;
	//string suffix = qj + ext;

	cp->fn_suffix = suffix;
}

string suffix_qj(RunParam * rp, ConfigParam * cp, int precision)
{
	stringstream suffix;

	//suffix << "_config(" << rp->sys_id << "_" << rp->task_id << "_" << rp->prop_id << ")";
	suffix << "_rnd(" << cp->seed << "_" << cp->mns << ")";

	return suffix.str();
}

string extension()
{
	stringstream suffix;
	suffix << ".txt";
	return suffix.str();
}
