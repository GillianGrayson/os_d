#include "outputter.h"

void DimerIOuputBehavior::suffix_param(ConfigParam * cp, int precision) const
{
	string qj = suffix_qj(cp, precision);

	stringstream params;

	params << "_N(" << int(cp->params.find("N")->second) << ")";

	params << "_diss("
		<< int(cp->params.find("diss_type")->second) << "_"
		<< setprecision(precision) << fixed << double(cp->params.find("diss_gamma")->second) << "_"
		<< setprecision(precision) << fixed << double(cp->params.find("diss_phase")->second) << ")";

	params << "_drv("
		<< int(cp->params.find("drv_type")->second) << "_"
		<< setprecision(precision) << fixed << double(cp->params.find("drv_ampl")->second) << "_"
		<< setprecision(precision) << fixed << double(cp->params.find("drv_freq")->second) << "_"
		<< setprecision(precision) << fixed << double(cp->params.find("drv_phase")->second) << ")";

	params << "_prm("
		<< setprecision(precision) << fixed << double(cp->params.find("prm_E")->second) << "_"
		<< setprecision(precision) << fixed << double(cp->params.find("prm_U")->second) << "_"
		<< setprecision(precision) << fixed << double(cp->params.find("prm_J")->second) << ")";

	params << "_start("
		<< int(cp->params.find("start_type")->second) << "_"
		<< int(cp->params.find("start_state")->second) << ")";

	string ext = extension();

	string suffix = qj + params.str() + ext;

	cp->fn_suffix = suffix;
}

string suffix_qj(ConfigParam * cp, int precision)
{
	stringstream suffix;

	suffix << "_rnd(" << cp->seed << "_" << cp->mns << ")";

	return suffix.str();
}

string extension()
{
	stringstream suffix;
	suffix << ".txt";
	return suffix.str();
}
