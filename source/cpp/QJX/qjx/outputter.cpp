#include "outputter.h"

void DimerIOuputBehavior::suffix_param(ConfigParam * cp, int precision) const
{
	string qj = suffix_qj(cp, precision);

	stringstream params;

	params << "_N(" << cp->params.find("N")->second << ")";

	params << "_diss("
		<< cp->params.find("diss_type")->second << "_"
		<< setprecision(precision) << cp->params.find("diss_gamma")->second << "_"
		<< setprecision(precision) << cp->params.find("diss_phase")->second << ")";

	params << "_drv("
		<< cp->params.find("drv_type")->second << "_"
		<< setprecision(precision) << cp->params.find("drv_ampl")->second << "_"
		<< setprecision(precision) << cp->params.find("drv_freq")->second << "_"
		<< setprecision(precision) << cp->params.find("drv_phase")->second << ")";

	params << "_start("
		<< cp->params.find("start_type")->second << "_"
		<< cp->params.find("start_state")->second << ")";

	params << "_prm("
		<< setprecision(precision) << cp->params.find("prm_E")->second << "_"
		<< setprecision(precision) << cp->params.find("prm_U")->second << "_"
		<< setprecision(precision) << cp->params.find("prm_J")->second << ")";

	string ext = extension();

	string suffix = qj + params.str() + ext;

	cp->fn_suffix = suffix;
}

string suffix_qj(ConfigParam * cp, int precision)
{
	stringstream suffix;

	suffix << "_qjrnd(" << cp->qj_seed << "_"
		<< cp->qj_mns << ")";

	return suffix.str();
}

string extension()
{
	stringstream suffix;
	suffix << ".txt";
	return suffix.str();
}
