#include "outputter.h"

string DimerIOuputBehavior::suffix(ConfigParam * cp, int precision) const
{
	stringstream suffix;

	suffix << "_N(" << cp->params.find("N")->second << ")";

	suffix << "_diss("
		<< cp->params.find("diss_type")->second << "_"
		<< setprecision(precision) << cp->params.find("diss_gamma")->second << "_"
		<< setprecision(precision) << cp->params.find("diss_phase")->second << ")";

	suffix << "_drv("
		<< cp->params.find("drv_type")->second << "_"
		<< setprecision(precision) << cp->params.find("drv_ampl")->second << "_"
		<< setprecision(precision) << cp->params.find("drv_freq")->second << "_"
		<< setprecision(precision) << cp->params.find("drv_phase")->second << ")";

	suffix << "_init("
		<< cp->params.find("init_type")->second << "_"
		<< cp->params.find("init_state")->second << ")";

	suffix << "_prm("
		<< setprecision(precision) << cp->params.find("prm_E")->second << "_"
		<< setprecision(precision) << cp->params.find("prm_U")->second << "_"
		<< setprecision(precision) << cp->params.find("prm_J")->second << ")";

	return suffix.str();
}
