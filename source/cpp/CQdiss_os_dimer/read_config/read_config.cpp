#include "read_config.h"


vector<string> split(const string& str, const string& delim)
{
	vector<string> tokens;
	size_t prev = 0, pos = 0;
	do
	{
		pos = str.find(delim, prev);

		if (pos == string::npos)
		{
			pos = str.length();
		}

		string token = str.substr(prev, pos - prev);

		if (!token.empty())
		{
			tokens.push_back(token);
		}

		prev = pos + delim.length();
	} while (pos < str.length() && prev < str.length());

	return tokens;
}

void set_param(RunParam &rp, ConfigParam &cp, string str, string val)
{
	if (str.compare("task") == 0)
	{
		rp.task = atoi(val.c_str());
	}

	if (str.compare("debug") == 0)
	{
		rp.debug = atoi(val.c_str());
	}

	if (str.compare("issmtx") == 0)
	{
		rp.issmtx = atoi(val.c_str());
	}

	if (str.compare("ipp") == 0)
	{
		rp.ipp = atoi(val.c_str());
	}

	if (str.compare("N") == 0)
	{
		cp.N = atoi(val.c_str());
	}

	if (str.compare("diss_type") == 0)
	{
		cp.diss_type = atoi(val.c_str());
	}
	if (str.compare("diss_gamma") == 0)
	{
		cp.diss_gamma = atof(val.c_str());
	}
	if (str.compare("diss_phase") == 0)
	{
		cp.diss_phase = atof(val.c_str());
	}


	if (str.compare("drv_type") == 0)
	{
		cp.drv_type = atoi(val.c_str());
	}
	if (str.compare("drv_ampl") == 0)
	{
		cp.drv_ampl = atof(val.c_str());
	}
	if (str.compare("drv_freq") == 0)
	{
		cp.drv_freq = atof(val.c_str());
	}
	if (str.compare("drv_phase") == 0)
	{
		cp.drv_phase = atof(val.c_str());
	}

	if (str.compare("prm_E") == 0)
	{
		cp.prm_E = atof(val.c_str());
	}
	if (str.compare("prm_U") == 0)
	{
		cp.prm_U = atof(val.c_str());
	}
	if (str.compare("prm_J") == 0)
	{
		cp.prm_J = atof(val.c_str());
	}


	if (str.compare("num_steps") == 0)
	{
		cp.num_steps = atoi(val.c_str());
	}
	if (str.compare("num_periods_trans") == 0)
	{
		cp.num_periods_trans = atoi(val.c_str());
	}
	if (str.compare("num_periods_obser") == 0)
	{
		cp.num_periods_obser = atoi(val.c_str());
	}

	if (str.compare("seed") == 0)
	{
		cp.seed = atoi(val.c_str());
	}
	if (str.compare("max_num_seeds") == 0)
	{
		cp.max_num_seeds = atoi(val.c_str());
	}

	if (str.compare("int_ist") == 0)
	{
		cp.int_ist = atoi(val.c_str());
	}
	if (str.compare("int_isi") == 0)
	{
		cp.int_isi = atoi(val.c_str());
	}
	if (str.compare("int_dt") == 0)
	{
		cp.int_dt = atoi(val.c_str());
	}
	if (str.compare("int_dn") == 0)
	{
		cp.int_dn = atoi(val.c_str());
	}
}

void init_params(RunParam &rp, ConfigParam &cp, char * file_name)
{
	string line;
	ifstream config_file(file_name);
	if (config_file.is_open())
	{
		while (getline(config_file, line))
		{
			vector<string> tokens = split(line, " ");

			if (tokens.size() == 2)
			{
				set_param(rp, cp, tokens[0], tokens[1]);
			}
		}
		config_file.close();
	}
	else
	{
		cout << "Unable to open file" << endl;
		cout << "Init with default params" << endl;
	}

	output_params(rp, cp);
}

void output_params(RunParam &rp, ConfigParam &cp)
{
	cout << "############# parameters #############" << endl;
	cout << "task = " << rp.task << endl;
	cout << "debug = " << rp.debug << endl;
	cout << "issmtx = " << rp.issmtx << endl;
	cout << "ipp = " << rp.ipp << endl;

	cout << "N = " << cp.N << endl;

	cout << "diss_type = " << cp.diss_type << endl;
	cout << "diss_gamma = " << cp.diss_gamma << endl;
	cout << "diss_phase = " << cp.diss_phase << endl;

	cout << "drv_type = " << cp.drv_type << endl;
	cout << "drv_ampl = " << cp.drv_ampl << endl;
	cout << "drv_freq = " << cp.drv_freq << endl;
	cout << "drv_phase = " << cp.drv_phase << endl;

	cout << "prm_E = " << cp.prm_E << endl;
	cout << "prm_U = " << cp.prm_U << endl;
	cout << "prm_J = " << cp.prm_J << endl;

	cout << "num_steps = " << cp.num_steps << endl;
	cout << "num_periods_trans = " << cp.num_periods_trans << endl;
	cout << "num_periods_obser = " << cp.num_periods_obser << endl;

	cout << "seed = " << cp.seed << endl;
	cout << "max_num_seeds = " << cp.max_num_seeds << endl;

	cout << "int_ist = " << cp.int_ist << endl;
	cout << "int_isi = " << cp.int_isi << endl;
	cout << "int_dt = " << cp.int_dt << endl;
	cout << "int_dn = " << cp.int_dn << endl;

	cout << "######################################" << endl;
}