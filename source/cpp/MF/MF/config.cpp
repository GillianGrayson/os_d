#include "config.h"

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
	} 
	while (pos < str.length() && prev < str.length());

	return tokens;
}

void set_param(RunParam &rp, ConfigParam &cp, string str, string val)
{
	if (str.compare("task") == 0)
	{
		rp.task = atoi(val.c_str());
	}
	if (str.compare("U_start") == 0)
	{
		rp.U_start = atof(val.c_str());
	}
	if (str.compare("U_shift") == 0)
	{
		rp.U_shift = atof(val.c_str());
	}
	if (str.compare("U_num") == 0)
	{
		rp.U_num = atoi(val.c_str());
	}
	if (str.compare("seed_start") == 0)
	{
		rp.seed_start = atoi(val.c_str());
	}
	if (str.compare("seed_num") == 0)
	{
		rp.seed_num = atoi(val.c_str());
	}
	if (str.compare("max_num_seeds") == 0)
	{
		rp.max_num_seeds = atoi(val.c_str());
	}
	if (str.compare("path") == 0)
	{
		rp.path = val;
	}


	if (str.compare("mt") == 0)
	{
		cp.mt = atoi(val.c_str());
	}
	if (str.compare("num_steps") == 0)
	{
		cp.num_steps = atoi(val.c_str());
	}
	if (str.compare("npt") == 0)
	{
		cp.npt = atoi(val.c_str());
	}
	if (str.compare("np") == 0)
	{
		cp.np = atoi(val.c_str());
	}
	if (str.compare("E") == 0)
	{
		cp.E = atof(val.c_str());
	}
	if (str.compare("A") == 0)
	{
		cp.A = atof(val.c_str());
	}
	if (str.compare("omega") == 0)
	{
		cp.omega = atof(val.c_str());
		cp.T = 2.0 * PI / cp.omega;
	}
	if (str.compare("phase") == 0)
	{
		cp.phase = atof(val.c_str());
	}
	if (str.compare("gamma") == 0)
	{
		cp.gamma = atof(val.c_str());
	}
	if (str.compare("J") == 0)
	{
		cp.J = atof(val.c_str());
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
		cout << "unable to open file" << endl;
		cout << "init with default params" << endl;
	}

	output_params(rp, cp);
}

void output_params(RunParam &rp, ConfigParam &cp)
{
	cout << "############# parameters #############" << endl;

	cout << "task = " << rp.task << endl;

	cout << "U_start = " << rp.U_start << endl;
	cout << "U_shift = " << rp.U_shift << endl;
	cout << "U_num = " << rp.U_num << endl;

	cout << "seed_start = " << rp.seed_start << endl;
	cout << "seed_num = " << rp.seed_num << endl;
	cout << "seed_num = " << rp.seed_num << endl;
	cout << "max_num_seeds = " << rp.max_num_seeds << endl;

	cout << "path = " << rp.path << endl;


	cout << "mt = " << cp.mt << endl;
	cout << "num_steps = " << cp.num_steps << endl;
	cout << "npt = " << cp.npt << endl;
	cout << "np = " << cp.np << endl;
	cout << "E = " << cp.E << endl;
	cout << "A = " << cp.A << endl;
	cout << "omega = " << cp.omega << endl;
	cout << "phase = " << cp.phase << endl;
	cout << "gamma = " << cp.gamma << endl;
	cout << "J = " << cp.J << endl;

	cout << "######################################" << endl;
}