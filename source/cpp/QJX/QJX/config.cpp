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
	} while (pos < str.length() && prev < str.length());

	return tokens;
}

void set_param(RunParam * rp, ConfigParam * cp, string str, string val)
{
	if (str.compare("sys_id") == 0)
	{
		rp->sys_id = atoi(val.c_str());
	}
	if (str.compare("task_id") == 0)
	{
		rp->task_id = atoi(val.c_str());
	}
	if (str.compare("prop_id") == 0)
	{
		rp->prop_id = atoi(val.c_str());
	}


	if (str.compare("is_debug") == 0)
	{
		rp->is_debug = atoi(val.c_str());
	}
	if (str.compare("is_pp") == 0)
	{
		rp->is_pp = atoi(val.c_str());
	}


	if (str.compare("init_fn") == 0)
	{
		rp->init_fn = val;
	}
	if (str.compare("path") == 0)
	{
		rp->path = val;
	}

	if (str.compare("num_threads") == 0)
	{
		rp->num_threads = atoi(val.c_str());
	}
	if (str.compare("num_trajectories") == 0)
	{
		cp->num_trajectories = atoi(val.c_str());
	}

	if (str.compare("seed") == 0)
	{
		cp->seed = atoi(val.c_str());
	}
	if (str.compare("mns") == 0)
	{
		cp->mns = atoi(val.c_str());
	}

	if (str.compare("num_tp_periods") == 0)
	{
		cp->num_tp_periods = atoi(val.c_str());
	}
	if (str.compare("num_obs_periods") == 0)
	{
		cp->num_obs_periods = atoi(val.c_str());
	}

	if (str.compare("ex_deep") == 0)
	{
		cp->ex_deep = atoi(val.c_str());
	}
	if (str.compare("rk_ns") == 0)
	{
		cp->rk_ns = atoi(val.c_str());
	}
}

void init_params(RunParam *rp, ConfigParam *cp, char * fn_config, char * fn_param)
{
	string line;
	ifstream config_file(fn_config);
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

	string key;
	double val;
	ifstream param_file(fn_param);
	if (param_file.is_open())
	{
		while (getline(param_file, line))
		{
			std::istringstream iss(line);
			iss >> key >> val;
			cp->params[key] = val;
		}

		param_file.close();
	}
	else
	{
		cout << "Unable to open file" << endl;
		cout << "Init with default params" << endl;
	}

	output_params(rp, cp);
}

void output_params(RunParam * rp, ConfigParam * cp)
{
	cout << "############# parameters #############" << endl;

	cout << "sys_id = " << rp->sys_id << endl;
	cout << "task_id = " << rp->task_id << endl;
	cout << "prop_id = " << rp->prop_id << endl;

	cout << "is_debug = " << rp->is_debug << endl;
	cout << "is_pp = " << rp->is_pp << endl;

	cout << "init_fn = " << rp->init_fn << endl;
	cout << "path = " << rp->path << endl;

	cout << "num_threads = " << rp->num_threads << endl;
	cout << "num_trajectories = " << cp->num_trajectories << endl;

	cout << "seed = " << cp->seed << endl;
	cout << "mns = " << cp->mns << endl;

	cout << "num_tp_periods = " << cp->num_tp_periods << endl;
	cout << "num_obs_periods = " << cp->num_obs_periods << endl;

	cout << "ex_deep = " << cp->ex_deep << endl;
	cout << "rk_ns = " << cp->rk_ns << endl;

	cout << endl;

	for (std::map<string, double>::iterator it = cp->params.begin(); it != cp->params.end(); ++it)
	{
		std::cout << it->first << ": " << it->second << endl;
	}
		
	cout << "######################################" << endl;
}