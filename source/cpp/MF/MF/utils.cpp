#include "utils.h"

void error(const string& err, const char* func, const char* file, int line)
{
	cout << "Error: " << err << " in " << func << " (" << file << ", line: " << line << ")" << endl;
	exit(EXIT_FAILURE);
}

void print_int_array(int * data, int N)
{
	cout << endl;
	for (int i = 0; i < N; i++)
	{
		cout << i << ": " << data[i] << endl;
	}
}

string file_name_suffix(RunParam &rp, ConfigParam &cp, int precision)
{
	stringstream fns;
	fns << "_t(" << rp.task << ")";
	if (rp.task >= 3)
	{
		fns << "_eps(" << setprecision(10) << fixed << cp.cd_eps << ")";
		fns << "_m(" << cp.cd_dim << ")";
	}
	fns << "_mt(" << cp.mt << ")";
	fns << "_omega(" << setprecision(precision) << fixed << cp.omega << ")";
	fns << "_phase(" << setprecision(precision) << fixed << cp.phase << ")";
	fns << "_g(" << setprecision(precision) << fixed << cp.gamma << ")";
	fns << "_J(" << setprecision(precision) << fixed << cp.J << ")";
	fns << "_E(" << setprecision(precision) << fixed << cp.E << ")";
	fns << "_A(" << setprecision(precision) << fixed << cp.A << ")";
	fns << "_U(" << setprecision(precision) << fixed << cp.U << ")";
	fns << "_seed(" << cp.seed << ")";
	fns << ".txt";

	return fns.str();
}

void write_double_data(string file_name, double * data, int size, int precision, bool append)
{
	if (append)
	{
		ofstream ofs = ofstream(file_name, ios::app);

		if (ofs.is_open())
		{
			ofs << setprecision(precision) << scientific;
			for (int i = 0; i < size; i++)
			{
				ofs << data[i] << endl;
			}

			ofs.close();
		}
		else
		{
			stringstream msg;
			msg << "unable to open file:" << endl << file_name << endl;
			Error(msg.str());
		}
	}
	else
	{
		ofstream ofs = ofstream(file_name);

		if (ofs.is_open())
		{
			ofs << setprecision(precision) << scientific;
			for (int i = 0; i < size; i++)
			{
				ofs << data[i] << endl;
			}

			ofs.close();
		}
		else
		{
			stringstream msg;
			msg << "unable to open file:" << endl << file_name << endl;
			Error(msg.str());
		}
	}
}

void write_2d_double_data(string file_name, double ** data, int size_1, int size_2, int precision, bool append)
{
	if (append)
	{
		ofstream ofs = ofstream(file_name, ios::app);

		if (ofs.is_open())
		{
			ofs << setprecision(precision) << scientific;

			for (int row = 0; row < size_2; row++)
			{
				for (int col = 0; col < size_1; col++)
				{
					ofs << data[col][row] << " ";
				}
				ofs << endl;
			}

			ofs.close();
		}
		else
		{
			stringstream msg;
			msg << "Unable to open file:" << endl << file_name << endl;
			Error(msg.str());
		}
	}
	else
	{
		ofstream ofs = ofstream(file_name);

		if (ofs.is_open())
		{
			ofs << setprecision(precision) << scientific;

			for (int row = 0; row < size_2; row++)
			{
				for (int col = 0; col < size_1; col++)
				{
					ofs << data[col][row] << " ";
				}
				ofs << endl;
			}

			ofs.close();
		}
		else
		{
			stringstream msg;
			msg << "Unable to open file:" << endl << file_name << endl;
			Error(msg.str());
		}
	}
}

void write_int_data(string file_name, int * data, int size, bool append)
{
	if (append)
	{
		ofstream ofs = ofstream(file_name, ios::app);

		if (ofs.is_open())
		{
			for (int i = 0; i < size; i++)
			{
				ofs << data[i] << endl;
			}

			ofs.close();
		}
		else
		{
			stringstream msg;
			msg << "unable to open file:" << endl << file_name << endl;
			Error(msg.str());
		}
	}
	else
	{
		ofstream ofs = ofstream(file_name);

		if (ofs.is_open())
		{
			for (int i = 0; i < size; i++)
			{
				ofs << data[i] << endl;
			}

			ofs.close();
		}
		else
		{
			stringstream msg;
			msg << "unable to open file:" << endl << file_name << endl;
			Error(msg.str());
		}
	}
}


void write_complex_data(string file_name, MKL_Complex16 * data, int size, int precision, bool append)
{
	if (append)
	{
		ofstream ofs = ofstream(file_name, ios::app);

		if (ofs.is_open())
		{
			ofs << setprecision(precision) << scientific;
			for (int i = 0; i < size; i++)
			{
				ofs << data[i].real << " " << data[i].imag << endl;
			}

			ofs.close();
		}
		else
		{
			stringstream msg;
			msg << "unable to open file:" << endl << file_name << endl;
			Error(msg.str());
		}
	}
	else
	{
		ofstream ofs = ofstream(file_name);

		if (ofs.is_open())
		{
			ofs << setprecision(precision) << scientific;
			for (int i = 0; i < size; i++)
			{
				ofs << data[i].real << " " << data[i].imag << endl;
			}

			ofs.close();
		}
		else
		{
			stringstream msg;
			msg << "unable to open file:" << endl << file_name << endl;
			Error(msg.str());
		}
	}
}

vector<int> sort_doubles_with_order(vector<double> &v)
{
	vector<int> order(v.size());

	iota(order.begin(), order.end(), 0);

	sort(order.begin(), order.end(),
		[&v](size_t i1, size_t i2) {return v[i1] < v[i2]; });

	return order;
}

double calc_norm(double * vec, int size)
{
	double sum = 0.0;
	for (int i = 0; i < size; i++)
	{
		sum += vec[i] * vec[i];
	}

	return sqrt(sum);
}

void delete_data(double * data)
{
	if (data)
	{
		delete[] data;
		data = NULL;
	}
}