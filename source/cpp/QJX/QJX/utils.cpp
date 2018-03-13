#include "utils.h"

void error(const string& err, const char* func, const char* file, int line)
{
	cout << "Error: " << err << " in " << func << " (" << file << ", line: " << line << ")" << endl;
	exit(EXIT_FAILURE);
}

int n_choose_k(int n, int k)
{
	int res = 1;

	if (k > n - k)
	{
		k = n - k;
	}

	for (int i = 0; i < k; ++i)
	{
		res *= (n - i);
		res /= (i + 1);
	}

	return res;
}

int bit_count(int value)
{
	int count = 0;

	while (value > 0)
	{
		if ((value & 1) == 1)
		{
			count++;
		}
		value >>= 1;
	}

	return count;
}

int bit_at(int value, int position)
{
	if ((value >> position) & 1)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

void print_int_array(int * data, int N)
{
	cout << endl;
	for (int i = 0; i < N; i++)
	{
		cout << i << ": " << data[i] << endl;
	}
}

void save_double_data(string file_name, double * data, int size, int precision, bool append)
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
			for (int i = 0; i < size; i++)
			{
				ofs << data[i] << endl;
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

void save_2d_double_data(string file_name, double * data, int num_rows, int num_cols, int precision, bool append)
{
	if (append)
	{
		ofstream ofs = ofstream(file_name, ios::app);

		if (ofs.is_open())
		{
			ofs << setprecision(precision) << scientific;

			for (int row = 0; row < num_rows; row++)
			{
				for (int col = 0; col < num_cols; col++)
				{
					int index = row * num_cols + col;
					ofs << data[index] << " ";
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

			for (int row = 0; row < num_rows; row++)
			{
				for (int col = 0; col < num_cols; col++)
				{
					int index = row * num_cols + col;
					ofs << data[index] << " ";
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

void save_2d_inv_double_data(string file_name, double * data, int num_rows, int num_cols, int precision, bool append)
{
	if (append)
	{
		ofstream ofs = ofstream(file_name, ios::app);

		if (ofs.is_open())
		{
			ofs << setprecision(precision) << scientific;

			for (int row = 0; row < num_rows; row++)
			{
				for (int col = 0; col < num_cols; col++)
				{
					int index = col * num_cols + row;
					ofs << data[index] << " ";
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

			for (int row = 0; row < num_rows; row++)
			{
				for (int col = 0; col < num_cols; col++)
				{
					int index = col * num_rows + row;
					ofs << data[index] << " ";
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

void save_double_vector(string file_name, vector<double> vec, int precision, bool append)
{
	if (append)
	{
		ofstream ofs = ofstream(file_name, ios::app);

		if (ofs.is_open())
		{
			ofs << setprecision(precision) << scientific;

			for (int st_id = 0; st_id < vec.size(); st_id++)
			{
				ofs << vec[st_id] << endl;
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

			for (int st_id = 0; st_id < vec.size(); st_id++)
			{
				ofs << vec[st_id] << endl;
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

void save_int_data(string file_name, int * data, int size, bool append)
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
			msg << "Unable to open file:" << endl << file_name << endl;
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
			msg << "Unable to open file:" << endl << file_name << endl;
			Error(msg.str());
		}
	}
}

void save_complex_data(string file_name, MKL_Complex16 * data, int size, int precision, bool append)
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
			for (int i = 0; i < size; i++)
			{
				ofs << data[i].real << " " << data[i].imag << endl;
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

vector<int> convert_int_to_vector_of_bits(int x, int size)
{
	vector<int> res;

	int id = 0;

	while (id < size)
	{
		if (x & 1)
		{
			res.push_back(1);
		}
		else
		{
			res.push_back(0);
		}

		x >>= 1;
		id++;
	}

	reverse(res.begin(), res.end());

	return res;
}

vector<int> sort_doubles_with_order(vector<double> &v)
{
	vector<int> order(v.size());

	iota(order.begin(), order.end(), 0);

	sort(order.begin(), order.end(),
		[&v](size_t i1, size_t i2) {return v[i1] < v[i2]; });

	return order;
}

void fill_cols_h_x(int * cols, int N, int i)
{
	int num_filled = 0;
	int chunk_size = pow(2, i);

	while (num_filled != pow(2, N))
	{
		cols[num_filled] = pow(2, i) + num_filled;
		int start_id = num_filled;

		num_filled++;

		if (chunk_size > 1)
		{
			for (int k = 1; k < chunk_size; k++)
			{
				cols[num_filled] = cols[num_filled - 1] + 1;
				num_filled++;
			}
		}

		cols[num_filled] = cols[start_id] - chunk_size;
		num_filled++;

		if (chunk_size > 1)
		{
			for (int k = 1; k < chunk_size; k++)
			{
				cols[num_filled] = cols[num_filled - 1] + 1;
				num_filled++;
			}
		}
	}
}

