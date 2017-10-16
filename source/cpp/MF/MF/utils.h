#pragma once
#include "config.h"

#define Error( msg ) error( msg, __FUNCTION__, __FILE__, __LINE__ )

void error(const std::string& err, const char* func, const char* file, int line);

void print_int_array(int * data, int N);

string file_name_suffix(ConfigParam &cp, int precision);

void write_double_data(string file_name, double * data, int size, int precision, bool append);

void write_int_data(string file_name, int * data, int size, bool append);

void write_complex_data(string file_name, MKL_Complex16 * data, int size, int precision, bool append);

vector<int> sort_doubles_with_order(vector<double> &v);

void delete_data(double * data);
