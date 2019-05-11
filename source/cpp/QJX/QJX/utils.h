#pragma once
#include "config.h"

#define Error( msg ) error( msg, __FUNCTION__, __FILE__, __LINE__ )

void error(const string& err, const char* func, const char* file, int line);

int n_choose_k(int n, int k);

int bit_count(int value);

int bit_at(int value, int position);

void print_int_array(int * data, int N);

void save_double_data(string file_name, double * data, int size, int precision, bool append);

void save_2d_double_data(string file_name, double * data, int num_rows, int num_cols, int precision, bool append);

void save_double_vector(string file_name, vector<double> vec, int precision, bool append);

void save_int_vector(string file_name, vector<int> vec, bool append);

void save_2d_inv_double_data(string file_name, double * data, int num_rows, int num_cols, int precision, bool append);

void save_2d_inv_complex_data(string file_name, MKL_Complex16 * data, int num_rows, int num_cols, int precision, bool append);

void save_int_data(string file_name, int * data, int size, bool append);

void save_complex_data(string file_name, MKL_Complex16 * data, int size, int precision, bool append);

vector<int> convert_int_to_vector_of_bits(int x, int size);

vector<int> sort_doubles_with_order(vector<double> &v);

void fill_cols_h_x(int * cols, int N, int i);