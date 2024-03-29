#pragma once
#include "Config.h"
#include <mkl.h>
#include "Matrix.h"

#define Error( msg ) error( msg, __FUNCTION__, __FILE__, __LINE__ )

void error(const string& err, const char* func, const char* file, int line);

int n_choose_k(int n, int k);

int bit_count(int value);

int bit_at(int value, int position);

void print_int_array(int * data, int N);

string file_name_suffix(ConfigParam &cp, int precision);

void save_double_data(string file_name, double * data, int size, int precision, bool append);

void save_2d_double_data(string file_name, double ** data, int size_1, int size_2, int precision, bool append);

void save_int_data(string file_name, int * data, int size, bool append);

void save_complex_data(string file_name, MKL_Complex16 * data, int size, int precision, bool append);

void save_sparse_complex_mtx(string file_name, crsMatrix *A, int precision, bool append);

vector<int> convert_int_to_vector_of_bits(int x, int size);

vector<int> sort_doubles_with_order(vector<double> &v);

void fill_cols_h_x(int * cols, int N, int i);

void scalar_mult(crsMatrix * res, double gamma);