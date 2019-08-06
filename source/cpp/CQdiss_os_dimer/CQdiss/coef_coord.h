#ifndef __COEF_COORD__
#define __COEF_COORD__

#include "Matrix.h"

void allocMemMat(Tensor_Coordinates_1 * mat);

void swap_row(Tensor_Coordinates_1 &mat, int i, int j);

void sort_matrix(Tensor_Coordinates * matrix, Tensor_Coordinates_1 * mat_res, int Nmat);

void fijk_coord(Tensor_Coordinates * f_ijk, int N);

void dijk_coord(Tensor_Coordinates * d_ijk, int N);

void print_matrix(Tensor_Coordinates * matrix, FILE * f);

#endif