#ifndef __COEF_COORD__
#define __COEF_COORD__

//#include "mkl.h"
#include "Matrix.h"


void sort_matrix(Tensor_Coordinates * matrix);
void fijk_coord(Tensor_Coordinates * f_ijk, int N);
void dijk_coord(Tensor_Coordinates * d_ijk, int N);
ulli fijk_coord_sym(crsMatrix *sel, int N);
ulli dijk_coord_sym(crsMatrix *sel, int N);
void fijk_coord_ch(Tensor_Coordinates * f_ijk, crsMatrix *sel, ulli NZ, int N);
void dijk_coord_ch(Tensor_Coordinates * f_ijk, crsMatrix *sel, ulli NZ, int N);


void sort_matrix(Tensor_Coordinates * matrix, Tensor_Coordinates_1 * mat_res, int Nmat);


#endif