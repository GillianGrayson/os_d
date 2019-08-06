#pragma once
#include "coef_coord.h"
#include "stdio.h"
#include "math.h"
#include "string.h"
#include "sortTensor.h"
#include "Matrix.h"

void fijk_coord_ch(Tensor_Coordinates * f_ijk, crsMatrix *sel, ulli NZ, int N);
void dijk_coord_ch(Tensor_Coordinates * d_ijk, crsMatrix *sel, ulli NZ, int N);
