#pragma once
#include "coef_coord.h"
#include "stdio.h"
#include "math.h"
#include "string.h"
#include "sortTensor.h"
#include "Matrix.h"

ulli fijk_coord_sym(crsMatrix *sel, int N);
ulli dijk_coord_sym(crsMatrix *sel, int N);
