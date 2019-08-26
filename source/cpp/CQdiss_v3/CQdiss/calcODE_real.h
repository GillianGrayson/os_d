#ifndef __CALC_REAL_ODE__
#define __CALC_REAL_ODE__

#include "Model.h"

void complex_to_real(dcomplex *mat, int N);
void real_to_complex(dcomplex *mat, int N);
void calcODE_real(Model *m, double h, int cntItr, double t = 0.0);

#endif