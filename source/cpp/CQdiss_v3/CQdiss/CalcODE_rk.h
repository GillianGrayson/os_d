#ifndef __CALC_ODE_RK__
#define __CALC_ODE_RK__

#include "Model.h"

dcomplex * initRhoODE_rk(Model *m);
void calcODE_rk(Model *m, crsMatrix * H, crsMatrix * He, crsMatrix * L,
	            dcomplex *rho, double h, int cntItr, double t = 0.0);

#endif