#ifndef __CALC_ODE__
#define __CALC_ODE__

#include "Model.h"
#include "data.h"
#include "CalcRho.h"

void init_conditions(dcomplex *mtx, RunParam &rp, ConfigParam &cp, MainData &md);

void multMatVec(crsMatrix *mat, dcomplex * x, dcomplex * res);

void calcVectValue_t0(double h, Model * m, dcomplex *x, dcomplex * res, dcomplex * tmp);

void calcVectValue_t1(double h, Model * m, dcomplex *x, dcomplex * res, dcomplex * tmp);

void initRho(Model *m, RunParam &rp, ConfigParam &cp, MainData &md);

dcomplex calcDiffIter(Model *m);

void calcODE_trans(Model *m, RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd);

void calcODE_std(Model *m, RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd);

void calcODE_deep(Model *m, RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd);

#endif