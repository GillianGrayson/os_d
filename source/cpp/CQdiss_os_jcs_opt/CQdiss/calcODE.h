#ifndef __CALC_ODE__
#define __CALC_ODE__

#include "Model.h"
#include "data.h"
#include "CalcRho.h"

void init_conditions(dcomplex *mtx, RunParam &rp, ConfigParam &cp, MainData &md);

void multMatVec(crsMatrix *mat, double * x, double * res);
void multMatVec_complex(crsMatrix *mat, dcomplex * x, dcomplex * res);

void calcVectValue_t0(double h, Model * m, double *x, double * res, double * tmp);
void calcVectValue_t1(double h, Model * m, double *x, double * res, double * tmp);

void initRhoODE(Model *m, RunParam &rp, ConfigParam &cp, MainData &md);

dcomplex calcDiffIter(Model *m);

void calcODE_trans(Model *m, RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd);
void calcODE_std(Model *m, RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd);
void calcODE_deep(Model *m, RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd);

void before(Model *m);
void after(Model *m);

void complex_to_real(dcomplex *mat, int N);
void real_to_complex(dcomplex *mat, int N);

#endif