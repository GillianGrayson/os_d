#ifndef __CALC_QS__
#define __CALC_QS__

#include "Model.h"

void calc_CooQs(int N_mat, Model *m, Tensor_Coordinates * f_ijk, crsMatrix *hMat, crsMatrix *&res);

void calc_Qs_base(Model * m);
void calc_Qs_drv(Model * m);

#endif