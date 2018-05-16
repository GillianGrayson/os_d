#ifndef __CALC_QS__
#define __CALC_QS__

#include "Model.h"

ulli countSelect(Tensor_Coordinates * f_ijk, crsMatrix *hMat, ulli from, ulli to);
void dataSelect(int N_mat, Tensor_Coordinates * f_ijk, crsMatrix *hMat, unsigned int from, unsigned int to, Tensor_Coordinates * res);

void calc_CooQs(int N_mat, Model *m, Tensor_Coordinates * f_ijk, crsMatrix *hMat, crsMatrix *&res);

void calc_Q_0(Model * m);
void calc_Q_1(Model * m);

#endif