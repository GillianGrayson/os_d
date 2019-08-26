#ifndef __INIT_A1_A2__
#define __INIT_A1_A2__

#include "Model.h"
#include "stdToCrs.h"

crsMatrix * createA1mat(int N);
crsMatrix * createA2mat(int N);

void init_a1_a2(Model * m);
void init_a1_a2_opt(Model * m);

#endif
