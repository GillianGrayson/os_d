#ifndef __INIT_H__
#define __INIT_H__

#include "Model.h"
#include "data.h"

void to_F_basis_for_zeros(crsMatrix * Mat, crsMatrix * vec);
void to_F_basis(crsMatrix * Mat, crsMatrix * vec);

crsMatrix * create_H_0_matrix(Model * m, RunParam &rp, ConfigParam &cp, MainData &md);
void init_h_0_vector(Model * m, RunParam &rp, ConfigParam &cp, MainData &md);

crsMatrix * create_H_1_matrix(Model * m, RunParam &rp, ConfigParam &cp, MainData &md);
void init_h_1_vector(Model * m, RunParam &rp, ConfigParam &cp, MainData &md);

void init_H0(Model * m, RunParam &rp, ConfigParam &cp, MainData &md);
void init_H1(Model * m, RunParam &rp, ConfigParam &cp, MainData &md);


#endif