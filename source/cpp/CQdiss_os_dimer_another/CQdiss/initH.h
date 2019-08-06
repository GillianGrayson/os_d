#ifndef __INIT_H__
#define __INIT_H__

#include "Model.h"
#include "data.h"

void to_F_basis_for_zeros(crsMatrix * Mat, crsMatrix * vec);
void to_F_basis(crsMatrix * Mat, crsMatrix * vec);

crsMatrix * create_H_base_matrix(Model * m, RunParam &rp, ConfigParam &cp, MainData &md);

crsMatrix * create_H_drv_matrix(Model * m, RunParam &rp, ConfigParam &cp, MainData &md);

void init_h_base_vector(Model * m, RunParam &rp, ConfigParam &cp, MainData &md);

void init_h_drv_vector(Model * m, RunParam &rp, ConfigParam &cp, MainData &md);

#endif