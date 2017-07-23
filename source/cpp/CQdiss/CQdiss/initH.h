#ifndef __INIT_H__
#define __INIT_H__

#include "Model.h"

crsMatrix * createHmatrix(Model * m, bool sub_trace = true);
crsMatrix * createHe_matrix(Model * m);
void init_h_vector(Model * m);
void init_h_vector_opt(Model * m);
void init_he_vector(Model * m);
void init_he_vector_opt(Model * m);

void to_F_basis(crsMatrix * Mat, crsMatrix * vec);

void print_h_vector(Model * m);

#endif