#pragma once
#include "config.h"
#include "utils.h"

struct MainData
{
	int sys_size;

	double * hamiltonian;		// Hamiltonian
};

void init_main_data(RunParam &rp, ConfigParam &cp, MainData &md);
void delete_main_data(MainData &md);

void f_basis_prop(RunParam &rp, ConfigParam &cp, MainData &md);