#pragma once
#include "config.h"
#include "utils.h"

struct MainData
{
	int sys_size;					// Number of states in system
	int num_diss;					// Number of dissipators

	double * hamiltonian;			// Hamiltonian
	double * hamiltonian_drv;		// Hamiltonian driving

	MKL_Complex16 ** dissipators;	// Dissipators
};

