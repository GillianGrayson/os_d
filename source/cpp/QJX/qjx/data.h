#pragma once
#include "config.h"
#include "split.h"
#include "utils.h"

struct MainData
{
	int sys_size;						// Number of states in system
	int num_diss;						// Number of dissipators
	int num_ham_qj;						// Number of qj hamiltonians

	double * hamiltonian;				// Hamiltonian
	double * hamiltonian_drv;			// Hamiltonian driving

	MKL_Complex16 ** dissipators;		// Dissipators

	MKL_Complex16 ** hamiltonians_qj;	// Non-Hermitian Hamiltonians for qj

	Split * structure;					// Split structure
};

