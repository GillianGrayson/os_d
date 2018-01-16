#pragma once
#include "config.h"

struct Split
{
	bool type;
	Split * next;
	Split * prev;
	double dt;
	int steps;
	int counter;
	int N;
	MKL_Complex16 * matrix;
	double * g;
};