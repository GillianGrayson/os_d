#include "CalcGs.h"

void calc_G_0_s(Model *m, ConfigParam &cp)
{
	int N_mat = m->N_mat;
	crsMatrix * G_0_s = new crsMatrix();

	crsMatrix * subSum = new crsMatrix();

	dcomplex sum_0;
	sum_0.re = cp.drv_ampl;
	sum_0.im = 0.0;

	dcomplex sum_1;
	sum_1.re = 1.0;
	sum_1.im = 0.0;

	SparseMKLAdd(*(m->Q_0), sum_0, *(m->Q_1), *subSum);
	SparseMKLAdd(*(subSum), sum_1, *(m->Rs), *G_0_s);

	m->G_0_s = G_0_s;

	delete subSum;
}

void calc_G_1_s(Model *m, ConfigParam &cp)
{
	int N_mat = m->N_mat;
	crsMatrix * G_1_s = new crsMatrix();

	dcomplex sum;
	sum.re = 1.0;
	sum.im = 0.0;

	SparseMKLAdd(*(m->Q_0), sum, *(m->Rs), *G_1_s);

	m->G_1_s = G_1_s;
}