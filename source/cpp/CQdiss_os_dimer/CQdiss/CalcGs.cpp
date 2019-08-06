#include "CalcGs.h"

void calc_Gs(Model *m, ConfigParam &cp)
{
	int N_mat = m->N_mat;
	crsMatrix* Gs = new crsMatrix();
	dcomplex sum;
	sum.re = 1.0;
	sum.im = 0.0;
	SparseMKLAdd(*(m->Qs_base), sum, *(m->Rs), *Gs);
	m->Gs = Gs;
}