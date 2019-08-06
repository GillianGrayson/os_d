#include "CalcGs.h"

void calcGs(Model *m)
{
  int N_mat = m->N_mat;
  crsMatrix * Gs = new crsMatrix();
  dcomplex sum;
  sum.re = 1.0;
  sum.im = 0.0;
  SparseMKLAdd(*(m->Qs),sum, *(m->Rs), *Gs);
  m->Gs = Gs;
  saveMatrix("Gs.txt", Gs);
}