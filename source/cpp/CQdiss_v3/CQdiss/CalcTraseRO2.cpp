#include "CalcTraseRO2.h"
#include "Matrix_op.h"
#include <stdio.h>

void calcTraseRO2(Model *m)
{
  crsMatrix * Rho = m->Rho;
  crsMatrix * Rho2 = new crsMatrix;
  
  int i, j;
  toOneBase(*Rho);

  SparseMKLMultOne(*Rho, *Rho, *Rho2);

  toZeroBase(*Rho);
  toZeroBase(*Rho2);
  
  dcomplex tr= trace(*Rho2);
  
  delete Rho2;
  
  FILE * f = fopen("traseRho2.txt", "w");
  
  fprintf(f, "%1.6lf %1.6lf\n", tr.re, tr.im);
  
  fclose(f);
}
