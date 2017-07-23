#ifndef __MODEL__
#define __MODEL__
#include "Matrix_op.h"
#include "Config.h"

#ifndef PI
#define PI 3.1415926535897932384626433832795
#endif
extern double number_of_allocs;
extern bool ismap;
extern bool ismap_print;
extern double map_size;

//extern double number_of_free;

struct FMatrixs
{
  crsMatrix **F;
  int countF;
};

struct Model
{
  int N;
  int N_mat;
  ConfigParam conf;
  FMatrixs *Fs;
  crsMatrix *h;
  crsMatrix *he;

  FILE * memlog;

//  crsMatrix ** f_mat;
//  crsMatrix ** f_H_mat;
//  crsMatrix ** d_mat;

//  crsMatrix * a_mat;
  
  crsMatrix * Qs;
  crsMatrix * QEs;
  dcomplex  * Ks;
  crsMatrix * Rs;
  crsMatrix * Gs;

  dcomplex  * prevRhoF;
  dcomplex  * RhoF;

  crsMatrix * Rho;

  crsMatrix * l_mat;
  Tensor_Coordinates *f_ijk;
//  Tensor_Coordinates *d_ijk;
};

Model * createModel(int N, ConfigParam conf);
void freeModel(Model * model);

void createFMatrixs(FMatrixs * Fs, int N);
void freeFMatrixs(FMatrixs * Fs);

#endif