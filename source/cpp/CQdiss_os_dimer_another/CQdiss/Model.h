#ifndef __MODEL__
#define __MODEL__
#include "Matrix_op.h"
#include "Config.h"

#ifndef PI
#define PI 3.1415926535897932384626433832795
#endif

struct Model
{
	int N;
	int N_mat;
	ConfigParam conf;

	crsMatrix *h_base;
	crsMatrix *h_drv;

	crsMatrix * H_base;
	crsMatrix * H_drv;

	crsMatrix * Qs_base;
	crsMatrix * Qs_drv;

	dcomplex  * Ks;
	crsMatrix * Rs;
	crsMatrix * Gs;

	dcomplex  * prevRhoF;
	dcomplex  * RhoF;

	crsMatrix * Rho;

	crsMatrix * l_mat;

	Tensor_Coordinates *f_ijk;
};

Model * createModel(int N, ConfigParam conf);
void freeModel(Model * model);

#endif