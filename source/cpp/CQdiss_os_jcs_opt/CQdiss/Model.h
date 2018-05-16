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

	crsMatrix *h_0;
	crsMatrix *h_1;

	crsMatrix * H_0;
	crsMatrix * H_1;


	crsMatrix * H0;
	crsMatrix * H1;

	crsMatrix * Q_0;
	crsMatrix * Q_1;

	dcomplex  * Ks;
	crsMatrix * Rs;
	crsMatrix * G_0_s;
	crsMatrix * G_1_s;

	dcomplex  * prevRhoF;
	dcomplex  * RhoF;

	crsMatrix * Rho;

	crsMatrix * l_mat;

	Tensor_Coordinates *f_ijk;
};

Model * createModel(int N, ConfigParam conf);
void freeModel(Model * model);

crsMatrix * create_a_std_matrix(Model * m);
crsMatrix * create_a_dag_matrix(Model * m);

#endif