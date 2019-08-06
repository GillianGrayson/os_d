#ifndef __MODEL__
#define __MODEL__
#include "Matrix_op.h"
#include "Config.h"

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

	dcomplex *h_0;
	dcomplex *h_1;

	crsMatrix * H_0;
	crsMatrix * H_1;

	crsMatrix * H0;
	crsMatrix * H1;

	crsMatrix ** f_mat;
	crsMatrix ** f_H_mat;
	crsMatrix ** d_mat;

	crsMatrix * a_mat;

	crsMatrix * Q_0;
	crsMatrix * Q_1;

	dcomplex  * Ks;
	crsMatrix * Rs;
	crsMatrix * G_0_s;
	crsMatrix * G_1_s;

	dcomplex  * prevRhoF;
	dcomplex  * RhoF;

	crsMatrix * Rho;
};

Model * createModel(int N, ConfigParam conf);
void freeModel(Model * model);

void createFMatrixs(FMatrixs * Fs, int N);
void freeFMatrixs(FMatrixs * Fs);

crsMatrix * create_a_std_matrix(Model * m);
crsMatrix * create_a_dag_matrix(Model * m);

#endif