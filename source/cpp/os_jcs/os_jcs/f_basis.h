#pragma once
#include "config.h"
#include "data.h"
#include "mkl_vsl.h"
#include <map>

int SparseMKLMult(crsMatrix &A, crsMatrix &B, crsMatrix &C, bool resize = false);
int SparseMKLMultOne(crsMatrix &A, crsMatrix &B, crsMatrix &C, bool resize = false);

int SparseMKLAdd(crsMatrix &A, dcomplex beta, crsMatrix &B, crsMatrix &C, bool resize = false);
int SparseMKLAddT(crsMatrix &A, dcomplex beta, crsMatrix &B, crsMatrix &C, bool resize = false);
int SparseMKLAddOne(crsMatrix &A, dcomplex beta, crsMatrix &B, crsMatrix &C, bool resize = false);

void toOneBase(crsMatrix &A);
void toZeroBase(crsMatrix &A);
dcomplex trace(crsMatrix &A);
int trace_struct(crsMatrix &A);
void printMatrix(crsMatrix *A);
void printMatrixVal(crsMatrix *A);
void saveAbsMatrixVal(char* file, crsMatrix *A);
void AbsMatrixDiagVal(crsMatrix *A, double * diag);
void saveAngleMatrixVal(char* file, crsMatrix *A);
void saveVectorVal(char* file, dcomplex *vec, int N, int M);
void printVectorVal(dcomplex *A, int N);
void Transpose(crsMatrix &Mat, crsMatrix &TMat, bool conj);
void saveMatrix(char* file, crsMatrix *A);

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

void initFs(FMatrixs *Fs, int N);
void outFs(FMatrixs *Fs);
crsMatrix * createFeyeType(int N);
crsMatrix * createFPairTypeRe(int N, int i, int j);
crsMatrix * createFPairTypeIm(int N, int i, int j);
crsMatrix * createLastType(int N, int i);

crsMatrix * create_a_std_matrix(Model * m);
crsMatrix * create_a_dag_matrix(Model * m);

crsMatrix * create_H_0_matrix(Model * m, RunParam &rp, ConfigParam &cp, MainData &md);
void init_h_0_vector(Model * m, RunParam &rp, ConfigParam &cp, MainData &md);

crsMatrix * create_H_1_matrix(Model * m, RunParam &rp, ConfigParam &cp, MainData &md);
void init_h_1_vector(Model * m, RunParam &rp, ConfigParam &cp, MainData &md);

void init_H0(Model * m, RunParam &rp, ConfigParam &cp, MainData &md);
void init_H1(Model * m, RunParam &rp, ConfigParam &cp, MainData &md);

crsMatrix * stdToCrs(vector<pair<int, dcomplex> > * mat, int N);
crsMatrix * stdToCrs(vector<map<int, dcomplex> > & mat, int N);

void scalar_mult(crsMatrix * res, double gamma);

vector<pair<int, dcomplex> > * create_a_std_matrix(crsMatrix * a1_mat, crsMatrix * a2_mat, int N_mat);

void init_diss_1(Model * m, RunParam &rp, ConfigParam &cp, MainData &md);

struct Tensor_Coordinates
{
	MKL_Complex16 * data;
	unsigned int * coord1;
	unsigned int * coord2;
	unsigned int * coord3;
	unsigned int * hash;
	unsigned int k;
};

struct Tensor_Coordinates_1
{
	MKL_Complex16 * data;
	unsigned int coord1;
	unsigned int * coord2;
	unsigned int * coord3;
	unsigned int N;
};

void sort_matrix(Tensor_Coordinates * matrix);
void fijk_coord(Tensor_Coordinates * f_ijk, int N);
void dijk_coord(Tensor_Coordinates * d_ijk, int N);
void free_matrix(Tensor_Coordinates * matrix);
void sort_matrix(Tensor_Coordinates * matrix, Tensor_Coordinates_1 * mat_res, int Nmat);

void init_f_d(Model *m);
void init_f_d_valentin(Model *m);

void transpFs(Model *m);

void calc_Q_0(Model * m);
void calc_Q_1(Model * m);

void calcKs(Model *m);
void calcRs(Model *m);
crsMatrix* calcSubRs(Model *m, int start, int finish);
void calc_G_0_s(Model *m);
void calc_G_1_s(Model *m);

void init_conditions(dcomplex *mtx, RunParam &rp, ConfigParam &cp, MainData &md);
void initRhoODE(Model *m, RunParam &rp, ConfigParam &cp, MainData &md);

dcomplex calcDiffIter(Model *m);

void calcODE_trans(Model *m, RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd);
void calcODE_std(Model *m, RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd);
void calcODE_deep(Model *m, RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd);

void calcRho(Model *m);

void characteristics_std(Model *m, RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd, int dump_id);
void characteristics_deep(Model *m, RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd, int dump_id);

void f_basis_init(Model* model, RunParam &rp, ConfigParam &cp, MainData &md);
