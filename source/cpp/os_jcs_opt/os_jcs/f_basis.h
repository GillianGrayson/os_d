#pragma once
#include "config.h"
#include "data.h"
#include "mkl_vsl.h"
#include <map>

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

ulli max_bit(ulli * mas, unsigned int n);
inline void swap(Tensor_Coordinates  * mas, unsigned int i, unsigned int j);
void msd_sort(Tensor_Coordinates * mas, unsigned int from, unsigned int to, ulli bit, int threads_level);
void sort_matrix(Tensor_Coordinates * matrix);

Tensor_Coordinates * create_matrix(int NZ);
void free_matrix(Tensor_Coordinates * matrix);

void fijk_coord(Tensor_Coordinates * f_ijk, int N);
void dijk_coord(Tensor_Coordinates * d_ijk, int N);

ulli fijk_coord_sym(crsMatrix *sel, int N);
ulli dijk_coord_sym(crsMatrix *sel, int N);

void fijk_coord_ch(Tensor_Coordinates * f_ijk, crsMatrix *sel, ulli NZ, int N);
void dijk_coord_ch(Tensor_Coordinates * d_ijk, crsMatrix *sel, ulli NZ, int N);

void calc_CooQs_new(int N_mat, Model *m, Tensor_Coordinates * f_ijk, crsMatrix *hMat, crsMatrix *&res);
void calc_CooQs(int N_mat, Model *m, Tensor_Coordinates * f_ijk, crsMatrix *hMat, crsMatrix *&res);
ulli calcZ_ijk(Tensor_Coordinates *f_ijk, Tensor_Coordinates *d_ijk, Tensor_Coordinates *&Z_ijk);
void createIndex(Tensor_Coordinates *m_ijk, ulli N_Mat, unsigned int * index);
ulli multTmpRsSTD(ulli N_mat, unsigned int * indexF, unsigned int * indexZ,
	Tensor_Coordinates *z_ijk, Tensor_Coordinates *f_ijk, crsMatrix * l_mat, bool swapInd, vector<map<int, dcomplex> > & mat);

void to_F_basis_for_zeros(crsMatrix * Mat, crsMatrix * vec);
void to_F_basis(crsMatrix * Mat, crsMatrix * vec);

int SparseMKLMult(crsMatrix &A, crsMatrix &B, crsMatrix &C, bool resize = false);
int SparseMKLMultOne(crsMatrix &A, crsMatrix &B, crsMatrix &C, bool resize = false);

int SparseMKLAdd(crsMatrix &A, dcomplex beta, crsMatrix &B, crsMatrix &C, bool resize = false);
int SparseMKLAddT(crsMatrix &A, dcomplex beta, crsMatrix &B, crsMatrix &C, bool resize = false);
int SparseMKLAddOne(crsMatrix &A, dcomplex beta, crsMatrix &B, crsMatrix &C, bool resize = false);

void complex_to_real(dcomplex *mat, int N);
void real_to_complex(dcomplex *mat, int N);

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

Model * createModel(int N, ConfigParam conf);
void freeModel(Model * model);

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

void calc_Q_0(Model * m);
void calc_Q_1(Model * m);

void calcKs(Model *m);
void calcRs(Model *m);
void calc_G_0_s(Model *m);
void calc_G_1_s(Model *m);

void multMatVec(crsMatrix *mat, double * x, double * res);
void multMatVec_complex(crsMatrix *mat, dcomplex * x, dcomplex * res);

void calcVectValue_t0(double h, Model * m, double *x, double * res, double * tmp);
void calcVectValue_t1(double h, Model * m, double *x, double * res, double * tmp);

void init_conditions(dcomplex *mtx, RunParam &rp, ConfigParam &cp, MainData &md);
void initRhoODE(Model *m, RunParam &rp, ConfigParam &cp, MainData &md);

dcomplex calcDiffIter(Model *m);

void before(Model *m);
void after(Model *m);

void calcODE_trans(Model *m, RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd);
void calcODE_std(Model *m, RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd);
void calcODE_deep(Model *m, RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd);

void calcRho(Model *m);

void characteristics_std(Model *m, RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd, int dump_id);
void characteristics_deep(Model *m, RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd, int dump_id);

void f_basis_init(Model* model, RunParam &rp, ConfigParam &cp, MainData &md);
