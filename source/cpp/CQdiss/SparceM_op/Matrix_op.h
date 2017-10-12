#ifndef __MATRIX_OP__
#define __MATRIX_OP__

#include "Matrix.h"

int SparseMKLMult(crsMatrix &A, crsMatrix &B, crsMatrix &C, bool resize = false);
int SparseMKLMultOne(crsMatrix &A, crsMatrix &B, crsMatrix &C, bool resize = false);
int SparseMKLMultOneT(crsMatrix &A, crsMatrix &B, crsMatrix &C, bool resize = false);

int SparseMKLAdd (crsMatrix &A, dcomplex beta, crsMatrix &B, crsMatrix &C, bool resize = false);
int SparseMKLAddT(crsMatrix &A, dcomplex beta, crsMatrix &B, crsMatrix &C, bool resize = false);

int SparseMKLAddOne (crsMatrix &A, dcomplex beta, crsMatrix &B, crsMatrix &C, bool resize = false);
//int SparseMKLAddOneT(crsMatrix &A, dcomplex beta, crsMatrix &B, crsMatrix &C, bool resize = false);


void toOneBase(crsMatrix &A);
void toZeroBase(crsMatrix &A);

void Transpose(crsMatrix &Mat, crsMatrix &TMat, bool conj = true);

dcomplex trace(crsMatrix &A);
int trace_struct(crsMatrix &A);

void printMatrix(crsMatrix *A);
void printMatrixVal(crsMatrix *A);
void printVectorVal(dcomplex *A, int N);

void saveMatrixVal(char* file, crsMatrix *A);
void saveMatrix_coor(char* file, crsMatrix *A);

void saveAngleMatrixVal(char* file, crsMatrix *A);
void saveAbsMatrixVal(char* file, crsMatrix *A);

void printMatrixVal_com(dcomplex *A, int N);
void saveMatrixVal_com(char* file, dcomplex *A, int N);
void saveAngleMatrixVal_com(char* file, dcomplex *A, int N);
void saveAbsMatrixVal_com(char* file, dcomplex *A, int N);
void save_mtx_diag(char* file, crsMatrix *A);

void AbsMatrixDiagVal(crsMatrix *A, double * diag);
void saveVectorVal(char* file, dcomplex *vec, int N, int M);
void saveMatrix(char* file, crsMatrix *A);

void save_dense_vector(char* file, dcomplex * vec, int size);

#endif