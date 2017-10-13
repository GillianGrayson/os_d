#include "Matrix_op.h"
#include <stdio.h>
#include <string.h>
#include <math.h>

void saveAbsMatrixVal_com(char* file, dcomplex *A, int N)
{
	FILE *f;
	f = fopen(file, "w");
	if (f == NULL) return;

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			fprintf(f, "%.8lf ", sqrt(A[i * N + j].re * A[i * N + j].re + A[i * N + j].im * A[i * N + j].im));
		}
		fprintf(f, "\n");
	}
	fclose(f);
}

void saveMatrixVal_com(char* file, dcomplex *A, int N)
{
	FILE *f;
	char fn_re[80];
	char fn_im[80];
	sprintf(fn_re, "re_%s", file);
	sprintf(fn_im, "im_%s", file);

	f = fopen(fn_re, "w");
	if (f == NULL) return;

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			fprintf(f, "%.8lf ", A[i * N + j].re);
		}
		fprintf(f, "\n");
	}
	fclose(f);
	f = fopen(fn_im, "w");
	if (f == NULL) return;

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			fprintf(f, "%.8lf ", A[i * N + j].im);
		}
		fprintf(f, "\n");
	}
	fclose(f);
}


void printMatrixVal_com(dcomplex *A, int N)
{
	printf("***********************\n");
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			printf("%.5lf ", A[i * N + j].re);
		}
		printf("\n");
	}
	printf("\n");
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			printf("%.5lf ", A[i * N + j].im);
		}
		printf("\n");
	}

	printf("***********************\n");
}

void toOneBase(crsMatrix &A)
{
	int i, j, n = A.N;
	for (i = 0; i < A.NZ; i++)
		A.Col[i]++;
	for (j = 0; j <= n; j++)
	{
		A.RowIndex[j]++;
	}

}
void toZeroBase(crsMatrix &A)
{
	int i, j, n = A.N;
	for (i = 0; i < A.NZ; i++)
		A.Col[i]--;
	for (j = 0; j <= n; j++)
	{
		A.RowIndex[j]--;
	}
}

dcomplex trace(crsMatrix &A)
{
	dcomplex res;
	res.re = 0.0;
	res.im = 0.0;

	for (int i = 0; i < A.N; i++)
	{
		for (int k = A.RowIndex[i]; k < A.RowIndex[i + 1]; k++)
		{
			int j = A.Col[k];
			if (i == j)
			{
				res.re += A.Value[k].re;
				res.im += A.Value[k].im;
			}
		}
	}

	return res;
}

int trace_struct(crsMatrix &A)
{
	int res = 0;

	for (int i = 0; i < A.N; i++)
	{
		for (int k = A.RowIndex[i]; k < A.RowIndex[i + 1]; k++)
		{
			int j = A.Col[k];
			if (i == j)
			{
				res++;
			}
		}
	}

	return res;
}

void printMatrix(crsMatrix *A)
{
	for (int i = 0; i < A->N; i++)
	{
		int j = 0;
		for (int k = A->RowIndex[i]; k < A->RowIndex[i + 1]; k++)
		{

			while (j < A->Col[k])
			{
				printf("0 ");
				j++;
			}
			printf("1 ");
			j++;
		}
		while (j < A->N)
		{
			printf("0 ");
			j++;
		}
		printf("\n");
	}
	printf("\n");
}



void printMatrixVal(crsMatrix *A)
{
	printf("####################################\n");
	for (int i = 0; i < A->N; i++)
	{
		int j = 0;
		for (int k = A->RowIndex[i]; k < A->RowIndex[i + 1]; k++)
		{

			while (j < A->Col[k])
			{
				printf(" 0.0 ");
				j++;
			}
			printf("%3.5lf ", A->Value[k].re);
			j++;
		}
		while (j < A->N)
		{
			printf(" 0.0 ");
			j++;
		}
		printf("\n");
	}
	printf("\n");

	for (int i = 0; i < A->N; i++)
	{
		int j = 0;
		for (int k = A->RowIndex[i]; k < A->RowIndex[i + 1]; k++)
		{

			while (j < A->Col[k])
			{
				printf(" 0.0 ");
				j++;
			}
			printf("%3.5lf ", A->Value[k].im);
			j++;
		}
		while (j < A->N)
		{
			printf(" 0.0 ");
			j++;
		}
		printf("\n");
	}
	printf("####################################\n");
	printf("\n");
}

void saveMatrix_coor(char* file, crsMatrix *A)
{
	FILE *f;

	f = fopen(file, "w");
	if (f == NULL) return;

	for (int i = 0; i < A->N; i++)
	{
		for (int k = A->RowIndex[i]; k < A->RowIndex[i + 1]; k++)
		{
			int c = A->Col[k];
			fprintf(f, "%d %d 1\n", i + 1, c + 1);
		}
	}
	fclose(f);
}

void saveMatrixVal(char* file, crsMatrix *A)
{
	FILE *f;
	char fn_re[80];
	char fn_im[80];
	sprintf(fn_re, "re_%s", file);
	sprintf(fn_im, "im_%s", file);

	f = fopen(fn_re, "w");
	if (f == NULL) return;

	for (int i = 0; i < A->N; i++)
	{
		int j = 0;
		for (int k = A->RowIndex[i]; k < A->RowIndex[i + 1]; k++)
		{

			while (j < A->Col[k])
			{
				fprintf(f, "0.00000000 ");
				j++;
			}
			fprintf(f, "%.8lf ", A->Value[k].re);
			j++;
		}
		while (j < A->N)
		{
			fprintf(f, "0.00000000 ");
			j++;
		}
		fprintf(f, "\n");
	}
	fclose(f);

	f = fopen(fn_im, "w");
	if (f == NULL) return;

	for (int i = 0; i < A->N; i++)
	{
		int j = 0;
		for (int k = A->RowIndex[i]; k < A->RowIndex[i + 1]; k++)
		{

			while (j < A->Col[k])
			{
				fprintf(f, "0.00000000 ");
				j++;
			}
			fprintf(f, "%.8lf ", A->Value[k].im);
			j++;
		}
		while (j < A->N)
		{
			fprintf(f, "0.00000000 ");
			j++;
		}
		fprintf(f, "\n");
	}
	fclose(f);
}

void save_mtx_diag(char* file, crsMatrix *A)
{

	double * diag = new double[A->N];

	for (int i = 0; i < A->N; i++)
	{
		for (int k = A->RowIndex[i]; k < A->RowIndex[i + 1]; k++)
		{
			int j = A->Col[k];
			if (i == j)
			{
				diag[i] = A->Value[k].re;
			}
		}
	}

	FILE *f;

	f = fopen(file, "w");
	if (f == NULL) return;

	for (int i = 0; i < A->N; i++)
	{
		fprintf(f, "%0.16le\n", diag[i]);
	}
	fclose(f);

	delete diag;
}


void saveAbsMatrixVal(char* file, crsMatrix *A)
{
	FILE *f;
	f = fopen(file, "w");
	if (f == NULL) return;

	for (int i = 0; i < A->N; i++)
	{
		int j = 0;
		for (int k = A->RowIndex[i]; k < A->RowIndex[i + 1]; k++)
		{

			while (j < A->Col[k])
			{
				fprintf(f, "0.00000000 ");
				j++;
			}
			fprintf(f, "%.8lf ", sqrt(A->Value[k].re * A->Value[k].re + A->Value[k].im * A->Value[k].im));
			j++;
		}
		while (j < A->N)
		{
			fprintf(f, "0.00000000 ");
			j++;
		}
		fprintf(f, "\n");
	}
	fclose(f);
}

void AbsMatrixDiagVal(crsMatrix *A, double * diag)
{

	for (int i = 0; i < A->N; i++)
	{
		diag[i] = 0.0;
		for (int k = A->RowIndex[i]; k < A->RowIndex[i + 1]; k++)
		{

			if (A->Col[k] == i)
			{
				diag[i] = sqrt(A->Value[k].re * A->Value[k].re + A->Value[k].im * A->Value[k].im);
			}
		}
	}
}

void saveAngleMatrixVal(char* file, crsMatrix *A)
{
	FILE *f;
	f = fopen(file, "w");
	if (f == NULL) return;

	for (int i = 0; i < A->N; i++)
	{
		int j = 0;
		for (int k = A->RowIndex[i]; k < A->RowIndex[i + 1]; k++)
		{

			while (j < A->Col[k])
			{
				fprintf(f, "0.00000000 ");
				j++;
			}
			fprintf(f, "%.8lf ", atan2(A->Value[k].im, A->Value[k].re) * 180.0 / 3.14159265);
			j++;
		}
		while (j < A->N)
		{
			fprintf(f, "0.00000000 ");
			j++;
		}
		fprintf(f, "\n");
	}
	fclose(f);
}

void save_complex_vector(char* file, dcomplex *vec, int N)
{
	FILE *f;
	f = fopen(file, "w");
	if (f == NULL) return;

	for (int i = 0; i < N; i++)
	{
		fprintf(f, "%0.16le %0.16le\n", vec[i].re, vec[i].im);
	}
	
	fclose(f);
}

void saveVectorVal(char* file, dcomplex *vec, int N, int M)
{
	FILE *f;
	f = fopen(file, "w");
	if (f == NULL) return;
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
		{
			fprintf(f, "%0.16le ", vec[i * M + j].re);
		}
		fprintf(f, "\n");
	}
	fprintf(f, "\n");
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
		{
			fprintf(f, "%.16le ", vec[i * M + j].im);
		}
		fprintf(f, "\n");
	}

	fclose(f);
}

void printVectorVal(dcomplex *A, int N)
{
	printf("####################################\n");
	for (int i = 0; i < N; i++)
	{
		printf("%3.4lf ", A[i].re);
	}
	printf("\n");
	for (int i = 0; i < N; i++)
	{
		printf("%3.1lf ", A[i].im);
	}
	printf("\n");
	printf("####################################\n");
	printf("\n");
}

/**
* API
*   int Transpose(int n, int* column, int* row, FLOAT_TYPE* val,
int** tColumn, int** tRow, FLOAT_TYPE** tVal);
*   ���������������� �������
* INPUT
*   int  n                 - ������ �������
*   int* column           - CRS �������� �������
*   int* row
*   FLOAT_TYPE* val
* OUTPUT
*   int** tColumn           - CRS �������� �������
*   int** tRow
*   FLOAT_TYPE* tVal
* RETURN
* ������������ ��� ������
**/
void Transpose(crsMatrix &Mat, crsMatrix &TMat, bool conj)
{
	int i, j, nz;
	int S;

	int n = Mat.N;
	int* column = Mat.Col;
	int* row = Mat.RowIndex;
	dcomplex* val = Mat.Value;

	nz = row[n];

	int* tColumn = TMat.Col;
	int* tRow = TMat.RowIndex;
	dcomplex* tVal = TMat.Value;

	memset(tRow, 0, (n + 1) * sizeof(int));
	for (i = 0; i < nz; i++)
		tRow[column[i] + 1]++;

	S = 0;
	for (i = 1; i <= n; i++)
	{
		int tmp = tRow[i];
		tRow[i] = S;
		S = S + tmp;
	}

	for (i = 0; i < n; i++)
	{
		int j1 = row[i];
		int j2 = row[i + 1];
		int Col = i; // ������� � AT - ������ � �
		for (j = j1; j < j2; j++)
		{
			dcomplex V = val[j];  // ��������
			int RIndex = column[j];  // ������ � AT
			int IIndex = tRow[RIndex + 1];
			tVal[IIndex] = V;
			tColumn[IIndex] = Col;
			tRow[RIndex + 1]++;
		}
	}
	if (conj)
	{
		for (i = 0; i < nz; i++)
		{
			tVal[i].im = -tVal[i].im;
		}
	}
}


void saveMatrix(char* file, crsMatrix *A)
{
	FILE *f;
	f = fopen(file, "w");
	if (f == NULL) return;

	for (int i = 0; i < A->N; i++)
	{
		for (int k = A->RowIndex[i]; k < A->RowIndex[i + 1]; k++)
		{
			fprintf(f, "%d %d %.16le %.16le \n", i + 1, A->Col[k] + 1, A->Value[k].re, A->Value[k].im);
		}
	}
	fclose(f);
}


void save_dense_vector(char* file, dcomplex * vec, int size)
{
	FILE *f;
	f = fopen(file, "w");
	if (f == NULL) return;

	for (int i = 0; i < size; i++)
	{
		fprintf(f, "%.16le %.16le \n", vec[i].re, vec[i].im);
	}
	fclose(f);
}