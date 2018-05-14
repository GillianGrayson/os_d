#include "f_basis.h"
#include <string.h>
#include <map>

#define IND(i, j, k) ((ulli)i)*((ulli)N * N - 1)*((ulli)N * N - 1) + ((ulli)j)*((ulli)N * N - 1) + ((ulli)k)

#define IndS(i, j) ((int)(((N - 1) + (N - i))/2.0 * (i) + (j - i - 1)) * 2)
#define IndJ(i, j) ((int)(((N - 1) + (N - i))/2.0 * (i) + (j - i - 1)) * 2 + 1)
#define IndD(l)    (N * (N-1) + l - 1)

int SparseMKLMult(crsMatrix &A, crsMatrix &B, crsMatrix &C, bool resize)
{
	int N = A.N;
	if (A.N != B.N)
		return 1;

	if (resize)
	{
		C.resize(N);
	}

	int n = A.N;

	// Настроим параметры для вызова функции MKL
	// Переиндексируем матрицы A и B с единицы
	int i, j;
	for (i = 0; i < A.NZ; i++)
		A.Col[i]++;
	for (i = 0; i < B.NZ; i++)
		B.Col[i]++;
	for (j = 0; j <= n; j++)
	{
		A.RowIndex[j]++;
		B.RowIndex[j]++;
	}

	// Используется функция, вычисляющая C = op(A) * B
	char trans = 'N'; // говорит о том, op(A) = A - не нужно транспонировать A

					  // Хитрый параметр, влияющий на то, как будет выделяться память
					  // request = 0: память для результирующей матрицы д.б. выделена заранее
					  // Если мы не знаем, сколько памяти необходимо для хранения результата,
					  // необходимо:
					  // 1) выделить память для массива индексов строк ic: "Кол-во строк+1" элементов;
					  // 2) вызвать функцию с параметром request = 1 - в массиве ic будет заполнен 
					  //                                                         последний элемент
					  // 3) выделить память для массивов c и jc 
					  //    (кол-во элементов = ic[Кол-во строк]-1)
					  // 4) вызвать функцию с параметром request = 2
	int request;

	// Еще один нетривиальный момент: есть возможность настроить, нужно ли 
	// упорядочивать матрицы A, B и C. У нас предполагается, что все матрицы
	// упорядочены, следовательно, выбираем вариант "No-No-Yes", который
	// соответствует любому значению, кроме целых чисел от 1 до 7 включительно
	int sort = 8;

	// Количество ненулевых элементов.
	// Используется только если request = 0
	int nzmax = -1;

	// Служебная информация
	int info;

	request = 1;
	if (!resize)
	{
		// Выделим память для индекса в матрице C
		C.RowIndex = new int[n + 1];
		// Сосчитаем количество ненулевых элементов в матрице C

		C.Value = 0;
		C.Col = 0;
	}
	mkl_zcsrmultcsr(&trans, &request, &sort, &n, &n, &n,
		(MKL_Complex16 *)A.Value, A.Col, A.RowIndex,
		(MKL_Complex16 *)B.Value, B.Col, B.RowIndex,
		(MKL_Complex16 *)C.Value, C.Col, C.RowIndex,
		&nzmax, &info);
	int nzc = C.RowIndex[n] - 1;
	if (!resize)
	{
		C.Value = new dcomplex[nzc];
		C.Col = new int[nzc];
	}
	else
	{
		C.resize(N, nzc);
	}
	// Сосчитаем C = A * B
	request = 2;
	mkl_zcsrmultcsr(&trans, &request, &sort, &n, &n, &n,
		(MKL_Complex16 *)A.Value, A.Col, A.RowIndex,
		(MKL_Complex16 *)B.Value, B.Col, B.RowIndex,
		(MKL_Complex16 *)C.Value, C.Col, C.RowIndex,
		&nzmax, &info);
	C.N = n;
	C.NZ = nzc;

	// Приведем к нормальному виду матрицы A, B и С
	for (i = 0; i < A.NZ; i++)
		A.Col[i]--;
	for (i = 0; i < B.NZ; i++)
		B.Col[i]--;
	for (i = 0; i < C.NZ; i++)
		C.Col[i]--;
	for (j = 0; j <= n; j++)
	{
		A.RowIndex[j]--;
		B.RowIndex[j]--;
		C.RowIndex[j]--;
	}

	return 0;
}
int SparseMKLMultOne(crsMatrix &A, crsMatrix &B, crsMatrix &C, bool resize)
{
	int N = A.N;
	if (A.N != B.N)
		return 1;

	if (resize)
	{
		C.resize(N);
	}

	int n = A.N;

	// Настроим параметры для вызова функции MKL
	// Переиндексируем матрицы A и B с единицы
	int i, j;

	// Используется функция, вычисляющая C = op(A) * B
	char trans = 'N'; // говорит о том, op(A) = A - не нужно транспонировать A

					  // Хитрый параметр, влияющий на то, как будет выделяться память
					  // request = 0: память для результирующей матрицы д.б. выделена заранее
					  // Если мы не знаем, сколько памяти необходимо для хранения результата,
					  // необходимо:
					  // 1) выделить память для массива индексов строк ic: "Кол-во строк+1" элементов;
					  // 2) вызвать функцию с параметром request = 1 - в массиве ic будет заполнен 
					  //                                                         последний элемент
					  // 3) выделить память для массивов c и jc 
					  //    (кол-во элементов = ic[Кол-во строк]-1)
					  // 4) вызвать функцию с параметром request = 2
	int request;

	// Еще один нетривиальный момент: есть возможность настроить, нужно ли 
	// упорядочивать матрицы A, B и C. У нас предполагается, что все матрицы
	// упорядочены, следовательно, выбираем вариант "No-No-Yes", который
	// соответствует любому значению, кроме целых чисел от 1 до 7 включительно
	int sort = 8;

	// Количество ненулевых элементов.
	// Используется только если request = 0
	int nzmax = -1;

	// Служебная информация
	int info;

	request = 1;
	if (!resize)
	{
		// Выделим память для индекса в матрице C
		C.RowIndex = new int[n + 1];
		// Сосчитаем количество ненулевых элементов в матрице C

		C.Value = 0;
		C.Col = 0;
	}
	mkl_zcsrmultcsr(&trans, &request, &sort, &n, &n, &n,
		(MKL_Complex16 *)A.Value, A.Col, A.RowIndex,
		(MKL_Complex16 *)B.Value, B.Col, B.RowIndex,
		(MKL_Complex16 *)C.Value, C.Col, C.RowIndex,
		&nzmax, &info);
	int nzc = C.RowIndex[n] - 1;
	if (!resize)
	{
		C.Value = new dcomplex[nzc];
		C.Col = new int[nzc];
	}
	else
	{
		C.resize(N, nzc);
	}
	// Сосчитаем C = A * B
	request = 2;
	mkl_zcsrmultcsr(&trans, &request, &sort, &n, &n, &n,
		(MKL_Complex16 *)A.Value, A.Col, A.RowIndex,
		(MKL_Complex16 *)B.Value, B.Col, B.RowIndex,
		(MKL_Complex16 *)C.Value, C.Col, C.RowIndex,
		&nzmax, &info);
	C.N = n;
	C.NZ = nzc;

	// Приведем к нормальному виду матрицы A, B и С

	return 0;
}

int SparseMKLAdd(crsMatrix &A, dcomplex beta, crsMatrix &B, crsMatrix &C, bool resize)
{
	int N = A.N;
	if (A.N != B.N)
		return 1;

	if (resize)
	{
		C.resize(N);
	}

	int n = A.N;

	// Настроим параметры для вызова функции MKL
	// Переиндексируем матрицы A и B с единицы
	int i, j;
	for (i = 0; i < A.NZ; i++)
		A.Col[i]++;
	for (i = 0; i < B.NZ; i++)
		B.Col[i]++;
	for (j = 0; j <= n; j++)
	{
		A.RowIndex[j]++;
		B.RowIndex[j]++;
	}

	// Используется функция, вычисляющая C = op(A) * B
	char trans;

	trans = 'N'; // говорит о том, op(A) = A - не нужно транспонировать A

				 // Хитрый параметр, влияющий на то, как будет выделяться память
				 // request = 0: память для результирующей матрицы д.б. выделена заранее
				 // Если мы не знаем, сколько памяти необходимо для хранения результата,
				 // необходимо:
				 // 1) выделить память для массива индексов строк ic: "Кол-во строк+1" элементов;
				 // 2) вызвать функцию с параметром request = 1 - в массиве ic будет заполнен 
				 //                                                         последний элемент
				 // 3) выделить память для массивов c и jc 
				 //    (кол-во элементов = ic[Кол-во строк]-1)
				 // 4) вызвать функцию с параметром request = 2
	int request;

	// Еще один нетривиальный момент: есть возможность настроить, нужно ли 
	// упорядочивать матрицы A, B и C. У нас предполагается, что все матрицы
	// упорядочены, следовательно, выбираем вариант "No-No-Yes", который
	// соответствует любому значению, кроме целых чисел от 1 до 7 включительно
	int sort = 8;

	// Количество ненулевых элементов.
	// Используется только если request = 0
	int nzmax = -1;

	// Служебная информация
	int info;
	request = 1;

	if (!resize)
	{
		// Выделим память для индекса в матрице C
		C.RowIndex = new int[n + 1];
		// Сосчитаем количество ненулевых элементов в матрице C
		C.Value = 0;
		C.Col = 0;
	}

	mkl_zcsradd(&trans, &request, &sort, &n, &n,
		(MKL_Complex16 *)A.Value, A.Col, A.RowIndex,
		(MKL_Complex16 *)(&beta),
		(MKL_Complex16 *)B.Value, B.Col, B.RowIndex,
		(MKL_Complex16 *)C.Value, C.Col, C.RowIndex,
		&nzmax, &info);
	int nzc = C.RowIndex[n] - 1;
	if (!resize)
	{
		C.Value = new dcomplex[nzc];
		C.Col = new int[nzc];
	}
	else
	{
		C.resize(N, nzc);
	}
	// Сосчитаем C = A * B
	request = 2;
	mkl_zcsradd(&trans, &request, &sort, &n, &n,
		(MKL_Complex16 *)A.Value, A.Col, A.RowIndex,
		(MKL_Complex16 *)(&beta),
		(MKL_Complex16 *)B.Value, B.Col, B.RowIndex,
		(MKL_Complex16 *)C.Value, C.Col, C.RowIndex,
		&nzmax, &info);
	C.N = n;
	C.NZ = nzc;

	// Приведем к нормальному виду матрицы A, B и С
	for (i = 0; i < A.NZ; i++)
		A.Col[i]--;
	for (i = 0; i < B.NZ; i++)
		B.Col[i]--;
	for (i = 0; i < C.NZ; i++)
		C.Col[i]--;
	for (j = 0; j <= n; j++)
	{
		A.RowIndex[j]--;
		B.RowIndex[j]--;
		C.RowIndex[j]--;
	}

	return 0;
}
int SparseMKLAddT(crsMatrix &A, dcomplex beta, crsMatrix &B, crsMatrix &C, bool resize)
{
	int N = A.N;
	if (A.N != B.N)
		return 1;

	if (resize)
	{
		C.resize(N);
	}

	int n = A.N;

	// Настроим параметры для вызова функции MKL
	// Переиндексируем матрицы A и B с единицы
	int i, j;
	for (i = 0; i < A.NZ; i++)
		A.Col[i]++;
	for (i = 0; i < B.NZ; i++)
		B.Col[i]++;
	for (j = 0; j <= n; j++)
	{
		A.RowIndex[j]++;
		B.RowIndex[j]++;
	}

	// Используется функция, вычисляющая C = op(A) * B
	char trans;

	trans = 'T'; // говорит о том, op(A) = A - не нужно транспонировать A

	int request;

	int sort = 8;

	int nzmax = -1;

	// Служебная информация
	int info;
	request = 1;

	if (!resize)
	{
		// Выделим память для индекса в матрице C
		C.RowIndex = new int[n + 1];
		// Сосчитаем количество ненулевых элементов в матрице C
		C.Value = 0;
		C.Col = 0;
	}

	mkl_zcsradd(&trans, &request, &sort, &n, &n,
		(MKL_Complex16 *)A.Value, A.Col, A.RowIndex,
		(MKL_Complex16 *)(&beta),
		(MKL_Complex16 *)B.Value, B.Col, B.RowIndex,
		(MKL_Complex16 *)C.Value, C.Col, C.RowIndex,
		&nzmax, &info);
	int nzc = C.RowIndex[n] - 1;
	if (!resize)
	{
		C.Value = new dcomplex[nzc];
		C.Col = new int[nzc];
	}
	else
	{
		C.resize(N, nzc);
	}
	// Сосчитаем C = A * B
	request = 2;
	mkl_zcsradd(&trans, &request, &sort, &n, &n,
		(MKL_Complex16 *)A.Value, A.Col, A.RowIndex,
		(MKL_Complex16 *)(&beta),
		(MKL_Complex16 *)B.Value, B.Col, B.RowIndex,
		(MKL_Complex16 *)C.Value, C.Col, C.RowIndex,
		&nzmax, &info);
	C.N = n;
	C.NZ = nzc;

	// Приведем к нормальному виду матрицы A, B и С
	for (i = 0; i < A.NZ; i++)
		A.Col[i]--;
	for (i = 0; i < B.NZ; i++)
		B.Col[i]--;
	for (i = 0; i < C.NZ; i++)
		C.Col[i]--;
	for (j = 0; j <= n; j++)
	{
		A.RowIndex[j]--;
		B.RowIndex[j]--;
		C.RowIndex[j]--;
	}

	return 0;
}
int SparseMKLAddOne(crsMatrix &A, dcomplex beta, crsMatrix &B, crsMatrix &C, bool resize)
{
	int N = A.N;
	if (A.N != B.N)
		return 1;

	if (resize)
	{
		C.resize(N);
	}

	int n = A.N;

	// Настроим параметры для вызова функции MKL
	// Переиндексируем матрицы A и B с единицы
	int i, j;

	// Используется функция, вычисляющая C = op(A) * B
	char trans = 'N'; // говорит о том, op(A) = A - не нужно транспонировать A

					  // Хитрый параметр, влияющий на то, как будет выделяться память
					  // request = 0: память для результирующей матрицы д.б. выделена заранее
					  // Если мы не знаем, сколько памяти необходимо для хранения результата,
					  // необходимо:
					  // 1) выделить память для массива индексов строк ic: "Кол-во строк+1" элементов;
					  // 2) вызвать функцию с параметром request = 1 - в массиве ic будет заполнен 
					  //                                                         последний элемент
					  // 3) выделить память для массивов c и jc 
					  //    (кол-во элементов = ic[Кол-во строк]-1)
					  // 4) вызвать функцию с параметром request = 2
	int request;

	// Еще один нетривиальный момент: есть возможность настроить, нужно ли 
	// упорядочивать матрицы A, B и C. У нас предполагается, что все матрицы
	// упорядочены, следовательно, выбираем вариант "No-No-Yes", который
	// соответствует любому значению, кроме целых чисел от 1 до 7 включительно
	int sort = 8;

	// Количество ненулевых элементов.
	// Используется только если request = 0
	int nzmax = -1;

	// Служебная информация
	int info;
	request = 1;

	if (!resize)
	{
		// Выделим память для индекса в матрице C
		C.RowIndex = new int[n + 1];
		// Сосчитаем количество ненулевых элементов в матрице C
		C.Value = 0;
		C.Col = 0;
	}

	mkl_zcsradd(&trans, &request, &sort, &n, &n,
		(MKL_Complex16 *)A.Value, A.Col, A.RowIndex,
		(MKL_Complex16 *)(&beta),
		(MKL_Complex16 *)B.Value, B.Col, B.RowIndex,
		(MKL_Complex16 *)C.Value, C.Col, C.RowIndex,
		&nzmax, &info);
	int nzc = C.RowIndex[n] - 1;
	if (!resize)
	{
		C.Value = new dcomplex[nzc];
		C.Col = new int[nzc];
	}
	else
	{
		C.resize(N, nzc);
	}
	// Сосчитаем C = A * B
	request = 2;
	mkl_zcsradd(&trans, &request, &sort, &n, &n,
		(MKL_Complex16 *)A.Value, A.Col, A.RowIndex,
		(MKL_Complex16 *)(&beta),
		(MKL_Complex16 *)B.Value, B.Col, B.RowIndex,
		(MKL_Complex16 *)C.Value, C.Col, C.RowIndex,
		&nzmax, &info);
	C.N = n;
	C.NZ = nzc;

	// Приведем к нормальному виду матрицы A, B и С

	return 0;
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
void saveVectorVal(char* file, dcomplex *vec, int N, int M)
{
	FILE *f;
	f = fopen(file, "w");
	if (f == NULL) return;
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
		{
			fprintf(f, "%.10lf ", vec[i * M + j].re);
		}
		fprintf(f, "\n");
	}
	fprintf(f, "\n");
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < M; j++)
		{
			fprintf(f, "%.10lf ", vec[i * M + j].im);
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
		int Col = i; // Столбец в AT - строка в А
		for (j = j1; j < j2; j++)
		{
			dcomplex V = val[j];  // Значение
			int RIndex = column[j];  // Строка в AT
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
			fprintf(f, "%d %d %.16lf %.16lf \n", i + 1, A->Col[k] + 1, A->Value[k].re, A->Value[k].im);
		}
	}
	fclose(f);
}

Model * createModel(int N, ConfigParam conf)
{
	int i;
	Model * model = new Model;

	model->N = N;
	model->N_mat = (N + 1) * (N + 1) - 1;
	model->conf = conf;

	model->Fs = new FMatrixs;
	createFMatrixs(model->Fs, N);

	model->h_0 = new dcomplex[model->N_mat];
	model->h_1 = new dcomplex[model->N_mat];

	model->H_0 = new crsMatrix(model->N_mat, model->N_mat);
	model->H_1 = new crsMatrix(model->N_mat, model->N_mat);
	
	model->H0 = NULL;
	model->H1 = NULL;

	model->f_mat = NULL;
	model->f_H_mat = NULL;
	model->d_mat = NULL;

	model->a_mat = NULL;

	model->Q_0 = NULL;
	model->Q_1 = NULL;

	model->Ks = new dcomplex[model->N_mat];
	model->Rs = NULL;
	model->G_0_s = NULL;
	model->G_1_s = NULL;

	model->prevRhoF = new dcomplex[model->N_mat];
	model->RhoF = new dcomplex[model->N_mat];
	memset(model->RhoF, 0, sizeof(dcomplex) * model->N_mat);
	memset(model->prevRhoF, 0, sizeof(dcomplex) * model->N_mat);
	for (i = 0; i < model->N_mat; i++)
	{
		model->RhoF[i].re = (double)rand() / (double)RAND_MAX;
		model->RhoF[i].im = 0.0;
	}

	model->Rho = NULL;

	return model;
}
void freeModel(Model * model)
{
	freeFMatrixs(model->Fs);
	delete model->Fs;

	delete[] model->h_0;
	delete[] model->h_1;
	
	if (model->H_0 != NULL)
	{
		delete model->H_0;
	}
	if (model->H_1 != NULL)
	{
		delete model->H_1;
	}

	if (model->H0 != NULL)
	{
		delete model->H0;
	}
	if (model->H1 != NULL)
	{
		delete model->H1;
	}

	if (model->f_mat != NULL)
	{
		int N = model->N;
		for (int i = 0; i < (N + 1) * (N + 1) - 1; i++)
		{
			delete model->f_mat[i];
		}
		delete[] model->f_mat;
	}

	if (model->f_H_mat != NULL)
	{
		int N = model->N;
		for (int i = 0; i < (N + 1) * (N + 1) - 1; i++)
		{
			delete model->f_H_mat[i];
		}
		delete[] model->f_H_mat;
	}

	if (model->d_mat != NULL)
	{
		int N = model->N;
		for (int i = 0; i < (N + 1) * (N + 1) - 1; i++)
		{
			delete model->d_mat[i];
		}
		delete[] model->d_mat;
	}

	if (model->a_mat != NULL)
	{
		delete model->a_mat;
	}

	if (model->Q_0 != NULL)
	{
		delete model->Q_0;
	}
	if (model->Q_1 != NULL)
	{
		delete model->Q_1;
	}

	if (model->Ks != NULL)
	{
		delete[] model->Ks;
	}
	if (model->Rs != NULL)
	{
		delete model->Rs;
	}
	if (model->G_0_s != NULL)
	{
		delete model->G_0_s;
	}
	if (model->G_1_s != NULL)
	{
		delete model->G_1_s;
	}

	if (model->RhoF != NULL)
	{
		delete[] model->RhoF;
	}

	if (model->prevRhoF != NULL)
	{
		delete[] model->prevRhoF;
	}

	if (model->Rho != NULL)
	{
		delete model->Rho;
	}
}

void createFMatrixs(FMatrixs * Fs, int N)
{
	Fs->countF = (2 + N) * N + 1;
	Fs->F = new crsMatrix *[Fs->countF];
	for (int i = 0; i < Fs->countF; i++)
	{
		Fs->F[i] = NULL;
	}
}
void freeFMatrixs(FMatrixs * Fs)
{
	for (int i = 0; i < Fs->countF; i++)
	{
		if (Fs->F[i] != NULL)
		{
			delete Fs->F[i];
		}
	}
	delete[] Fs->F;
}

void initFs(FMatrixs *Fs, int N)
{
	int i, j, k;
	Fs->F[0] = createFeyeType(N + 1);
	k = 1;
	for (i = 0; i < N + 1; i++)
	{
		for (j = i + 1; j < N + 1; j++)
		{
			Fs->F[k] = createFPairTypeRe(N + 1, i, j); k++;
			Fs->F[k] = createFPairTypeIm(N + 1, i, j); k++;
		}
	}

	for (i = 0; i < N; i++)
	{
		if (k < Fs->countF)
		{
			Fs->F[k] = createLastType(N + 1, i); k++;
		}
		else
		{
			throw("error count calc (no mem)");
		}
	}

	if (k != Fs->countF)
	{
		throw("error count calc (countF > k)");
	}

	//outFs(Fs);
}
void outFs(FMatrixs *Fs)
{
	for (int i = 0; i < Fs->countF; i++)
	{
		printMatrixVal(Fs->F[i]);
	}
}
crsMatrix * createFeyeType(int N)
{
	crsMatrix * mat;
	mat = new crsMatrix(N, N);
	for (int i = 0; i < N; i++)
	{
		mat->Col[i] = i;
		mat->RowIndex[i] = i;
		mat->Value[i].re = 1.0;
	}
	mat->RowIndex[N] = N;

	return mat;
}
crsMatrix * createFPairTypeRe(int N, int i, int j)
{
	double val = 1.0 / sqrt(2.0);
	crsMatrix * mat;
	mat = new crsMatrix(N, 2);
	mat->Value[0].re = val;
	mat->Value[1].re = val;
	mat->Col[0] = j;
	mat->Col[1] = i;

	for (int ii = 0; ii < i + 1; ii++)
	{
		mat->RowIndex[ii] = 0;
	}
	for (int ii = i + 1; ii < j + 1; ii++)
	{
		mat->RowIndex[ii] = 1;
	}
	for (int ii = j + 1; ii <= N; ii++)
	{
		mat->RowIndex[ii] = 2;
	}

	return mat;
}
crsMatrix * createFPairTypeIm(int N, int i, int j)
{
	double val = -1.0 / sqrt(2.0);
	crsMatrix * mat;
	mat = new crsMatrix(N, 2);
	mat->Value[0].im = val;
	mat->Value[1].im = -val;
	mat->Col[0] = j;
	mat->Col[1] = i;

	for (int ii = 0; ii < i + 1; ii++)
	{
		mat->RowIndex[ii] = 0;
	}
	for (int ii = i + 1; ii < j + 1; ii++)
	{
		mat->RowIndex[ii] = 1;
	}
	for (int ii = j + 1; ii <= N; ii++)
	{
		mat->RowIndex[ii] = 2;
	}

	return mat;
}
crsMatrix * createLastType(int N, int i)
{
	crsMatrix * mat;
	mat = new crsMatrix(N, i + 2);
	int ii;
	double val = 1.0 / sqrt((double)((i + 1)* (i + 2)));
	mat->RowIndex[0] = 0;
	for (ii = 0; ii <= i; ii++)
	{
		mat->Col[ii] = ii;
		mat->RowIndex[ii + 1] = mat->RowIndex[ii] + 1;
		mat->Value[ii].re = val;
	}
	mat->Col[ii] = ii;
	mat->RowIndex[ii + 1] = mat->RowIndex[ii] + 1;
	mat->Value[ii].re = -(i + 1) * val;
	ii++;

	for (; ii < N; ii++)
	{
		mat->RowIndex[ii + 1] = mat->RowIndex[ii];
	}

	return mat;
}

crsMatrix * create_a_std_matrix(Model * m)
{
	int N = m->N;

	crsMatrix * a_mtx = new crsMatrix(N + 1, N);

	for (int i = 0; i < N; i++)
	{
		a_mtx->RowIndex[i] = i;
		a_mtx->Col[i] = i+1;
		a_mtx->Value[i].re = sqrt(double(i + 1));
	}
	a_mtx->RowIndex[N] = N;
	a_mtx->RowIndex[N + 1] = N;

	return a_mtx;
}
crsMatrix * create_a_dag_matrix(Model * m)
{
	int N = m->N;

	crsMatrix * a_mtx = new crsMatrix(N + 1, N);

	a_mtx->RowIndex[0] = 0;
	for (int i = 0; i < N; i++)
	{
		a_mtx->RowIndex[i+1] = i;
		a_mtx->Col[i] = i;
		a_mtx->Value[i].re = sqrt(double(i + 1));
	}
	a_mtx->RowIndex[N + 1] = N;

	return a_mtx;
}

crsMatrix * create_H_0_matrix(Model * m, RunParam &rp, ConfigParam &cp, MainData &md)
{
	int N = m->N;

	crsMatrix * a_std = create_a_std_matrix(m);
	crsMatrix * a_std_copy = new crsMatrix(*a_std);
	crsMatrix * a_dag = create_a_dag_matrix(m);
	crsMatrix * a_dag_copy = new crsMatrix(*a_dag);

	if (rp.debug == 1)
	{
		string fn = rp.path + "a_std" + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, a_std, 16, false);

		fn = rp.path + "a_dag" + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, a_dag, 16, false);
	}

	crsMatrix * tmp_0 = new crsMatrix;
	crsMatrix * tmp_1 = new crsMatrix;
	crsMatrix * tmp_2 = new crsMatrix;

	SparseMKLMult(*a_dag, *a_dag_copy, *tmp_0);
	SparseMKLMult(*a_std, *a_std_copy, *tmp_1);
	SparseMKLMult(*tmp_0, *tmp_1, *tmp_2);

	double coeff = 0.5 * 1.0 / pow(cp.prm_alpha, 3);

	scalar_mult(tmp_2, coeff);

	crsMatrix * H_0 = new crsMatrix(*tmp_2);

	if (rp.debug == 1)
	{
		string fn = rp.path + "H_0" + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, H_0, 16, false);

		fn = rp.path + "a_std" + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, a_std, 16, false);

		fn = rp.path + "a_dag" + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, a_dag, 16, false);
	}

	delete a_std;
	delete a_std_copy;
	delete a_dag;
	delete a_dag_copy;
	
	delete tmp_0;
	delete tmp_1;
	delete tmp_2;

	return H_0;
}
void init_h_0_vector(Model * m, RunParam &rp, ConfigParam &cp, MainData &md)
{
	crsMatrix * H_0 = create_H_0_matrix(m, rp, cp, md), *res;
	m->H_0 = H_0;
	int N = m->N;
	FMatrixs *Fs = m->Fs;
	dcomplex * h_0 = m->h_0;

	int i = 0;
	for (i = 0; i < (N + 1) * (N + 1) - 1; i++)
	{
		res = new crsMatrix;
		SparseMKLMult(*H_0, *(Fs->F[i + 1]), *res);
		h_0[i] = trace(*res);

		delete res;
	}
}

crsMatrix * create_H_1_matrix(Model * m, RunParam &rp, ConfigParam &cp, MainData &md)
{
	int N = m->N;

	crsMatrix * a_std = create_a_std_matrix(m);
	crsMatrix * a_dag = create_a_dag_matrix(m);

	crsMatrix * tmp_0 = new crsMatrix();

	dcomplex beta = { -1.0, 0.0 };

	SparseMKLAdd(*a_dag, beta, *a_std, *tmp_0);
	
	int NZ = tmp_0->NZ;
	dcomplex * Value = tmp_0->Value;
	double real = 0.0;
	double imag = 0.0;
	for (int nz_id = 0; nz_id < NZ; nz_id++)
	{
		real = Value[nz_id].re;
		imag = Value[nz_id].im;

		Value[nz_id].re = -imag;
		Value[nz_id].im = real;
	}

	crsMatrix * H_1 = new crsMatrix(*tmp_0);

	if (rp.debug == 1)
	{
		string fn = rp.path + "H_1" + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, H_1, 16, false);
	}

	delete a_std;
	delete a_dag;

	delete tmp_0;

	return H_1;
}
void init_h_1_vector(Model * m, RunParam &rp, ConfigParam &cp, MainData &md)
{
	crsMatrix * H_1 = create_H_1_matrix(m, rp, cp, md), *res;
	m->H_1 = H_1;
	int N = m->N;
	FMatrixs *Fs = m->Fs;
	dcomplex * h_1 = m->h_1;

	int i = 0;
	for (i = 0; i < (N + 1) * (N + 1) - 1; i++)
	{
		res = new crsMatrix;
		SparseMKLMult(*H_1, *(Fs->F[i + 1]), *res);
		h_1[i] = trace(*res);

		delete res;
	}
}

void init_H0(Model * m, RunParam &rp, ConfigParam &cp, MainData &md)
{
	int N_mat = m->N_mat;
	crsMatrix * H0 = create_H_0_matrix(m, rp, cp, md);
	m->H0 = H0;
}
void init_H1(Model * m, RunParam &rp, ConfigParam &cp, MainData &md)
{
	int N_mat = m->N_mat;
	crsMatrix * H1 = create_H_1_matrix(m, rp, cp, md);
	m->H1 = H1;
}

crsMatrix * stdToCrs(vector<map<int, dcomplex> > & mat, int N)
{
	crsMatrix * res = NULL;

	int NZ = 0;
	int i, j, k, Nl;

	for (i = 0; i < N; i++)
	{
		NZ += mat[i].size();
	}
	res = new crsMatrix(N, NZ);
	k = 0;
	res->RowIndex[0] = 0;
	for (i = 0; i < N; i++)
	{
		Nl = mat[i].size();
		map<int, dcomplex>::iterator itr;
		for (itr = mat[i].begin(); itr != mat[i].end(); itr++)
		{
			res->Col[k] = itr->first;
			res->Value[k] = itr->second;
			k++;
		}
		res->RowIndex[i + 1] = k;
	}

	return res;
}
crsMatrix * stdToCrs(vector<pair<int, dcomplex> > * mat, int N)
{
	crsMatrix * res = NULL;

	int NZ = 0;
	int i, j, k, Nl;

	for (i = 0; i < N; i++)
	{
		NZ += mat[i].size();
	}
	res = new crsMatrix(N, NZ);
	k = 0;
	res->RowIndex[0] = 0;
	for (i = 0; i < N; i++)
	{
		Nl = mat[i].size();
		for (j = 0; j < Nl; j++)
		{
			res->Col[k] = mat[i][j].first;
			res->Value[k] = mat[i][j].second;
			k++;
		}
		res->RowIndex[i + 1] = k;
	}

	return res;
}

void scalar_mult(crsMatrix * res, double gamma)
{
	int NZ = res->NZ;
	dcomplex * Value = res->Value;

	for (int nz_id = 0; nz_id < NZ; nz_id++)
	{
		Value[nz_id].re *= gamma;
		Value[nz_id].im *= gamma;
	}
}

crsMatrix * create_A1_diss1_matrix(Model * m, int diss_id, RunParam &rp, ConfigParam &cp)
{
	int N = m->N;

	crsMatrix * a_std = create_a_std_matrix(m);

	return a_std;
}

crsMatrix * create_A1_diss1_matrix_new(Model * m, int diss_id, RunParam &rp, ConfigParam &cp)
{
	int N = m->N;

	crsMatrix * mat;

	mat = new crsMatrix((N + 1), (N + 1));

	int * RowIndex = mat->RowIndex;
	int * Col = mat->Col;
	dcomplex * Value = mat->Value;

	for (int state_id_1 = 0; state_id_1 < N + 1; state_id_1++)
	{
		Col[state_id_1] = state_id_1;
		RowIndex[state_id_1] = state_id_1;
		Value[state_id_1].re = 0.0;
		Value[state_id_1].im = 0.0;
	}
	RowIndex[N + 1] = N + 1;

	return mat;
}

crsMatrix * create_A2_diss1_matrix(Model * m, int diss_id, RunParam &rp, ConfigParam &cp)
{
	int N = m->N;

	crsMatrix * mat;

	mat = new crsMatrix((N + 1), (N + 1));

	int * RowIndex = mat->RowIndex;
	int * Col = mat->Col;
	dcomplex * Value = mat->Value;

	for (int state_id_1 = 0; state_id_1 < N + 1; state_id_1++)
	{
		Col[state_id_1] = state_id_1;
		RowIndex[state_id_1] = state_id_1;
		Value[state_id_1].re = 0.0;
	}
	RowIndex[N + 1] = N + 1;

	return mat;
}

crsMatrix * create_A2_diss1_matrix_new(Model * m, int diss_id, RunParam &rp, ConfigParam &cp)
{
	int N = m->N;

	crsMatrix * a_mtx = new crsMatrix(N + 1, N);

	for (int i = 0; i < N; i++)
	{
		a_mtx->RowIndex[i] = i;
		a_mtx->Col[i] = i + 1;
		a_mtx->Value[i].re = 0.0;
		a_mtx->Value[i].im = -sqrt(double(i + 1));
	}
	a_mtx->RowIndex[N] = N;
	a_mtx->RowIndex[N + 1] = N;

	return a_mtx;
}

vector<pair<int, dcomplex> > * create_a_std_matrix(crsMatrix * a1_mat, crsMatrix * a2_mat, int N_mat)
{
	vector<pair<int, dcomplex> > * a_std;
	a_std = new vector<pair<int, dcomplex> >[N_mat];

	for (int i = 0; i < N_mat; i++)
	{
		for (int j = 0; j < N_mat; j++)
		{
			dcomplex v1, v2, v3, v4, res1, res2, res;
			int c1, c2, c3, c4;
			c1 = a1_mat->RowIndex[i + 1] - a1_mat->RowIndex[i];
			c2 = a1_mat->RowIndex[j + 1] - a1_mat->RowIndex[j];
			c3 = a2_mat->RowIndex[i + 1] - a2_mat->RowIndex[i];
			c4 = a2_mat->RowIndex[j + 1] - a2_mat->RowIndex[j];

			if ((c1 + c3) * (c2 + c4) > 0)
			{

				v1.re = a1_mat->Value[a1_mat->RowIndex[i]].re * c1;
				v3.re = a2_mat->Value[a2_mat->RowIndex[i]].re * c3;
				v2.re = a1_mat->Value[a1_mat->RowIndex[j]].re * c2;
				v4.re = a2_mat->Value[a2_mat->RowIndex[j]].re * c4;

				v1.im = a1_mat->Value[a1_mat->RowIndex[i]].im * c1;
				v3.im = a2_mat->Value[a2_mat->RowIndex[i]].im * c3;
				v2.im = a1_mat->Value[a1_mat->RowIndex[j]].im * c2;
				v4.im = a2_mat->Value[a2_mat->RowIndex[j]].im * c4;

				res1.re = v1.re + v3.im;
				res1.im = v1.im - v3.re;

				res2.re = v2.re + v4.im;
				res2.im = v2.im - v4.re;

				res2.im = -res2.im;

				res.re = res1.re * res2.re - res1.im * res2.im;
				res.im = res1.re * res2.im + res1.im * res2.re;
				a_std[i].push_back(make_pair(j, res));
			}
		}
	}

	return a_std;
}

void init_diss_1(Model * m, RunParam &rp, ConfigParam &cp, MainData &md)
{
	int N = m->N;
	FMatrixs *Fs = m->Fs;

	int N_mat = m->N_mat;

	crsMatrix * result_a_matrix = NULL;

	double time = omp_get_wtime();
	double init_time = time;

	crsMatrix * A1 = create_A1_diss1_matrix(m, 0, rp, cp);
	crsMatrix * A2 = create_A2_diss1_matrix(m, 0, rp, cp);

	crsMatrix * a1_mat = new crsMatrix(N_mat, N_mat);
	crsMatrix * a2_mat = new crsMatrix(N_mat, N_mat);

	cout << "Dissipation" << endl;

	int k = 0;

	crsMatrix * res;
	int cnt;
	a1_mat->RowIndex[0] = 0;
	for (int i = 0; i < N_mat; i++)
	{
		res = new crsMatrix;
		SparseMKLMult(*A1, *(Fs->F[i + 1]), *res);
		cnt = trace_struct(*res);
		if (cnt > 0)
		{
			a1_mat->Value[k] = trace(*res);
			a1_mat->Col[k] = i;
			k++;
		}
		a1_mat->RowIndex[i + 1] = k;
		delete res;
	}
	a1_mat->NZ = k;

	k = 0;
	a2_mat->RowIndex[0] = 0;
	for (int i = 0; i < N_mat; i++)
	{
		res = new crsMatrix;
		SparseMKLMult(*A2, *(Fs->F[i + 1]), *res);
		cnt = trace_struct(*res);
		if (cnt > 0)
		{
			a2_mat->Value[k] = trace(*res);
			a2_mat->Col[k] = i;
			k++;
		}
		a2_mat->RowIndex[i + 1] = k;
		delete res;
	}
	a2_mat->NZ = k;

	vector<pair<int, dcomplex> > * a_std = create_a_std_matrix(a1_mat, a2_mat, N_mat);

	delete a1_mat;
	delete a2_mat;

	delete A1;
	delete A2;

	crsMatrix * a1_i_a2_mat = NULL;
	a1_i_a2_mat = stdToCrs(a_std, N_mat);

	delete[] a_std;

	scalar_mult(a1_i_a2_mat, cp.g);
	result_a_matrix = new crsMatrix(*a1_i_a2_mat);

	m->a_mat = new crsMatrix(*a1_i_a2_mat);

	delete a1_i_a2_mat;

	time = omp_get_wtime() - init_time;
	cout << "time of a_" << "0 : " << time << endl << endl;
}

void sort_matrix(Tensor_Coordinates * matrix)
{
	unsigned int key;
	unsigned int pos;
	MKL_Complex16 data_tmp;
	unsigned int int_tmp;

	for (unsigned int i = 0; i < matrix[0].k; i++)
	{
		key = matrix[0].hash[i];
		pos = i;
		for (unsigned int j = i + 1; j < matrix[0].k; j++)
		{
			if (matrix[0].hash[j] < key)
			{
				key = matrix[0].hash[j];
				pos = j;
			}
		}
		if (pos != i)
		{
			data_tmp = matrix[0].data[i];
			matrix[0].data[i] = matrix[0].data[pos];
			matrix[0].data[pos] = data_tmp;

			int_tmp = matrix[0].coord1[i];
			matrix[0].coord1[i] = matrix[0].coord1[pos];
			matrix[0].coord1[pos] = int_tmp;

			int_tmp = matrix[0].coord2[i];
			matrix[0].coord2[i] = matrix[0].coord2[pos];
			matrix[0].coord2[pos] = int_tmp;

			int_tmp = matrix[0].coord3[i];
			matrix[0].coord3[i] = matrix[0].coord3[pos];
			matrix[0].coord3[pos] = int_tmp;

			int_tmp = matrix[0].hash[i];
			matrix[0].hash[i] = matrix[0].hash[pos];
			matrix[0].hash[pos] = int_tmp;
		}
	}
}
void allocMemMat(Tensor_Coordinates_1 * mat)
{
	mat->coord2 = new unsigned  int[mat->N];
	mat->coord3 = new unsigned  int[mat->N];
	mat->data = new MKL_Complex16[mat->N];
}
void swap_row(Tensor_Coordinates_1 &mat, int i, int j)
{
	unsigned int t;
	MKL_Complex16 v;
	t = mat.coord2[i]; mat.coord2[i] = mat.coord2[j]; mat.coord2[j] = t;
	t = mat.coord3[i]; mat.coord3[i] = mat.coord3[j]; mat.coord3[j] = t;
	v = mat.data[i]; mat.data[i] = mat.data[j]; mat.data[j] = v;
}
void sort_matrix(Tensor_Coordinates * matrix, Tensor_Coordinates_1 * mat_res, int Nmat)
{
	int i = 0;
	int NN = matrix->k;
	for (i = 0; i < Nmat; i++)
	{
		mat_res[i].coord1 = i;
		mat_res[i].N = 0;
	}

	for (i = 0; i < NN; i++)
	{
		int jj = matrix->coord1[i];
		if (jj < Nmat)
		{
			mat_res[jj].N++;
		}
		else
		{
			printf("*");
		}
	}

	for (i = 0; i < Nmat; i++)
	{
		allocMemMat(mat_res + i);
		mat_res[i].N = 0;
	}

	int jj, ii, k;
	for (i = 0; i < NN; i++)
	{
		jj = matrix->coord1[i];
		ii = mat_res[jj].N;
		mat_res[jj].coord2[ii] = matrix->coord2[i];
		mat_res[jj].coord3[ii] = matrix->coord3[i];
		mat_res[jj].data[ii] = matrix->data[i];
		mat_res[jj].N++;
	}

	for (i = 0; i < Nmat; i++)
	{
		Tensor_Coordinates_1 tmp = mat_res[i];
		for (ii = 0; ii < tmp.N - 1; ii++)
		{
			for (jj = 0; jj < tmp.N - 1; jj++)
			{
				if (tmp.coord2[jj] > tmp.coord2[jj + 1])
				{
					swap_row(tmp, jj, jj + 1);
				}
				else
				{
					if (tmp.coord2[jj] == tmp.coord2[jj + 1])
					{
						if (tmp.coord3[jj] > tmp.coord3[jj + 1])
						{
							swap_row(tmp, jj, jj + 1);
						}
					}
				}
			}
		}
	}
}
void fijk_coord(Tensor_Coordinates * f_ijk, int N)
{
	unsigned int size = 5 * N * N * N - 9 * N * N - 2 * N + 6;
	f_ijk[0].data = new MKL_Complex16[size];
	memset(f_ijk[0].data, 0, size * sizeof(MKL_Complex16));
	f_ijk[0].coord1 = new unsigned int[size];
	memset(f_ijk[0].coord1, 0, size * sizeof(unsigned int));
	f_ijk[0].coord2 = new unsigned int[size];
	memset(f_ijk[0].coord2, 0, size * sizeof(unsigned int));
	f_ijk[0].coord3 = new unsigned int[size];
	memset(f_ijk[0].coord3, 0, size * sizeof(unsigned int));
	f_ijk[0].hash = new unsigned int[size];
	memset(f_ijk[0].hash, 0, size * sizeof(unsigned int));
	f_ijk[0].k = 0;

	for (int i = 1; i < N; i++)
	{
		for (int j = i + 1; j < N; j++)
		{
			f_ijk[0].data[f_ijk[0].k].real = (i) / sqrt((double)(i) * (i + 1));
			f_ijk[0].coord1[f_ijk[0].k] = IndD(i);
			f_ijk[0].coord2[f_ijk[0].k] = IndJ(i, j);
			f_ijk[0].coord3[f_ijk[0].k] = IndS(i, j);
			f_ijk[0].hash[f_ijk[0].k] = IND(IndD(i), IndJ(i, j), IndS(i, j));
			f_ijk[0].k++;

			f_ijk[0].data[f_ijk[0].k].real = -(i) / sqrt((double)(i) * (i + 1));
			f_ijk[0].coord1[f_ijk[0].k] = IndD(i);
			f_ijk[0].coord2[f_ijk[0].k] = IndS(i, j);
			f_ijk[0].coord3[f_ijk[0].k] = IndJ(i, j);
			f_ijk[0].hash[f_ijk[0].k] = IND(IndD(i), IndS(i, j), IndJ(i, j));
			f_ijk[0].k++;

			f_ijk[0].data[f_ijk[0].k].real = -(i) / sqrt((double)(i) * (i + 1));
			f_ijk[0].coord1[f_ijk[0].k] = IndJ(i, j);
			f_ijk[0].coord2[f_ijk[0].k] = IndD(i);
			f_ijk[0].coord3[f_ijk[0].k] = IndS(i, j);
			f_ijk[0].hash[f_ijk[0].k] = IND(IndJ(i, j), IndD(i), IndS(i, j));
			f_ijk[0].k++;

			f_ijk[0].data[f_ijk[0].k].real = (i) / sqrt((double)(i) * (i + 1));
			f_ijk[0].coord1[f_ijk[0].k] = IndS(i, j);
			f_ijk[0].coord2[f_ijk[0].k] = IndD(i);
			f_ijk[0].coord3[f_ijk[0].k] = IndJ(i, j);
			f_ijk[0].hash[f_ijk[0].k] = IND(IndS(i, j), IndD(i), IndJ(i, j));
			f_ijk[0].k++;

			f_ijk[0].data[f_ijk[0].k].real = (i) / sqrt((double)(i) * (i + 1));
			f_ijk[0].coord1[f_ijk[0].k] = IndJ(i, j);
			f_ijk[0].coord2[f_ijk[0].k] = IndS(i, j);
			f_ijk[0].coord3[f_ijk[0].k] = IndD(i);
			f_ijk[0].hash[f_ijk[0].k] = IND(IndJ(i, j), IndS(i, j), IndD(i));
			f_ijk[0].k++;

			f_ijk[0].data[f_ijk[0].k].real = -(i) / sqrt((double)(i) * (i + 1));
			f_ijk[0].coord1[f_ijk[0].k] = IndS(i, j);
			f_ijk[0].coord2[f_ijk[0].k] = IndJ(i, j);
			f_ijk[0].coord3[f_ijk[0].k] = IndD(i);
			f_ijk[0].hash[f_ijk[0].k] = IND(IndS(i, j), IndJ(i, j), IndD(i));
			f_ijk[0].k++;
		}
	}

	for (int i = 0; i < N; i++)
	{
		for (int j = i + 1; j < N; j++)
		{
			for (int m = i + 1; m < j; m++)
			{
				f_ijk[0].data[f_ijk[0].k].real = -1.0 / sqrt((double)(m) * (m + 1));
				f_ijk[0].coord1[f_ijk[0].k] = IndD(m);
				f_ijk[0].coord2[f_ijk[0].k] = IndJ(i, j);
				f_ijk[0].coord3[f_ijk[0].k] = IndS(i, j);
				f_ijk[0].hash[f_ijk[0].k] = IND(IndD(m), IndJ(i, j), IndS(i, j));
				f_ijk[0].k++;

				f_ijk[0].data[f_ijk[0].k].real = 1.0 / sqrt((double)(m) * (m + 1));
				f_ijk[0].coord1[f_ijk[0].k] = IndD(m);
				f_ijk[0].coord2[f_ijk[0].k] = IndS(i, j);
				f_ijk[0].coord3[f_ijk[0].k] = IndJ(i, j);
				f_ijk[0].hash[f_ijk[0].k] = IND(IndD(m), IndS(i, j), IndJ(i, j));
				f_ijk[0].k++;

				f_ijk[0].data[f_ijk[0].k].real = 1.0 / sqrt((double)(m) * (m + 1));
				f_ijk[0].coord1[f_ijk[0].k] = IndJ(i, j);
				f_ijk[0].coord2[f_ijk[0].k] = IndD(m);
				f_ijk[0].coord3[f_ijk[0].k] = IndS(i, j);
				f_ijk[0].hash[f_ijk[0].k] = IND(IndJ(i, j), IndD(m), IndS(i, j));
				f_ijk[0].k++;

				f_ijk[0].data[f_ijk[0].k].real = -1.0 / sqrt((double)(m) * (m + 1));
				f_ijk[0].coord1[f_ijk[0].k] = IndS(i, j);
				f_ijk[0].coord2[f_ijk[0].k] = IndD(m);
				f_ijk[0].coord3[f_ijk[0].k] = IndJ(i, j);
				f_ijk[0].hash[f_ijk[0].k] = IND(IndS(i, j), IndD(m), IndJ(i, j));
				f_ijk[0].k++;

				f_ijk[0].data[f_ijk[0].k].real = -1.0 / sqrt((double)(m) * (m + 1));
				f_ijk[0].coord1[f_ijk[0].k] = IndJ(i, j);
				f_ijk[0].coord2[f_ijk[0].k] = IndS(i, j);
				f_ijk[0].coord3[f_ijk[0].k] = IndD(m);
				f_ijk[0].hash[f_ijk[0].k] = IND(IndJ(i, j), IndS(i, j), IndD(m));
				f_ijk[0].k++;

				f_ijk[0].data[f_ijk[0].k].real = 1.0 / sqrt((double)(m) * (m + 1));
				f_ijk[0].coord1[f_ijk[0].k] = IndS(i, j);
				f_ijk[0].coord2[f_ijk[0].k] = IndJ(i, j);
				f_ijk[0].coord3[f_ijk[0].k] = IndD(m);
				f_ijk[0].hash[f_ijk[0].k] = IND(IndS(i, j), IndJ(i, j), IndD(m));
				f_ijk[0].k++;
			}
		}
	}

	for (int i = 0; i < N; i++)
	{
		for (int j = i + 1; j < N; j++)
		{
			f_ijk[0].data[f_ijk[0].k].real = -(1 + j) / sqrt((double)(j) * (j + 1));
			f_ijk[0].coord1[f_ijk[0].k] = IndD(j);
			f_ijk[0].coord2[f_ijk[0].k] = IndJ(i, j);
			f_ijk[0].coord3[f_ijk[0].k] = IndS(i, j);
			f_ijk[0].hash[f_ijk[0].k] = IND(IndD(j), IndJ(i, j), IndS(i, j));
			f_ijk[0].k++;

			f_ijk[0].data[f_ijk[0].k].real = (1 + j) / sqrt((double)(j) * (j + 1));
			f_ijk[0].coord1[f_ijk[0].k] = IndD(j);
			f_ijk[0].coord2[f_ijk[0].k] = IndS(i, j);
			f_ijk[0].coord3[f_ijk[0].k] = IndJ(i, j);
			f_ijk[0].hash[f_ijk[0].k] = IND(IndD(j), IndS(i, j), IndJ(i, j));
			f_ijk[0].k++;

			f_ijk[0].data[f_ijk[0].k].real = (1 + j) / sqrt((double)(j) * (j + 1));
			f_ijk[0].coord1[f_ijk[0].k] = IndJ(i, j);
			f_ijk[0].coord2[f_ijk[0].k] = IndD(j);
			f_ijk[0].coord3[f_ijk[0].k] = IndS(i, j);
			f_ijk[0].hash[f_ijk[0].k] = IND(IndJ(i, j), IndD(j), IndS(i, j));
			f_ijk[0].k++;

			f_ijk[0].data[f_ijk[0].k].real = -(1 + j) / sqrt((double)(j) * (j + 1));
			f_ijk[0].coord1[f_ijk[0].k] = IndS(i, j);
			f_ijk[0].coord2[f_ijk[0].k] = IndD(j);
			f_ijk[0].coord3[f_ijk[0].k] = IndJ(i, j);
			f_ijk[0].hash[f_ijk[0].k] = IND(IndS(i, j), IndD(j), IndJ(i, j));
			f_ijk[0].k++;

			f_ijk[0].data[f_ijk[0].k].real = -(1 + j) / sqrt((double)(j) * (j + 1));
			f_ijk[0].coord1[f_ijk[0].k] = IndJ(i, j);
			f_ijk[0].coord2[f_ijk[0].k] = IndS(i, j);
			f_ijk[0].coord3[f_ijk[0].k] = IndD(j);
			f_ijk[0].hash[f_ijk[0].k] = IND(IndJ(i, j), IndS(i, j), IndD(j));
			f_ijk[0].k++;

			f_ijk[0].data[f_ijk[0].k].real = (1 + j) / sqrt((double)(j) * (j + 1));
			f_ijk[0].coord1[f_ijk[0].k] = IndS(i, j);
			f_ijk[0].coord2[f_ijk[0].k] = IndJ(i, j);
			f_ijk[0].coord3[f_ijk[0].k] = IndD(j);
			f_ijk[0].hash[f_ijk[0].k] = IND(IndS(i, j), IndJ(i, j), IndD(j));
			f_ijk[0].k++;
		}
	}

	for (int i = 0; i < N - 1; i++)
	{
		for (int j = i + 1; j < N; j++)
		{
			for (int k = 0; k < N; k++)
			{
				if (k > j)
				{
					if (k > i)
					{
						f_ijk[0].data[f_ijk[0].k].real = 1.0 / sqrt(2.0); //-i*Jjk*[Sij,Sik]  //and symmetric Ski, Jkj
						f_ijk[0].coord1[f_ijk[0].k] = IndJ(j, k);
						f_ijk[0].coord2[f_ijk[0].k] = IndS(i, j);
						f_ijk[0].coord3[f_ijk[0].k] = IndS(i, k);
						f_ijk[0].hash[f_ijk[0].k] = IND(IndJ(j, k), IndS(i, j), IndS(i, k));
						f_ijk[0].k++;

						f_ijk[0].data[f_ijk[0].k].real = 1.0 / sqrt(2.0); //-i*Jik*[Sij,Sjk]  //and symmetric Skj, Jki
						f_ijk[0].coord1[f_ijk[0].k] = IndJ(i, k);
						f_ijk[0].coord2[f_ijk[0].k] = IndS(i, j);
						f_ijk[0].coord3[f_ijk[0].k] = IndS(j, k);
						f_ijk[0].hash[f_ijk[0].k] = IND(IndJ(i, k), IndS(i, j), IndS(j, k));
						f_ijk[0].k++;

						f_ijk[0].data[f_ijk[0].k].real = -1.0 / sqrt(2.0); //-i*Sjk*[Sij,Jik]  //and symmetric Skj, Jki
						f_ijk[0].coord1[f_ijk[0].k] = IndS(j, k);
						f_ijk[0].coord2[f_ijk[0].k] = IndS(i, j);
						f_ijk[0].coord3[f_ijk[0].k] = IndJ(i, k);
						f_ijk[0].hash[f_ijk[0].k] = IND(IndS(j, k), IndS(i, j), IndJ(i, k));
						f_ijk[0].k++;

						f_ijk[0].data[f_ijk[0].k].real = -1.0 / sqrt(2.0); //-i*Sik*[Sij,Jjk]  //and symmetric Ski, Jkj
						f_ijk[0].coord1[f_ijk[0].k] = IndS(i, k);
						f_ijk[0].coord2[f_ijk[0].k] = IndS(i, j);
						f_ijk[0].coord3[f_ijk[0].k] = IndJ(j, k);
						f_ijk[0].hash[f_ijk[0].k] = IND(IndS(i, k), IndS(i, j), IndJ(j, k));
						f_ijk[0].k++;

						f_ijk[0].data[f_ijk[0].k].real = 1.0 / sqrt(2.0); //-i*Sjk*[Jij,Sik]  //and symmetric Ski, Skj
						f_ijk[0].coord1[f_ijk[0].k] = IndS(j, k);
						f_ijk[0].coord2[f_ijk[0].k] = IndJ(i, j);
						f_ijk[0].coord3[f_ijk[0].k] = IndS(i, k);
						f_ijk[0].hash[f_ijk[0].k] = IND(IndS(j, k), IndJ(i, j), IndS(i, k));
						f_ijk[0].k++;

						f_ijk[0].data[f_ijk[0].k].real = -1.0 / sqrt(2.0); //-i*Sik*[Jij,Sjk]  //and symmetric Skj, Ski
						f_ijk[0].coord1[f_ijk[0].k] = IndS(i, k);
						f_ijk[0].coord2[f_ijk[0].k] = IndJ(i, j);
						f_ijk[0].coord3[f_ijk[0].k] = IndS(j, k);
						f_ijk[0].hash[f_ijk[0].k] = IND(IndS(i, k), IndJ(i, j), IndS(j, k));
						f_ijk[0].k++;

						f_ijk[0].data[f_ijk[0].k].real = 1.0 / sqrt(2.0); //-i*Jjk*[Jij,Jik]  //and symmetric Jkj, Jki
						f_ijk[0].coord1[f_ijk[0].k] = IndJ(j, k);
						f_ijk[0].coord2[f_ijk[0].k] = IndJ(i, j);
						f_ijk[0].coord3[f_ijk[0].k] = IndJ(i, k);
						f_ijk[0].hash[f_ijk[0].k] = IND(IndJ(j, k), IndJ(i, j), IndJ(i, k));
						f_ijk[0].k++;

						f_ijk[0].data[f_ijk[0].k].real = -1.0 / sqrt(2.0); //-i*Jik*[Jij,Jjk]  //and symmetric Jki, Jkj
						f_ijk[0].coord1[f_ijk[0].k] = IndJ(i, k);
						f_ijk[0].coord2[f_ijk[0].k] = IndJ(i, j);
						f_ijk[0].coord3[f_ijk[0].k] = IndJ(j, k);
						f_ijk[0].hash[f_ijk[0].k] = IND(IndJ(i, k), IndJ(i, j), IndJ(j, k));
						f_ijk[0].k++;
					}
				}
				if (k < j)
				{
					if (k > i)
					{
						f_ijk[0].data[f_ijk[0].k].real = -1.0 / sqrt(2.0);
						f_ijk[0].coord1[f_ijk[0].k] = IndJ(k, j);
						f_ijk[0].coord2[f_ijk[0].k] = IndS(i, j);
						f_ijk[0].coord3[f_ijk[0].k] = IndS(i, k);
						f_ijk[0].hash[f_ijk[0].k] = IND(IndJ(k, j), IndS(i, j), IndS(i, k));
						f_ijk[0].k++;

						f_ijk[0].data[f_ijk[0].k].real = 1.0 / sqrt(2.0);
						f_ijk[0].coord1[f_ijk[0].k] = IndJ(i, k);
						f_ijk[0].coord2[f_ijk[0].k] = IndS(i, j);
						f_ijk[0].coord3[f_ijk[0].k] = IndS(k, j);
						f_ijk[0].hash[f_ijk[0].k] = IND(IndJ(i, k), IndS(i, j), IndS(k, j));
						f_ijk[0].k++;

						f_ijk[0].data[f_ijk[0].k].real = -1.0 / sqrt(2.0);
						f_ijk[0].coord1[f_ijk[0].k] = IndS(k, j);
						f_ijk[0].coord2[f_ijk[0].k] = IndS(i, j);
						f_ijk[0].coord3[f_ijk[0].k] = IndJ(i, k);
						f_ijk[0].hash[f_ijk[0].k] = IND(IndS(k, j), IndS(i, j), IndJ(i, k));
						f_ijk[0].k++;

						f_ijk[0].data[f_ijk[0].k].real = 1.0 / sqrt(2.0);
						f_ijk[0].coord1[f_ijk[0].k] = IndS(i, k);
						f_ijk[0].coord2[f_ijk[0].k] = IndS(i, j);
						f_ijk[0].coord3[f_ijk[0].k] = IndJ(k, j);
						f_ijk[0].hash[f_ijk[0].k] = IND(IndS(i, k), IndS(i, j), IndJ(k, j));
						f_ijk[0].k++;

						f_ijk[0].data[f_ijk[0].k].real = 1.0 / sqrt(2.0);
						f_ijk[0].coord1[f_ijk[0].k] = IndS(k, j);
						f_ijk[0].coord2[f_ijk[0].k] = IndJ(i, j);
						f_ijk[0].coord3[f_ijk[0].k] = IndS(i, k);
						f_ijk[0].hash[f_ijk[0].k] = IND(IndS(k, j), IndJ(i, j), IndS(i, k));
						f_ijk[0].k++;

						f_ijk[0].data[f_ijk[0].k].real = -1.0 / sqrt(2.0);
						f_ijk[0].coord1[f_ijk[0].k] = IndS(i, k);
						f_ijk[0].coord2[f_ijk[0].k] = IndJ(i, j);
						f_ijk[0].coord3[f_ijk[0].k] = IndS(k, j);
						f_ijk[0].hash[f_ijk[0].k] = IND(IndS(i, k), IndJ(i, j), IndS(k, j));
						f_ijk[0].k++;

						f_ijk[0].data[f_ijk[0].k].real = -1.0 / sqrt(2.0);
						f_ijk[0].coord1[f_ijk[0].k] = IndJ(k, j);
						f_ijk[0].coord2[f_ijk[0].k] = IndJ(i, j);
						f_ijk[0].coord3[f_ijk[0].k] = IndJ(i, k);
						f_ijk[0].hash[f_ijk[0].k] = IND(IndJ(k, j), IndJ(i, j), IndJ(i, k));
						f_ijk[0].k++;

						f_ijk[0].data[f_ijk[0].k].real = 1.0 / sqrt(2.0);
						f_ijk[0].coord1[f_ijk[0].k] = IndJ(i, k);
						f_ijk[0].coord2[f_ijk[0].k] = IndJ(i, j);
						f_ijk[0].coord3[f_ijk[0].k] = IndJ(k, j);
						f_ijk[0].hash[f_ijk[0].k] = IND(IndJ(i, k), IndJ(i, j), IndJ(k, j));
						f_ijk[0].k++;
					}
					if (k < i)
					{
						f_ijk[0].data[f_ijk[0].k].real = -1.0 / sqrt(2.0);
						f_ijk[0].coord1[f_ijk[0].k] = IndJ(k, j);
						f_ijk[0].coord2[f_ijk[0].k] = IndS(i, j);
						f_ijk[0].coord3[f_ijk[0].k] = IndS(k, i);
						f_ijk[0].hash[f_ijk[0].k] = IND(IndJ(k, j), IndS(i, j), IndS(k, i));
						f_ijk[0].k++;

						f_ijk[0].data[f_ijk[0].k].real = -1.0 / sqrt(2.0);
						f_ijk[0].coord1[f_ijk[0].k] = IndJ(k, i);
						f_ijk[0].coord2[f_ijk[0].k] = IndS(i, j);
						f_ijk[0].coord3[f_ijk[0].k] = IndS(k, j);
						f_ijk[0].hash[f_ijk[0].k] = IND(IndJ(k, i), IndS(i, j), IndS(k, j));
						f_ijk[0].k++;

						f_ijk[0].data[f_ijk[0].k].real = 1.0 / sqrt(2.0);
						f_ijk[0].coord1[f_ijk[0].k] = IndS(k, j);
						f_ijk[0].coord2[f_ijk[0].k] = IndS(i, j);
						f_ijk[0].coord3[f_ijk[0].k] = IndJ(k, i);
						f_ijk[0].hash[f_ijk[0].k] = IND(IndS(k, j), IndS(i, j), IndJ(k, i));
						f_ijk[0].k++;

						f_ijk[0].data[f_ijk[0].k].real = 1.0 / sqrt(2.0);
						f_ijk[0].coord1[f_ijk[0].k] = IndS(k, i);
						f_ijk[0].coord2[f_ijk[0].k] = IndS(i, j);
						f_ijk[0].coord3[f_ijk[0].k] = IndJ(k, j);
						f_ijk[0].hash[f_ijk[0].k] = IND(IndS(k, i), IndS(i, j), IndJ(k, j));
						f_ijk[0].k++;

						f_ijk[0].data[f_ijk[0].k].real = 1.0 / sqrt(2.0);
						f_ijk[0].coord1[f_ijk[0].k] = IndS(k, j);
						f_ijk[0].coord2[f_ijk[0].k] = IndJ(i, j);
						f_ijk[0].coord3[f_ijk[0].k] = IndS(k, i);
						f_ijk[0].hash[f_ijk[0].k] = IND(IndS(k, j), IndJ(i, j), IndS(k, i));
						f_ijk[0].k++;

						f_ijk[0].data[f_ijk[0].k].real = -1.0 / sqrt(2.0);
						f_ijk[0].coord1[f_ijk[0].k] = IndS(k, i);
						f_ijk[0].coord2[f_ijk[0].k] = IndJ(i, j);
						f_ijk[0].coord3[f_ijk[0].k] = IndS(k, j);
						f_ijk[0].hash[f_ijk[0].k] = IND(IndS(k, i), IndJ(i, j), IndS(k, j));
						f_ijk[0].k++;

						f_ijk[0].data[f_ijk[0].k].real = 1.0 / sqrt(2.0);
						f_ijk[0].coord1[f_ijk[0].k] = IndJ(k, j);
						f_ijk[0].coord2[f_ijk[0].k] = IndJ(i, j);
						f_ijk[0].coord3[f_ijk[0].k] = IndJ(k, i);
						f_ijk[0].hash[f_ijk[0].k] = IND(IndJ(k, j), IndJ(i, j), IndJ(k, i));
						f_ijk[0].k++;

						f_ijk[0].data[f_ijk[0].k].real = -1.0 / sqrt(2.0);
						f_ijk[0].coord1[f_ijk[0].k] = IndJ(k, i);
						f_ijk[0].coord2[f_ijk[0].k] = IndJ(i, j);
						f_ijk[0].coord3[f_ijk[0].k] = IndJ(k, j);
						f_ijk[0].hash[f_ijk[0].k] = IND(IndJ(k, i), IndJ(i, j), IndJ(k, j));
						f_ijk[0].k++;
					}
				}
			}
		}
	}
}
void dijk_coord(Tensor_Coordinates * d_ijk, int N)
{
	unsigned int size = 6 * N * N * N - (N * (21 * N + 7)) / 2 + 1;
	d_ijk[0].data = new MKL_Complex16[size];
	memset(d_ijk[0].data, 0, size * sizeof(MKL_Complex16));
	d_ijk[0].coord1 = new unsigned int[size];
	memset(d_ijk[0].coord1, 0, size * sizeof(unsigned int));
	d_ijk[0].coord2 = new unsigned int[size];
	memset(d_ijk[0].coord2, 0, size * sizeof(unsigned int));
	d_ijk[0].coord3 = new unsigned int[size];
	memset(d_ijk[0].coord3, 0, size * sizeof(unsigned int));
	d_ijk[0].hash = new unsigned int[size];
	memset(d_ijk[0].hash, 0, size * sizeof(unsigned int));
	d_ijk[0].k = 0;

	for (int i = 1; i < N; i++)
	{
		for (int j = i + 1; j < N; j++)
		{
			d_ijk[0].data[d_ijk[0].k].real = -(i) / sqrt((double)(i) * (i + 1));
			d_ijk[0].coord1[d_ijk[0].k] = IndD(i);
			d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
			d_ijk[0].coord3[d_ijk[0].k] = IndS(i, j);
			d_ijk[0].hash[d_ijk[0].k] = IND(IndD(i), IndS(i, j), IndS(i, j));
			d_ijk[0].k++;

			d_ijk[0].data[d_ijk[0].k].real = -(i) / sqrt((double)(i) * (i + 1));
			d_ijk[0].coord1[d_ijk[0].k] = IndS(i, j);
			d_ijk[0].coord2[d_ijk[0].k] = IndD(i);
			d_ijk[0].coord3[d_ijk[0].k] = IndS(i, j);
			d_ijk[0].hash[d_ijk[0].k] = IND(IndS(i, j), IndD(i), IndS(i, j));
			d_ijk[0].k++;

			d_ijk[0].data[d_ijk[0].k].real = -(i) / sqrt((double)(i) * (i + 1));
			d_ijk[0].coord1[d_ijk[0].k] = IndD(i);
			d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
			d_ijk[0].coord3[d_ijk[0].k] = IndJ(i, j);
			d_ijk[0].hash[d_ijk[0].k] = IND(IndD(i), IndJ(i, j), IndJ(i, j));
			d_ijk[0].k++;

			d_ijk[0].data[d_ijk[0].k].real = -(i) / sqrt((double)(i) * (i + 1));
			d_ijk[0].coord1[d_ijk[0].k] = IndJ(i, j);
			d_ijk[0].coord2[d_ijk[0].k] = IndD(i);
			d_ijk[0].coord3[d_ijk[0].k] = IndJ(i, j);
			d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(i, j), IndD(i), IndJ(i, j));
			d_ijk[0].k++;

			d_ijk[0].data[d_ijk[0].k].real = -(i) / sqrt((double)(i) * (i + 1));
			d_ijk[0].coord1[d_ijk[0].k] = IndS(i, j);
			d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
			d_ijk[0].coord3[d_ijk[0].k] = IndD(i);
			d_ijk[0].hash[d_ijk[0].k] = IND(IndS(i, j), IndS(i, j), IndD(i));
			d_ijk[0].k++;

			d_ijk[0].data[d_ijk[0].k].real = -(i) / sqrt((double)(i) * (i + 1));
			d_ijk[0].coord1[d_ijk[0].k] = IndJ(i, j);
			d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
			d_ijk[0].coord3[d_ijk[0].k] = IndD(i);
			d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(i, j), IndJ(i, j), IndD(i));
			d_ijk[0].k++;
		}
	}

	for (int i = 0; i < N; i++)
	{
		for (int j = i + 1; j < N; j++)
		{
			for (int m = i + 1; m < j; m++)
			{
				d_ijk[0].data[d_ijk[0].k].real = 1.0 / sqrt((double)(m) * (m + 1));
				d_ijk[0].coord1[d_ijk[0].k] = IndD(m);
				d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
				d_ijk[0].coord3[d_ijk[0].k] = IndS(i, j);
				d_ijk[0].hash[d_ijk[0].k] = IND(IndD(m), IndS(i, j), IndS(i, j));
				d_ijk[0].k++;

				d_ijk[0].data[d_ijk[0].k].real = 1.0 / sqrt((double)(m) * (m + 1));
				d_ijk[0].coord1[d_ijk[0].k] = IndS(i, j);
				d_ijk[0].coord2[d_ijk[0].k] = IndD(m);
				d_ijk[0].coord3[d_ijk[0].k] = IndS(i, j);
				d_ijk[0].hash[d_ijk[0].k] = IND(IndS(i, j), IndD(m), IndS(i, j));
				d_ijk[0].k++;

				d_ijk[0].data[d_ijk[0].k].real = 1.0 / sqrt((double)(m) * (m + 1));
				d_ijk[0].coord1[d_ijk[0].k] = IndD(m);
				d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
				d_ijk[0].coord3[d_ijk[0].k] = IndJ(i, j);
				d_ijk[0].hash[d_ijk[0].k] = IND(IndD(m), IndJ(i, j), IndJ(i, j));
				d_ijk[0].k++;

				d_ijk[0].data[d_ijk[0].k].real = 1.0 / sqrt((double)(m) * (m + 1));
				d_ijk[0].coord1[d_ijk[0].k] = IndJ(i, j);
				d_ijk[0].coord2[d_ijk[0].k] = IndD(m);
				d_ijk[0].coord3[d_ijk[0].k] = IndJ(i, j);
				d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(i, j), IndD(m), IndJ(i, j));
				d_ijk[0].k++;

				d_ijk[0].data[d_ijk[0].k].real = 1.0 / sqrt((double)(m) * (m + 1));
				d_ijk[0].coord1[d_ijk[0].k] = IndS(i, j);
				d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
				d_ijk[0].coord3[d_ijk[0].k] = IndD(m);
				d_ijk[0].hash[d_ijk[0].k] = IND(IndS(i, j), IndS(i, j), IndD(m));
				d_ijk[0].k++;

				d_ijk[0].data[d_ijk[0].k].real = 1.0 / sqrt((double)(m) * (m + 1));
				d_ijk[0].coord1[d_ijk[0].k] = IndJ(i, j);
				d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
				d_ijk[0].coord3[d_ijk[0].k] = IndD(m);
				d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(i, j), IndJ(i, j), IndD(m));
				d_ijk[0].k++;
			}
		}
	}

	for (int j = 2; j < N; j++)
	{
		int i = 0;
		d_ijk[0].data[d_ijk[0].k].real = (1 - j) / sqrt((double)(j) * (j + 1));
		d_ijk[0].coord1[d_ijk[0].k] = IndD(j);
		d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
		d_ijk[0].coord3[d_ijk[0].k] = IndS(i, j);
		d_ijk[0].hash[d_ijk[0].k] = IND(IndD(j), IndS(i, j), IndS(i, j));
		d_ijk[0].k++;

		d_ijk[0].data[d_ijk[0].k].real = (1 - j) / sqrt((double)(j) * (j + 1));
		d_ijk[0].coord1[d_ijk[0].k] = IndS(i, j);
		d_ijk[0].coord2[d_ijk[0].k] = IndD(j);
		d_ijk[0].coord3[d_ijk[0].k] = IndS(i, j);
		d_ijk[0].hash[d_ijk[0].k] = IND(IndS(i, j), IndD(j), IndS(i, j));
		d_ijk[0].k++;

		d_ijk[0].data[d_ijk[0].k].real = (1 - j) / sqrt((double)(j) * (j + 1));
		d_ijk[0].coord1[d_ijk[0].k] = IndD(j);
		d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
		d_ijk[0].coord3[d_ijk[0].k] = IndJ(i, j);
		d_ijk[0].hash[d_ijk[0].k] = IND(IndD(j), IndJ(i, j), IndJ(i, j));
		d_ijk[0].k++;

		d_ijk[0].data[d_ijk[0].k].real = (1 - j) / sqrt((double)(j) * (j + 1));
		d_ijk[0].coord1[d_ijk[0].k] = IndJ(i, j);
		d_ijk[0].coord2[d_ijk[0].k] = IndD(j);
		d_ijk[0].coord3[d_ijk[0].k] = IndJ(i, j);
		d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(i, j), IndD(j), IndJ(i, j));
		d_ijk[0].k++;

		d_ijk[0].data[d_ijk[0].k].real = (1 - j) / sqrt((double)(j) * (j + 1));
		d_ijk[0].coord1[d_ijk[0].k] = IndS(i, j);
		d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
		d_ijk[0].coord3[d_ijk[0].k] = IndD(j);
		d_ijk[0].hash[d_ijk[0].k] = IND(IndS(i, j), IndS(i, j), IndD(j));
		d_ijk[0].k++;

		d_ijk[0].data[d_ijk[0].k].real = (1 - j) / sqrt((double)(j) * (j + 1));
		d_ijk[0].coord1[d_ijk[0].k] = IndJ(i, j);
		d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
		d_ijk[0].coord3[d_ijk[0].k] = IndD(j);
		d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(i, j), IndJ(i, j), IndD(j));
		d_ijk[0].k++;
	}

	for (int i = 1; i < N; i++)
	{
		for (int j = i + 1; j < N; j++)
		{
			d_ijk[0].data[d_ijk[0].k].real = (1 - j) / sqrt((double)(j) * (j + 1));
			d_ijk[0].coord1[d_ijk[0].k] = IndD(j);
			d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
			d_ijk[0].coord3[d_ijk[0].k] = IndS(i, j);
			d_ijk[0].hash[d_ijk[0].k] = IND(IndD(j), IndS(i, j), IndS(i, j));
			d_ijk[0].k++;

			d_ijk[0].data[d_ijk[0].k].real = (1 - j) / sqrt((double)(j) * (j + 1));
			d_ijk[0].coord1[d_ijk[0].k] = IndS(i, j);
			d_ijk[0].coord2[d_ijk[0].k] = IndD(j);
			d_ijk[0].coord3[d_ijk[0].k] = IndS(i, j);
			d_ijk[0].hash[d_ijk[0].k] = IND(IndS(i, j), IndD(j), IndS(i, j));
			d_ijk[0].k++;

			d_ijk[0].data[d_ijk[0].k].real = (1 - j) / sqrt((double)(j) * (j + 1));
			d_ijk[0].coord1[d_ijk[0].k] = IndD(j);
			d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
			d_ijk[0].coord3[d_ijk[0].k] = IndJ(i, j);
			d_ijk[0].hash[d_ijk[0].k] = IND(IndD(j), IndJ(i, j), IndJ(i, j));
			d_ijk[0].k++;

			d_ijk[0].data[d_ijk[0].k].real = (1 - j) / sqrt((double)(j) * (j + 1));
			d_ijk[0].coord1[d_ijk[0].k] = IndJ(i, j);
			d_ijk[0].coord2[d_ijk[0].k] = IndD(j);
			d_ijk[0].coord3[d_ijk[0].k] = IndJ(i, j);
			d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(i, j), IndD(j), IndJ(i, j));
			d_ijk[0].k++;

			d_ijk[0].data[d_ijk[0].k].real = (1 - j) / sqrt((double)(j) * (j + 1));
			d_ijk[0].coord1[d_ijk[0].k] = IndS(i, j);
			d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
			d_ijk[0].coord3[d_ijk[0].k] = IndD(j);
			d_ijk[0].hash[d_ijk[0].k] = IND(IndS(i, j), IndS(i, j), IndD(j));
			d_ijk[0].k++;

			d_ijk[0].data[d_ijk[0].k].real = (1 - j) / sqrt((double)(j) * (j + 1));
			d_ijk[0].coord1[d_ijk[0].k] = IndJ(i, j);
			d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
			d_ijk[0].coord3[d_ijk[0].k] = IndD(j);
			d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(i, j), IndJ(i, j), IndD(j));
			d_ijk[0].k++;
		}
	}

	for (int i = 0; i < N; i++)
	{
		for (int j = i + 1; j < N; j++)
		{
			for (int z = j + 1; z < N; z++)
			{
				d_ijk[0].data[d_ijk[0].k].real = 2.0 / sqrt((double)(z) * (z + 1));
				d_ijk[0].coord1[d_ijk[0].k] = IndD(z);
				d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
				d_ijk[0].coord3[d_ijk[0].k] = IndS(i, j);
				d_ijk[0].hash[d_ijk[0].k] = IND(IndD(z), IndS(i, j), IndS(i, j));
				d_ijk[0].k++;

				d_ijk[0].data[d_ijk[0].k].real = 2.0 / sqrt((double)(z) * (z + 1));
				d_ijk[0].coord1[d_ijk[0].k] = IndS(i, j);
				d_ijk[0].coord2[d_ijk[0].k] = IndD(z);
				d_ijk[0].coord3[d_ijk[0].k] = IndS(i, j);
				d_ijk[0].hash[d_ijk[0].k] = IND(IndS(i, j), IndD(z), IndS(i, j));
				d_ijk[0].k++;

				d_ijk[0].data[d_ijk[0].k].real = 2.0 / sqrt((double)(z) * (z + 1));
				d_ijk[0].coord1[d_ijk[0].k] = IndD(z);
				d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
				d_ijk[0].coord3[d_ijk[0].k] = IndJ(i, j);
				d_ijk[0].hash[d_ijk[0].k] = IND(IndD(z), IndJ(i, j), IndJ(i, j));
				d_ijk[0].k++;

				d_ijk[0].data[d_ijk[0].k].real = 2.0 / sqrt((double)(z) * (z + 1));
				d_ijk[0].coord1[d_ijk[0].k] = IndJ(i, j);
				d_ijk[0].coord2[d_ijk[0].k] = IndD(z);
				d_ijk[0].coord3[d_ijk[0].k] = IndJ(i, j);
				d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(i, j), IndD(z), IndJ(i, j));
				d_ijk[0].k++;

				d_ijk[0].data[d_ijk[0].k].real = 2.0 / sqrt((double)(z) * (z + 1));
				d_ijk[0].coord1[d_ijk[0].k] = IndS(i, j);
				d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
				d_ijk[0].coord3[d_ijk[0].k] = IndD(z);
				d_ijk[0].hash[d_ijk[0].k] = IND(IndS(i, j), IndS(i, j), IndD(z));
				d_ijk[0].k++;

				d_ijk[0].data[d_ijk[0].k].real = 2.0 / sqrt((double)(z) * (z + 1));
				d_ijk[0].coord1[d_ijk[0].k] = IndJ(i, j);
				d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
				d_ijk[0].coord3[d_ijk[0].k] = IndD(z);
				d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(i, j), IndJ(i, j), IndD(z));
				d_ijk[0].k++;
			}
		}
	}

	for (int i = 0; i < N - 1; i++)
	{
		for (int j = i + 1; j < N; j++)
		{
			for (int k = 0; k < N; k++)
			{
				if (k > j)
				{
					if (k > i)
					{
						d_ijk[0].data[d_ijk[0].k].real = 1.0 / sqrt(2.0);  //Sjk*{Sij,Sik}  //and symmetric Ski, Skj
						d_ijk[0].coord1[d_ijk[0].k] = IndS(j, k);
						d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
						d_ijk[0].coord3[d_ijk[0].k] = IndS(i, k);
						d_ijk[0].hash[d_ijk[0].k] = IND(IndS(j, k), IndS(i, j), IndS(i, k));
						d_ijk[0].k++;

						d_ijk[0].data[d_ijk[0].k].real = 1.0 / sqrt(2.0);  //Sik*{Sij,Sjk}  //and symmetric Skj, Ski
						d_ijk[0].coord1[d_ijk[0].k] = IndS(i, k);
						d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
						d_ijk[0].coord3[d_ijk[0].k] = IndS(j, k);
						d_ijk[0].hash[d_ijk[0].k] = IND(IndS(i, k), IndS(i, j), IndS(j, k));
						d_ijk[0].k++;

						d_ijk[0].data[d_ijk[0].k].real = 1.0 / sqrt(2.0);  //Jjk*{Sij,Jik}  //and symmetric Jkj, Jki
						d_ijk[0].coord1[d_ijk[0].k] = IndJ(j, k);
						d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
						d_ijk[0].coord3[d_ijk[0].k] = IndJ(i, k);
						d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(j, k), IndS(i, j), IndJ(i, k));
						d_ijk[0].k++;

						d_ijk[0].data[d_ijk[0].k].real = 1.0 / sqrt(2.0);  //Jik*{Sij,Jjk}  //and symmetric Jki, Jkj
						d_ijk[0].coord1[d_ijk[0].k] = IndJ(i, k);
						d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
						d_ijk[0].coord3[d_ijk[0].k] = IndJ(j, k);
						d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(i, k), IndS(i, j), IndJ(j, k));
						d_ijk[0].k++;

						d_ijk[0].data[d_ijk[0].k].real = -1.0 / sqrt(2.0);  //Jjk*{Jij,Sik}  //and symmetric Ski, Jkj
						d_ijk[0].coord1[d_ijk[0].k] = IndJ(j, k);
						d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
						d_ijk[0].coord3[d_ijk[0].k] = IndS(i, k);
						d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(j, k), IndJ(i, j), IndS(i, k));
						d_ijk[0].k++;

						d_ijk[0].data[d_ijk[0].k].real = 1.0 / sqrt(2.0);  //Jik*{Jij,Sjk}  //and symmetric Skj, Jki
						d_ijk[0].coord1[d_ijk[0].k] = IndJ(i, k);
						d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
						d_ijk[0].coord3[d_ijk[0].k] = IndS(j, k);
						d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(i, k), IndJ(i, j), IndS(j, k));
						d_ijk[0].k++;

						d_ijk[0].data[d_ijk[0].k].real = 1.0 / sqrt(2.0);  //Sjk*{Jij,Jik}  //and symmetric Skj, Jki
						d_ijk[0].coord1[d_ijk[0].k] = IndS(j, k);
						d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
						d_ijk[0].coord3[d_ijk[0].k] = IndJ(i, k);
						d_ijk[0].hash[d_ijk[0].k] = IND(IndS(j, k), IndJ(i, j), IndJ(i, k));
						d_ijk[0].k++;

						d_ijk[0].data[d_ijk[0].k].real = -1.0 / sqrt(2.0);  //Sik*{Jij,Jjk}  //and symmetric Ski, Jkj
						d_ijk[0].coord1[d_ijk[0].k] = IndS(i, k);
						d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
						d_ijk[0].coord3[d_ijk[0].k] = IndJ(j, k);
						d_ijk[0].hash[d_ijk[0].k] = IND(IndS(i, k), IndJ(i, j), IndJ(j, k));
						d_ijk[0].k++;
					}
				}
				if (k < j)
				{
					if (k > i)
					{
						d_ijk[0].data[d_ijk[0].k].real = 1.0 / sqrt(2.0);
						d_ijk[0].coord1[d_ijk[0].k] = IndS(k, j);
						d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
						d_ijk[0].coord3[d_ijk[0].k] = IndS(i, k);
						d_ijk[0].hash[d_ijk[0].k] = IND(IndS(k, j), IndS(i, j), IndS(i, k));
						d_ijk[0].k++;

						d_ijk[0].data[d_ijk[0].k].real = 1.0 / sqrt(2.0);
						d_ijk[0].coord1[d_ijk[0].k] = IndS(i, k);
						d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
						d_ijk[0].coord3[d_ijk[0].k] = IndS(k, j);
						d_ijk[0].hash[d_ijk[0].k] = IND(IndS(i, k), IndS(i, j), IndS(k, j));
						d_ijk[0].k++;

						d_ijk[0].data[d_ijk[0].k].real = -1.0 / sqrt(2.0);
						d_ijk[0].coord1[d_ijk[0].k] = IndJ(k, j);
						d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
						d_ijk[0].coord3[d_ijk[0].k] = IndJ(i, k);
						d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(k, j), IndS(i, j), IndJ(i, k));
						d_ijk[0].k++;

						d_ijk[0].data[d_ijk[0].k].real = -1.0 / sqrt(2.0);
						d_ijk[0].coord1[d_ijk[0].k] = IndJ(i, k);
						d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
						d_ijk[0].coord3[d_ijk[0].k] = IndJ(k, j);
						d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(i, k), IndS(i, j), IndJ(k, j));
						d_ijk[0].k++;

						d_ijk[0].data[d_ijk[0].k].real = 1.0 / sqrt(2.0);
						d_ijk[0].coord1[d_ijk[0].k] = IndJ(k, j);
						d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
						d_ijk[0].coord3[d_ijk[0].k] = IndS(i, k);
						d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(k, j), IndJ(i, j), IndS(i, k));
						d_ijk[0].k++;

						d_ijk[0].data[d_ijk[0].k].real = 1.0 / sqrt(2.0);
						d_ijk[0].coord1[d_ijk[0].k] = IndJ(i, k);
						d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
						d_ijk[0].coord3[d_ijk[0].k] = IndS(k, j);
						d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(i, k), IndJ(i, j), IndS(k, j));
						d_ijk[0].k++;

						d_ijk[0].data[d_ijk[0].k].real = 1.0 / sqrt(2.0);
						d_ijk[0].coord1[d_ijk[0].k] = IndS(k, j);
						d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
						d_ijk[0].coord3[d_ijk[0].k] = IndJ(i, k);
						d_ijk[0].hash[d_ijk[0].k] = IND(IndS(k, j), IndJ(i, j), IndJ(i, k));
						d_ijk[0].k++;

						d_ijk[0].data[d_ijk[0].k].real = 1.0 / sqrt(2.0);
						d_ijk[0].coord1[d_ijk[0].k] = IndS(i, k);
						d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
						d_ijk[0].coord3[d_ijk[0].k] = IndJ(k, j);
						d_ijk[0].hash[d_ijk[0].k] = IND(IndS(i, k), IndJ(i, j), IndJ(k, j));
						d_ijk[0].k++;
					}
					if (k < i)
					{
						d_ijk[0].data[d_ijk[0].k].real = 1.0 / sqrt(2.0);
						d_ijk[0].coord1[d_ijk[0].k] = IndS(k, j);
						d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
						d_ijk[0].coord3[d_ijk[0].k] = IndS(k, i);
						d_ijk[0].hash[d_ijk[0].k] = IND(IndS(k, j), IndS(i, j), IndS(k, i));
						d_ijk[0].k++;

						d_ijk[0].data[d_ijk[0].k].real = 1.0 / sqrt(2.0);
						d_ijk[0].coord1[d_ijk[0].k] = IndS(k, i);
						d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
						d_ijk[0].coord3[d_ijk[0].k] = IndS(k, j);
						d_ijk[0].hash[d_ijk[0].k] = IND(IndS(k, i), IndS(i, j), IndS(k, j));
						d_ijk[0].k++;

						d_ijk[0].data[d_ijk[0].k].real = 1.0 / sqrt(2.0);
						d_ijk[0].coord1[d_ijk[0].k] = IndJ(k, j);
						d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
						d_ijk[0].coord3[d_ijk[0].k] = IndJ(k, i);
						d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(k, j), IndS(i, j), IndJ(k, i));
						d_ijk[0].k++;

						d_ijk[0].data[d_ijk[0].k].real = 1.0 / sqrt(2.0);
						d_ijk[0].coord1[d_ijk[0].k] = IndJ(k, i);
						d_ijk[0].coord2[d_ijk[0].k] = IndS(i, j);
						d_ijk[0].coord3[d_ijk[0].k] = IndJ(k, j);
						d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(k, i), IndS(i, j), IndJ(k, j));
						d_ijk[0].k++;

						d_ijk[0].data[d_ijk[0].k].real = 1.0 / sqrt(2.0);
						d_ijk[0].coord1[d_ijk[0].k] = IndJ(k, j);
						d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
						d_ijk[0].coord3[d_ijk[0].k] = IndS(k, i);
						d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(k, j), IndJ(i, j), IndS(k, i));
						d_ijk[0].k++;

						d_ijk[0].data[d_ijk[0].k].real = -1.0 / sqrt(2.0);
						d_ijk[0].coord1[d_ijk[0].k] = IndJ(k, i);
						d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
						d_ijk[0].coord3[d_ijk[0].k] = IndS(k, j);
						d_ijk[0].hash[d_ijk[0].k] = IND(IndJ(k, i), IndJ(i, j), IndS(k, j));
						d_ijk[0].k++;

						d_ijk[0].data[d_ijk[0].k].real = -1.0 / sqrt(2.0);
						d_ijk[0].coord1[d_ijk[0].k] = IndS(k, j);
						d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
						d_ijk[0].coord3[d_ijk[0].k] = IndJ(k, i);
						d_ijk[0].hash[d_ijk[0].k] = IND(IndS(k, j), IndJ(i, j), IndJ(k, i));
						d_ijk[0].k++;

						d_ijk[0].data[d_ijk[0].k].real = 1.0 / sqrt(2.0);
						d_ijk[0].coord1[d_ijk[0].k] = IndS(k, i);
						d_ijk[0].coord2[d_ijk[0].k] = IndJ(i, j);
						d_ijk[0].coord3[d_ijk[0].k] = IndJ(k, j);
						d_ijk[0].hash[d_ijk[0].k] = IND(IndS(k, i), IndJ(i, j), IndJ(k, j));
						d_ijk[0].k++;
					}
				}
			}
		}
	}


	for (int j = 2; j < N; j++)
	{
		int i = 1;
		d_ijk[0].data[d_ijk[0].k].real = 2.0 / sqrt((double)(j) * (j + 1));
		d_ijk[0].coord1[d_ijk[0].k] = IndD(j);
		d_ijk[0].coord2[d_ijk[0].k] = IndD(i);
		d_ijk[0].coord3[d_ijk[0].k] = IndD(i);
		d_ijk[0].hash[d_ijk[0].k] = IND(IndD(j), IndD(i), IndD(i));
		d_ijk[0].k++;

		d_ijk[0].data[d_ijk[0].k].real = 2.0 / sqrt((double)(j) * (j + 1));
		d_ijk[0].coord1[d_ijk[0].k] = IndD(i);
		d_ijk[0].coord2[d_ijk[0].k] = IndD(i);
		d_ijk[0].coord3[d_ijk[0].k] = IndD(j);
		d_ijk[0].hash[d_ijk[0].k] = IND(IndD(i), IndD(i), IndD(j));
		d_ijk[0].k++;

		d_ijk[0].data[d_ijk[0].k].real = 2.0 / sqrt((double)(j) * (j + 1));
		d_ijk[0].coord1[d_ijk[0].k] = IndD(i);
		d_ijk[0].coord2[d_ijk[0].k] = IndD(j);
		d_ijk[0].coord3[d_ijk[0].k] = IndD(i);
		d_ijk[0].hash[d_ijk[0].k] = IND(IndD(i), IndD(j), IndD(i));
		d_ijk[0].k++;
	}

	for (int i = 2; i < N; i++)
	{
		d_ijk[0].data[d_ijk[0].k].real = 2.0 * (1 - i) / sqrt((double)(i) * (i + 1));
		d_ijk[0].coord1[d_ijk[0].k] = IndD(i);
		d_ijk[0].coord2[d_ijk[0].k] = IndD(i);
		d_ijk[0].coord3[d_ijk[0].k] = IndD(i);
		d_ijk[0].hash[d_ijk[0].k] = IND(IndD(i), IndD(i), IndD(i));
		d_ijk[0].k++;

		for (int j = i + 1; j < N; j++)
		{
			d_ijk[0].data[d_ijk[0].k].real = 2.0 / sqrt((double)(j) * (j + 1));
			d_ijk[0].coord1[d_ijk[0].k] = IndD(j);
			d_ijk[0].coord2[d_ijk[0].k] = IndD(i);
			d_ijk[0].coord3[d_ijk[0].k] = IndD(i);
			d_ijk[0].hash[d_ijk[0].k] = IND(IndD(j), IndD(i), IndD(i));
			d_ijk[0].k++;

			d_ijk[0].data[d_ijk[0].k].real = 2.0 / sqrt((double)(j) * (j + 1));
			d_ijk[0].coord1[d_ijk[0].k] = IndD(i);
			d_ijk[0].coord2[d_ijk[0].k] = IndD(i);
			d_ijk[0].coord3[d_ijk[0].k] = IndD(j);
			d_ijk[0].hash[d_ijk[0].k] = IND(IndD(i), IndD(i), IndD(j));
			d_ijk[0].k++;

			d_ijk[0].data[d_ijk[0].k].real = 2.0 / sqrt((double)(j) * (j + 1));
			d_ijk[0].coord1[d_ijk[0].k] = IndD(i);
			d_ijk[0].coord2[d_ijk[0].k] = IndD(j);
			d_ijk[0].coord3[d_ijk[0].k] = IndD(i);
			d_ijk[0].hash[d_ijk[0].k] = IND(IndD(i), IndD(j), IndD(i));
			d_ijk[0].k++;
		}
	}
}
void print_matrix(Tensor_Coordinates * matrix, FILE * f)
{
	for (unsigned int i = 0; i < matrix[0].k; i++)
	{
		fprintf(f, "%2.3i %2.3i %2.3i %lf\n", matrix[0].coord1[i] + 1, matrix[0].coord2[i] + 1, matrix[0].coord3[i] + 1, matrix[0].data[i].real);
	}
}
void free_matrix(Tensor_Coordinates * matrix)
{
	delete(matrix[0].data);
	delete(matrix[0].coord1);
	delete(matrix[0].coord2);
	delete(matrix[0].coord3);
	delete(matrix[0].hash);
	matrix[0].k = 0;
}

void init_f_d(Model *m)
{
	int N = m->N;
	int N_mat = m->N_mat;

	FMatrixs * Fs = m->Fs;

	crsMatrix * Fjk, *Fkj, *Fk, *Fsum, *Fsub, *Fres;

	dcomplex add, sub;
	add.re = 1.0; add.im = 0.0;
	sub.re = -1.0; sub.im = 0.0;

	vector<pair<int, dcomplex> > * f_std;
	vector<pair<int, dcomplex> > * d_std;

	m->d_mat = new crsMatrix *[N_mat];
	m->f_mat = new crsMatrix *[N_mat];

	int i, j, k, cnt;

	Fjk = new crsMatrix;
	Fkj = new crsMatrix;

	Fsum = new crsMatrix;
	Fsub = new crsMatrix;

	Fres = new crsMatrix;

	for (i = 0; i < N_mat; i++)
	{
		f_std = new vector<pair<int, dcomplex> >[N_mat];
		d_std = new vector<pair<int, dcomplex> >[N_mat];

		for (j = 0; j < N_mat; j++)
		{
			for (k = 0; k < N_mat; k++)
			{

				Fk = new crsMatrix(*(Fs->F[k + 1]));

				SparseMKLMult(*(Fs->F[j + 1]), *Fk, *Fjk, true);
				SparseMKLMult(*Fk, *(Fs->F[j + 1]), *Fkj, true);
				delete Fk;

				SparseMKLAdd(*Fjk, add, *Fkj, *Fsum, true);
				SparseMKLAdd(*Fjk, sub, *Fkj, *Fsub, true);


				SparseMKLMult(*(Fs->F[i + 1]), *Fsub, *Fres, true);
				cnt = trace_struct(*Fres);
				if (cnt > 0)
				{
					dcomplex val, tr;
					tr = trace(*Fres);
					val.re = tr.im;
					val.im = -tr.re;
					f_std[j].push_back(make_pair(k, val));
				}

				SparseMKLMult(*(Fs->F[i + 1]), *Fsum, *Fres, true);
				cnt = trace_struct(*Fres);
				if (cnt > 0)
				{
					dcomplex val;
					val = trace(*Fres);
					d_std[j].push_back(make_pair(k, val));
				}
			}
		}

		m->f_mat[i] = stdToCrs(f_std, N_mat);
		m->d_mat[i] = stdToCrs(d_std, N_mat);

		//    printMatrixVal(m->f_mat[i]);
		//    printMatrixVal(m->d_mat[i]);

		delete[] f_std;
		delete[] d_std;
	}
	delete Fjk;
	delete Fkj;

	delete Fsum;
	delete Fsub;

	delete Fres;

}
crsMatrix * TensorToCrs(Tensor_Coordinates &t_ijk, int s, int f, int N)
{
	crsMatrix *  mat = new crsMatrix(N, f - s);
	int i = 0;
	dcomplex val;

	for (i = 0; i < f - s; i++)
	{
		val.re = t_ijk.data[s + i].real;
		val.im = t_ijk.data[s + i].imag;
		mat->Value[i] = val;
		mat->Col[i] = t_ijk.coord3[s + i];
	}

	for (i = 0; i < N + 1; i++)
	{
		mat->RowIndex[i] = 0;
	}

	for (i = 0; i < f - s; i++)
	{
		mat->RowIndex[t_ijk.coord2[s + i] + 1]++;
	}

	int cnt = 0;

	for (i = 0; i < N; i++)
	{
		mat->RowIndex[i] = cnt;
		cnt += mat->RowIndex[i + 1];
	}
	mat->RowIndex[i] = cnt;

	return mat;
}
crsMatrix * TensorToCrs(Tensor_Coordinates_1 &t_ijk, int s, int f, int N)
{
	crsMatrix *  mat = new crsMatrix(N, f - s);
	int i = 0;
	dcomplex val;

	for (i = 0; i < f - s; i++)
	{
		val.re = t_ijk.data[s + i].real;
		val.im = t_ijk.data[s + i].imag;
		mat->Value[i] = val;
		mat->Col[i] = t_ijk.coord3[s + i];
	}

	for (i = 0; i < N + 1; i++)
	{
		mat->RowIndex[i] = 0;
	}

	for (i = 0; i < f - s; i++)
	{
		mat->RowIndex[t_ijk.coord2[s + i] + 1]++;
	}

	int cnt = 0;

	for (i = 0; i < N; i++)
	{
		mat->RowIndex[i] = cnt;
		cnt += mat->RowIndex[i + 1];
	}
	mat->RowIndex[i] = cnt;

	return mat;
}
void init_f_d_valentin(Model *m)
{
	int N = m->N;
	int N_mat = m->N_mat;
	int i;

	m->d_mat = new crsMatrix *[N_mat];
	m->f_mat = new crsMatrix *[N_mat];

	//  vector<pair<int , dcomplex> > * f_std;
	//  vector<pair<int , dcomplex> > * d_std;

	Tensor_Coordinates f_ijk;
	Tensor_Coordinates d_ijk;
	Tensor_Coordinates_1 * f_1 = new Tensor_Coordinates_1[N_mat];
	Tensor_Coordinates_1 * d_1 = new Tensor_Coordinates_1[N_mat];

	fijk_coord(&f_ijk, N + 1);
	sort_matrix(&f_ijk, f_1, N_mat);
	//  sort_matrix(&f_ijk);

	dijk_coord(&d_ijk, N + 1);
	sort_matrix(&d_ijk, d_1, N_mat);
	//  sort_matrix(&d_ijk);

	int k, s, f;

	//  k = 0;
	for (i = 0; i < N_mat; i++)
	{
		//    s = k;
		//    while(f_ijk.coord1[k] == i)
		//    {
		//      k++;
		//    }
		//    f = k;
		//    m->f_mat[i] = TensorToCrs(f_ijk, s, f, N_mat);
		m->f_mat[i] = TensorToCrs(f_1[i], 0, f_1[i].N, N_mat);
	}

	//  k = 0;
	for (i = 0; i < N_mat; i++)
	{
		//    s = k;
		//    while(d_ijk.coord1[k] == i)
		//    {
		//      k++;
		//    }
		//    f = k;
		//    m->d_mat[i] = TensorToCrs(d_ijk, s, f, N_mat);
		m->d_mat[i] = TensorToCrs(d_1[i], 0, d_1[i].N, N_mat);
	}

	for (i = 0; i < N_mat; i++)
	{
		delete[] f_1[i].coord2;
		delete[] f_1[i].coord3;
		delete[] f_1[i].data;
		delete[] d_1[i].coord2;
		delete[] d_1[i].coord3;
		delete[] d_1[i].data;
	}
	delete[] f_1;
	delete[] d_1;

	free_matrix(&f_ijk);
	free_matrix(&d_ijk);
}

void transpFs(Model *m)
{
	int N_mat = m->N_mat;
	crsMatrix **f_mat = m->f_mat;
	crsMatrix **f_H_mat = new crsMatrix*[N_mat];
	for (int i = 0; i < N_mat; i++)
	{
		f_H_mat[i] = new crsMatrix(*(f_mat[i]));
		Transpose(*(f_mat[i]), *(f_H_mat[i]), false);
		//    printMatrixVal(f_H_mat[i]);
	}
	m->f_H_mat = f_H_mat;
}

void calc_Q_0(Model * m)
{
	int N_mat = m->N_mat;
	dcomplex * h_0 = m->h_0;

	crsMatrix ** f_mat = m->f_H_mat;
	crsMatrix * Q_0 = new crsMatrix(m->N_mat, 1);
	crsMatrix * resSum;

	Q_0->Col[0] = 0;
	Q_0->RowIndex[0] = 0;
	for (int i = 1; i <= N_mat; i++)
	{
		Q_0->RowIndex[i] = 0;
	}

	for (int i = 0; i < N_mat; i++)
	{
		if ((h_0[i].re != 0.0) || (h_0[i].im != 0.0))
		{
			resSum = new crsMatrix;
			SparseMKLAdd(*Q_0, h_0[i], *(f_mat[i]), *resSum);
			delete Q_0;
			Q_0 = resSum;
		}
	}

	m->Q_0 = Q_0;
}
void calc_Q_1(Model * m)
{
	int N_mat = m->N_mat;
	dcomplex * h_1 = m->h_1;

	crsMatrix ** f_mat = m->f_H_mat;
	crsMatrix * Q_1 = new crsMatrix(m->N_mat, 1);
	crsMatrix * resSum;

	Q_1->Col[0] = 0;
	Q_1->RowIndex[0] = 0;
	for (int i = 1; i <= N_mat; i++)
	{
		Q_1->RowIndex[i] = 0;
	}

	for (int i = 0; i < N_mat; i++)
	{
		if ((h_1[i].re != 0.0) || (h_1[i].im != 0.0))
		{
			resSum = new crsMatrix;
			SparseMKLAdd(*Q_1, h_1[i], *(f_mat[i]), *resSum);
			delete Q_1;
			Q_1 = resSum;
		}
	}

	m->Q_1 = Q_1;
}

void calcKs(Model *m)
{
	int N = m->N;
	int N_mat = m->N_mat;
	crsMatrix *Ks_tmp;
	crsMatrix *FsT;
	dcomplex  *Ks = m->Ks;
	crsMatrix *As = m->a_mat;
	crsMatrix **Fs = m->f_mat;
	crsMatrix *AsT;

	for (int i = 0; i < N_mat; i++)
	{
		Ks[i].re = 0.0;
		Ks[i].im = 0.0;
	}

	//  printMatrixVal(As);
	for (int i = 0; i < N_mat; i++)
	{
		AsT = As;
		//AsT = new crsMatrix(*(As));
		//Transpose(*(As), *AsT);
		FsT = Fs[i];
		FsT = new crsMatrix(*(Fs[i]));
		Transpose(*(Fs[i]), *FsT, false);
		//    printMatrixVal(FsT);
		for (int j = 0; j < N_mat; j++)
		{
			int ii, jj, m1, m2;
			ii = As->RowIndex[i];
			m1 = As->RowIndex[i + 1];
			jj = FsT->RowIndex[j];
			m2 = FsT->RowIndex[j + 1];

			while ((ii < m1) && (jj < m2))
			{
				if (AsT->Col[ii] < FsT->Col[jj])
				{
					ii++;
				}
				else
				{
					if (AsT->Col[ii] > FsT->Col[jj])
					{
						jj++;
					}
					else
					{
						dcomplex as, fs;
						as = AsT->Value[ii];
						fs = FsT->Value[jj];
						Ks[j].re += as.re * fs.re - as.im * fs.im;
						Ks[j].im += as.re * fs.im + as.im * fs.re;

						ii++;
						jj++;
					}
				}
			}

		}
		delete FsT;
		//delete AsT;
	}

	dcomplex val;
	for (int i = 0; i < N_mat; i++)
	{
		Ks[i].re *= (-1.0) / (N + 1);
		Ks[i].im *= (1.0) / (N + 1);
		val = Ks[i];
		Ks[i].re = val.im;
		Ks[i].im = val.re;
	}
}
void calcRs(Model *m)
{
	dcomplex sum;
	sum.re = 1.0;
	sum.im = 0.0;

	int N_mat = m->N_mat;

	crsMatrix * Rs_tmp;
	crsMatrix ** RsTh;
	crsMatrix * Rs = new crsMatrix(m->N_mat, 1);

	Rs->Col[0] = 1;
	Rs->RowIndex[0] = 1;
	for (int i = 1; i <= N_mat; i++)
	{
		Rs->RowIndex[i] = 1;
	}

	int nThread = 1;
	//  omp_set_num_threads(3);
#pragma omp parallel
	{
#pragma omp single
		nThread = omp_get_num_threads();
	}
	printf("nThread %d \n", nThread);

	int SubTask = 1;
	int CntTask = nThread * SubTask;

	RsTh = new crsMatrix *[CntTask];

	crsMatrix ** f_mat = m->f_mat;
	crsMatrix ** f_H_mat = m->f_H_mat;
	crsMatrix ** d_mat = m->d_mat;

	for (int i = 0; i < N_mat; i++)
	{
		toOneBase(*(f_mat[i]));
		toOneBase(*(d_mat[i]));
		toOneBase(*(f_H_mat[i]));
	}

	int size = N_mat / CntTask;
	//  int *start = new int[CntTask];
	//  int *finish = new int[CntTask];

	//  start[0] = 0;
	//  finish[0] = size;
	//  if((N_mat % CntTask) != 0) finish[0]++;
	printf("N_mat %d \n", N_mat);
	//  printf("%d %d \n", start[0], finish[0]);
	//  for(int i = 1; i < CntTask; i++)
	//  {
	//    start[i] = finish[i - 1];
	//    finish[i] = start[i] + size;
	//    if((N_mat % CntTask) >  i)finish[i]++;
	//    printf("%d %d \n", start[i], finish[i]);
	//  }


	int id = 0;
	int portion = 100;
#pragma omp parallel
	{
		int tid = omp_get_thread_num();
		int localTskID;
		while (true)
		{
#pragma omp critical
			{
				localTskID = id;
				id++;
			}
			int start = localTskID * portion;
			int finish = (localTskID + 1) * portion;
			if (start > N_mat) break;
			if (finish > N_mat) {
				finish = N_mat;
				printf("%d %d \n", start, finish);
			}
			RsTh[tid] = calcSubRs(m, start, finish);
#pragma omp critical
			{
				Rs_tmp = new crsMatrix;
				SparseMKLAddOne(*Rs, sum, *RsTh[tid], *Rs_tmp);
				delete RsTh[tid];
				delete Rs;
				Rs = Rs_tmp;
			}
		}

	}
	//delete [] start;
	//delete [] finish;

	//  for(int i = 0; i < CntTask; i++)
	//  {
	//      Rs_tmp = new crsMatrix;
	//      SparseMKLAddOne(*Rs, sum, *RsTh[i], *Rs_tmp);
	//      delete RsTh[i];
	//      delete Rs;
	//      Rs = Rs_tmp;
	//  }
	delete[] RsTh;

	//Rs = calcSubRs(m, 0, N_mat);
	toZeroBase(*Rs);
	for (int i = 0; i < N_mat; i++)
	{
		toZeroBase(*(f_mat[i]));
		toZeroBase(*(d_mat[i]));
		toZeroBase(*(f_H_mat[i]));
	}
	//  delete FkT;
	//  delete FiT;

	double mm = -1.0 / 4.0;
	for (int i = 0; i < Rs->NZ; i++)
	{
		Rs->Value[i].re *= mm;
		Rs->Value[i].im *= mm;
	}
	m->Rs = Rs;
}
crsMatrix* calcSubRs(Model *m, int start, int finish)
{
	int N_mat = m->N_mat;
	dcomplex sum_i;
	sum_i.re = 0.0;
	sum_i.im = 1.0;

	dcomplex sum;
	sum.re = 1.0;
	sum.im = 0.0;

	crsMatrix * tmp;
	crsMatrix * MatSDi, *MatSDk;
	crsMatrix * M1, *M2, *Rsum;

	crsMatrix * Rs_tmp;
	crsMatrix * Rs = new crsMatrix(m->N_mat, 1);

	Rs->Col[0] = 1;
	Rs->RowIndex[0] = 1;
	for (int i = 1; i <= N_mat; i++)
	{
		Rs->RowIndex[i] = 1;
	}

	crsMatrix ** f_mat = m->f_mat;
	crsMatrix ** f_H_mat = m->f_H_mat;
	crsMatrix ** d_mat = m->d_mat;
	crsMatrix *  As = m->a_mat;

	M1 = new crsMatrix;
	M2 = new crsMatrix;
	MatSDi = new crsMatrix;
	MatSDk = new crsMatrix;
	Rsum = new crsMatrix;

	//  for(int i = 0; i < N_mat; i++)
	//  {
	//    toOneBase(*(f_mat[i]));
	//    toOneBase(*(d_mat[i]));
	//    toOneBase(*(f_H_mat[i]));
	//  }

	for (int i = start; i < finish; i++)
	{
		SparseMKLAddOne(*(f_mat[i]), sum_i, *(d_mat[i]), *MatSDi, true);

		int st = As->RowIndex[i];
		int fn = As->RowIndex[i + 1];

		for (int j = st; j < fn; j++)
		{
			int k = As->Col[j];
			SparseMKLAddOne(*(f_mat[k]), sum_i, *(d_mat[k]), *MatSDk, true);

			for (int conj_i = 0; conj_i < MatSDk->NZ; conj_i++)
				MatSDk->Value[conj_i].im = -MatSDk->Value[conj_i].im;

			SparseMKLMultOne(*(f_H_mat[k]), *MatSDi, *M1, true);
			SparseMKLMultOne(*(f_H_mat[i]), *MatSDk, *M2, true);

			SparseMKLAddOne(*M1, sum, *M2, *Rsum, true);

			dcomplex val1 = As->Value[j];
			for (int ii = 0; ii < Rsum->NZ; ii++)
			{
				dcomplex val2 = Rsum->Value[ii];
				Rsum->Value[ii].re = val1.re *  val2.re - val1.im *  val2.im;
				Rsum->Value[ii].im = val1.re *  val2.im + val1.im *  val2.re;
			}

			Rs_tmp = new crsMatrix;
			SparseMKLAddOne(*Rs, sum, *Rsum, *Rs_tmp);
			delete Rs;
			Rs = Rs_tmp;
		}
	}

	//  toZeroBase(*Rs);
	//  for(int i = 0; i < N_mat; i++)
	//  {
	//    toZeroBase(*(f_mat[i]));
	//    toZeroBase(*(d_mat[i]));
	//    toZeroBase(*(f_H_mat[i]));
	//  }
	//  delete FkT;
	//  delete FiT;

	delete Rsum;

	delete MatSDi;
	delete MatSDk;

	delete M1;
	delete M2;

	//double mm = -m->conf.g / 4.0 / m->N;
	//for(int i = 0; i < Rs->NZ; i++)
	//{
	//  Rs->Value[i].re *= mm;
	//  Rs->Value[i].im *= mm;
	//}
	return Rs;
}

void calc_G_0_s(Model *m, ConfigParam &cp)
{
	int N_mat = m->N_mat;
	crsMatrix * G_0_s = new crsMatrix();

	crsMatrix * subSum = new crsMatrix();

	dcomplex sum_0;
	sum_0.re = cp.drv_ampl;
	sum_0.im = 0.0;

	dcomplex sum_1;
	sum_1.re = 1.0;
	sum_1.im = 0.0;

	SparseMKLAdd(*(m->Q_0), sum_0, *(m->Q_1), *subSum);
	SparseMKLAdd(*(subSum), sum_1, *(m->Rs), *G_0_s);

	m->G_0_s = G_0_s;

	delete subSum;
}
void calc_G_1_s(Model *m, ConfigParam &cp)
{
	int N_mat = m->N_mat;
	crsMatrix * G_1_s = new crsMatrix();

	dcomplex sum;
	sum.re = 1.0;
	sum.im = 0.0;

	SparseMKLAdd(*(m->Q_0), sum, *(m->Rs), *G_1_s);

	m->G_1_s = G_1_s;
}

void init_conditions(dcomplex *mtx, RunParam &rp, ConfigParam &cp, MainData &md)
{
	for (int s_id_1 = 0; s_id_1 < md.size; s_id_1++)
	{
		for (int s_id_2 = 0; s_id_2 < md.size; s_id_2++)
		{
			mtx[s_id_1 * md.size + s_id_2].re = 0.0;
			mtx[s_id_1 * md.size + s_id_2].im = 0.0;
		}
	}

	if (cp.int_ist == 0)
	{
		mtx[cp.int_isi * md.size + cp.int_isi].re = 1.0;
		mtx[cp.int_isi * md.size + cp.int_isi].im = 0.0;
	}
	else if (cp.int_ist == 1)
	{
		for (int s_id_1 = 0; s_id_1 < md.size; s_id_1++)
		{
			mtx[s_id_1 * md.size + s_id_1].re = 1.0 / md.size;
			mtx[s_id_1 * md.size + s_id_1].im = 0.0;
		}
	}
	else
	{
		stringstream msg;
		msg << "wrong int_ist value: " << cp.int_ist << endl;
		Error(msg.str());
	}
}

void multMatVec(crsMatrix *mat, dcomplex * x, dcomplex * res)
{
	int i, j, s, f;
	for (i = 0; i < mat->N; i++)
	{
		s = mat->RowIndex[i];
		f = mat->RowIndex[i + 1];
		res[i].re = 0.0;
		res[i].im = 0.0;
		for (j = s; j < f; j++)
		{
			dcomplex v1 = mat->Value[j];
			dcomplex v2 = x[mat->Col[j]];
			res[i].re += v1.re * v2.re;// - v1.im * v2.im;
			res[i].im = 0.0;//+= v1.re * v2.im + v1.im * v2.re;
		}
	}
}

void calcVectValue_t0(double h, Model * m, dcomplex *x, dcomplex * res, dcomplex * tmp)
{
	int i;
	int N_mat = m->N_mat;
	crsMatrix * G_0_s = m->G_0_s;
	dcomplex  * Ks = m->Ks;

	multMatVec(G_0_s, x, tmp);

	for (i = 0; i < N_mat; i++)
	{
		res[i].re = (tmp[i].re - Ks[i].re) * h;
		res[i].im = 0.0;
	}
}

void calcVectValue_t1(double h, Model * m, dcomplex *x, dcomplex * res, dcomplex * tmp)
{
	int i;
	int N_mat = m->N_mat;
	crsMatrix * G_1_s = m->G_1_s;
	dcomplex  * Ks = m->Ks;

	multMatVec(G_1_s, x, tmp);

	for (i = 0; i < N_mat; i++)
	{
		res[i].re = (tmp[i].re - Ks[i].re) * h;
		res[i].im = 0.0;
	}
}

void initRhoODE(Model *m, RunParam &rp, ConfigParam &cp, MainData &md)
{
	int N = m->N;
	int N_mat = m->N_mat;
	dcomplex  * RhoF = m->RhoF;
	dcomplex * psi = new dcomplex[(N + 1)*(N + 1)];
	init_conditions(psi, rp, cp, md);

	if (rp.debug == 1)
	{
		string fn = rp.path + "rho_ini" + file_name_suffix(cp, 4);
		save_complex_data(fn, (MKL_Complex16 *)psi, (N + 1)*(N + 1), 16, false);
	}

	for (int i = 0; i < N_mat; i++)
	{
		RhoF[i].re = 0.0;
		RhoF[i].im = 0.0;
	}
	int k = 0;
	double val = 1.0 / sqrt(2.0);
	for (int i = 0; i < N + 1; i++)
	{
		for (int j = i + 1; j < N + 1; j++)
		{
			RhoF[k].re += psi[(i)* (N + 1) + (j)].re * val;
			RhoF[k].re += psi[(j)* (N + 1) + (i)].re * val;
			k++;

			RhoF[k].re += psi[(i)* (N + 1) + (j)].im * (-val);
			RhoF[k].re += psi[(j)* (N + 1) + (i)].im * (+val);
			k++;
		}
	}

	for (int i = 0; i < N; i++)
	{
		val = 1.0 / sqrt((double)((i + 1)* (i + 2)));
		for (int j = 0; j <= i; j++)
		{
			RhoF[k].re += psi[j * (N + 1) + j].re *val;
		}
		RhoF[k].re -= psi[(i + 1) * (N + 1) + (i + 1)].re * val * (i + 1);
		k++;
	}

	delete[] psi;
}

dcomplex calcDiffIter(Model *m)
{
	dcomplex max_diff, diff;
	max_diff.re = 0.0;
	max_diff.im = 0.0;
	for (int i = 0; i < m->N_mat; i++)
	{
		diff.re = fabs(m->prevRhoF[i].re - m->RhoF[i].re);
		diff.im = fabs(m->prevRhoF[i].im - m->RhoF[i].im);
		if (max_diff.re < diff.re)max_diff.re = diff.re;
		if (max_diff.im < diff.im)max_diff.im = diff.im;
	}

	return max_diff;
}

void calcODE_trans(Model *m, RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd)
{
	int i;
	int N_mat = m->N_mat;
	dcomplex * RhoF = m->RhoF;

	double eps = 1.0e-10;

	for (i = 0; i < N_mat; i++)
	{
		m->prevRhoF[i] = RhoF[i];
	}

	dcomplex * k1 = pd.k1;
	dcomplex * k2 = pd.k2;
	dcomplex * k3 = pd.k3;
	dcomplex * k4 = pd.k4;
	dcomplex * val = pd.val;
	dcomplex * tmp = pd.tmp;

	double step_t_0 = md.step_t_0;
	double step_t_1 = md.step_t_1;


	for (int period = 1; period <= cp.num_periods_trans; period++)
	{
		for (int t_0_step_id = 0; t_0_step_id < cp.num_steps_t_0; t_0_step_id++)
		{
			calcVectValue_t0(step_t_0, m, RhoF, k1, tmp);
			for (i = 0; i < N_mat; i++)
			{
				val[i].re = RhoF[i].re + k1[i].re / 2.0;
				val[i].im = RhoF[i].im + k1[i].im / 2.0;
			}
			calcVectValue_t0(step_t_0, m, val, k2, tmp);
			for (i = 0; i < N_mat; i++)
			{
				val[i].re = RhoF[i].re + k2[i].re / 2.0;
				val[i].im = RhoF[i].im + k2[i].im / 2.0;
			}
			calcVectValue_t0(step_t_0, m, val, k3, tmp);
			for (i = 0; i < N_mat; i++)
			{
				val[i].re = RhoF[i].re + k3[i].re;
				val[i].im = RhoF[i].im + k3[i].im;
			}
			calcVectValue_t0(step_t_0, m, val, k4, tmp);

			for (i = 0; i < N_mat; i++)
			{
				RhoF[i].re += (k1[i].re + 2.0 * k2[i].re + 2.0 * k3[i].re + k4[i].re) / 6.0;
				RhoF[i].im += (k1[i].im + 2.0 * k2[i].im + 2.0 * k3[i].im + k4[i].im) / 6.0;
			}
		}

		for (int t_1_step_id = 0; t_1_step_id < cp.num_steps_t_1; t_1_step_id++)
		{
			calcVectValue_t1(step_t_1, m, RhoF, k1, tmp);
			for (i = 0; i < N_mat; i++)
			{
				val[i].re = RhoF[i].re + k1[i].re / 2.0;
				val[i].im = RhoF[i].im + k1[i].im / 2.0;
			}
			calcVectValue_t1(step_t_1, m, val, k2, tmp);
			for (i = 0; i < N_mat; i++)
			{
				val[i].re = RhoF[i].re + k2[i].re / 2.0;
				val[i].im = RhoF[i].im + k2[i].im / 2.0;
			}
			calcVectValue_t1(step_t_1, m, val, k3, tmp);
			for (i = 0; i < N_mat; i++)
			{
				val[i].re = RhoF[i].re + k3[i].re;
				val[i].im = RhoF[i].im + k3[i].im;
			}
			calcVectValue_t1(step_t_1, m, val, k4, tmp);

			for (i = 0; i < N_mat; i++)
			{
				RhoF[i].re += (k1[i].re + 2.0 * k2[i].re + 2.0 * k3[i].re + k4[i].re) / 6.0;
				RhoF[i].im += (k1[i].im + 2.0 * k2[i].im + 2.0 * k3[i].im + k4[i].im) / 6.0;
			}
		}
	}
}

void calcODE_std(Model *m, RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd)
{
	int i;
	int N_mat = m->N_mat;
	dcomplex * RhoF = m->RhoF;

	double eps = 1.0e-10;

	for (i = 0; i < N_mat; i++)
	{
		m->prevRhoF[i] = RhoF[i];
	}

	dcomplex * k1 = pd.k1;
	dcomplex * k2 = pd.k2;
	dcomplex * k3 = pd.k3;
	dcomplex * k4 = pd.k4;
	dcomplex * val = pd.val;
	dcomplex * tmp = pd.tmp;

	double step_t_0 = md.step_t_0;
	double step_t_1 = md.step_t_1;

	double curr_time = 0.0;
	int dump_id = 0;

	calcRho(m);

	characteristics_std(m, rp, cp, md, pd, dump_id);
	dump_id++;

	for (int period = 1; period <= cp.num_periods_obser; period++)
	{
		for (int t_0_step_id = 0; t_0_step_id < cp.num_steps_t_0; t_0_step_id++)
		{
			calcVectValue_t0(step_t_0, m, RhoF, k1, tmp);
			for (i = 0; i < N_mat; i++)
			{
				val[i].re = RhoF[i].re + k1[i].re / 2.0;
				val[i].im = RhoF[i].im + k1[i].im / 2.0;
			}
			calcVectValue_t0(step_t_0, m, val, k2, tmp);
			for (i = 0; i < N_mat; i++)
			{
				val[i].re = RhoF[i].re + k2[i].re / 2.0;
				val[i].im = RhoF[i].im + k2[i].im / 2.0;
			}
			calcVectValue_t0(step_t_0, m, val, k3, tmp);
			for (i = 0; i < N_mat; i++)
			{
				val[i].re = RhoF[i].re + k3[i].re;
				val[i].im = RhoF[i].im + k3[i].im;
			}
			calcVectValue_t0(step_t_0, m, val, k4, tmp);

			for (i = 0; i < N_mat; i++)
			{
				RhoF[i].re += (k1[i].re + 2.0 * k2[i].re + 2.0 * k3[i].re + k4[i].re) / 6.0;
				RhoF[i].im += (k1[i].im + 2.0 * k2[i].im + 2.0 * k3[i].im + k4[i].im) / 6.0;
			}
		}

		for (int t_1_step_id = 0; t_1_step_id < cp.num_steps_t_1; t_1_step_id++)
		{
			calcVectValue_t1(step_t_1, m, RhoF, k1, tmp);
			for (i = 0; i < N_mat; i++)
			{
				val[i].re = RhoF[i].re + k1[i].re / 2.0;
				val[i].im = RhoF[i].im + k1[i].im / 2.0;
			}
			calcVectValue_t1(step_t_1, m, val, k2, tmp);
			for (i = 0; i < N_mat; i++)
			{
				val[i].re = RhoF[i].re + k2[i].re / 2.0;
				val[i].im = RhoF[i].im + k2[i].im / 2.0;
			}
			calcVectValue_t1(step_t_1, m, val, k3, tmp);
			for (i = 0; i < N_mat; i++)
			{
				val[i].re = RhoF[i].re + k3[i].re;
				val[i].im = RhoF[i].im + k3[i].im;
			}
			calcVectValue_t1(step_t_1, m, val, k4, tmp);

			for (i = 0; i < N_mat; i++)
			{
				RhoF[i].re += (k1[i].re + 2.0 * k2[i].re + 2.0 * k3[i].re + k4[i].re) / 6.0;
				RhoF[i].im += (k1[i].im + 2.0 * k2[i].im + 2.0 * k3[i].im + k4[i].im) / 6.0;
			}
		}

		curr_time = period * md.T;

		if (period == pd.dump_periods[dump_id])
		{

			if (rp.ipp == 1)
			{
				cout << endl << "Dump period: " << period << endl;
			}

			calcRho(m);
			characteristics_std(m, rp, cp, md, pd, dump_id);

			dump_id++;
		}
	}
}

void calcODE_deep(Model *m, RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd)
{
	string fn;

	int dupm_step_t_0 = cp.num_steps_t_0 / cp.int_dn;
	int dupm_step_t_1 = cp.num_steps_t_1 / cp.int_dn;

	int i;
	int N_mat = m->N_mat;
	dcomplex * RhoF = m->RhoF;

	double eps = 1.0e-10;

	for (i = 0; i < N_mat; i++)
	{
		m->prevRhoF[i] = RhoF[i];
	}

	dcomplex * k1 = pd.k1;
	dcomplex * k2 = pd.k2;
	dcomplex * k3 = pd.k3;
	dcomplex * k4 = pd.k4;
	dcomplex * val = pd.val;
	dcomplex * tmp = pd.tmp;

	double step_t_0 = md.step_t_0;
	double step_t_1 = md.step_t_1;

	double curr_time = 0.0;
	int dump_id = 0;

	calcRho(m);

	if (rp.debug == 1)
	{
		fn = rp.path + "rho" + "_" + to_string(dump_id) + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, m->Rho, 16, false);
	}

	characteristics_deep(m, rp, cp, md, pd, dump_id);
	dump_id++;

	for (int period = 1; period <= cp.num_periods_obser; period++)
	{
		for (int t_0_step_id = 0; t_0_step_id < cp.num_steps_t_0; t_0_step_id++)
		{
			calcVectValue_t0(step_t_0, m, RhoF, k1, tmp);
			for (i = 0; i < N_mat; i++)
			{
				val[i].re = RhoF[i].re + k1[i].re / 2.0;
				val[i].im = RhoF[i].im + k1[i].im / 2.0;
			}
			calcVectValue_t0(step_t_0, m, val, k2, tmp);
			for (i = 0; i < N_mat; i++)
			{
				val[i].re = RhoF[i].re + k2[i].re / 2.0;
				val[i].im = RhoF[i].im + k2[i].im / 2.0;
			}
			calcVectValue_t0(step_t_0, m, val, k3, tmp);
			for (i = 0; i < N_mat; i++)
			{
				val[i].re = RhoF[i].re + k3[i].re;
				val[i].im = RhoF[i].im + k3[i].im;
			}
			calcVectValue_t0(step_t_0, m, val, k4, tmp);

			for (i = 0; i < N_mat; i++)
			{
				RhoF[i].re += (k1[i].re + 2.0 * k2[i].re + 2.0 * k3[i].re + k4[i].re) / 6.0;
				RhoF[i].im += (k1[i].im + 2.0 * k2[i].im + 2.0 * k3[i].im + k4[i].im) / 6.0;
			}

			if (rp.debug == 1)
			{
				calcRho(m);

				fn = rp.path + "rho" + "_" + to_string(dump_id) + file_name_suffix(cp, 4);
				save_sparse_complex_mtx(fn, m->Rho, 16, false);
			}

			if (t_0_step_id % dupm_step_t_0 == 0)
			{
				calcRho(m);
				characteristics_deep(m, rp, cp, md, pd, dump_id);
				dump_id++;
			}		
		}

		for (int t_1_step_id = 0; t_1_step_id < cp.num_steps_t_1; t_1_step_id++)
		{
			calcVectValue_t1(step_t_1, m, RhoF, k1, tmp);
			for (i = 0; i < N_mat; i++)
			{
				val[i].re = RhoF[i].re + k1[i].re / 2.0;
				val[i].im = RhoF[i].im + k1[i].im / 2.0;
			}
			calcVectValue_t1(step_t_1, m, val, k2, tmp);
			for (i = 0; i < N_mat; i++)
			{
				val[i].re = RhoF[i].re + k2[i].re / 2.0;
				val[i].im = RhoF[i].im + k2[i].im / 2.0;
			}
			calcVectValue_t1(step_t_1, m, val, k3, tmp);
			for (i = 0; i < N_mat; i++)
			{
				val[i].re = RhoF[i].re + k3[i].re;
				val[i].im = RhoF[i].im + k3[i].im;
			}
			calcVectValue_t1(step_t_1, m, val, k4, tmp);

			for (i = 0; i < N_mat; i++)
			{
				RhoF[i].re += (k1[i].re + 2.0 * k2[i].re + 2.0 * k3[i].re + k4[i].re) / 6.0;
				RhoF[i].im += (k1[i].im + 2.0 * k2[i].im + 2.0 * k3[i].im + k4[i].im) / 6.0;
			}

			if (rp.debug == 1)
			{
				calcRho(m);

				fn = rp.path + "rho" + "_" + to_string(dump_id) + file_name_suffix(cp, 4);
				save_sparse_complex_mtx(fn, m->Rho, 16, false);
			}

			if (t_1_step_id % dupm_step_t_1 == 0)
			{
				calcRho(m);
				characteristics_deep(m, rp, cp, md, pd, dump_id);
				dump_id++;
			}
		}

		curr_time = period * md.T;

		if (rp.ipp == 1)
		{
			cout << endl << "Dump period: " << period << endl;
		}
	}
}

void calcRho(Model *m)
{
	int N_mat = m->N_mat;
	int N = m->N;
	dcomplex  * RhoF = m->RhoF;
	crsMatrix * mRho = new crsMatrix(N + 1, (N + 1) * (N + 1));
	dcomplex  * Rho = mRho->Value;

	int i, j;
	for (i = 0; i < (N + 1); i++)
	{
		mRho->RowIndex[i] = i * (N + 1);
		for (j = 0; j < (N + 1); j++)
		{
			mRho->Col[i * (N + 1) + j] = j;
		}
	}
	mRho->RowIndex[i] = i * (N + 1);


	for (i = 0; i < (N + 1) * (N + 1); i++)
	{
		Rho[i].re = 0.0;
		Rho[i].im = 0.0;
	}
	for (i = 0; i < (N + 1); i++)
	{
		Rho[i * (N + 1) + i].re = 1.0 / (N + 1);
	}

	int k = 0;
	double val = 1.0 / sqrt(2.0);
	for (i = 0; i < N + 1; i++)
	{
		for (j = i + 1; j < N + 1; j++)
		{
			Rho[(i) * (N + 1) + (j)].re += RhoF[k].re * val;
			Rho[(j) * (N + 1) + (i)].re += RhoF[k].re * val;
			k++;

			Rho[(i) * (N + 1) + (j)].im += RhoF[k].re * (+val);
			Rho[(j) * (N + 1) + (i)].im += RhoF[k].re * (-val);
			k++;
		}
	}

	for (i = 0; i < N; i++)
	{
		val = 1.0 / sqrt((double)((i + 1)* (i + 2)));
		for (j = 0; j <= i; j++)
		{
			Rho[j * (N + 1) + j].re += RhoF[k].re * val;
		}
		Rho[j * (N + 1) + j].re -= RhoF[k].re * val * (i + 1);
		k++;
	}

	if (m->Rho != NULL)
	{
		delete m->Rho;
	}
	m->Rho = mRho;
}

void characteristics_std(Model *m, RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd, int dump_id)
{
	// Add here regular characteristics
}
void characteristics_deep(Model *m, RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd, int dump_id)
{
	MKL_Complex16 * rho_in_d = new MKL_Complex16[md.size * md.size];
	MKL_Complex16 * evals = new MKL_Complex16[md.size];

	for (int state_id_1 = 0; state_id_1 < md.size; state_id_1++)
	{
		for (int state_id_2 = 0; state_id_2 < md.size; state_id_2++)
		{
			rho_in_d[state_id_1 * md.size + state_id_2].real = 0.0;
			rho_in_d[state_id_1 * md.size + state_id_2].imag = 0.0;
		}
	}

	for (int i = 0; i < m->Rho->N; i++)
	{
		for (int k = m->Rho->RowIndex[i]; k < m->Rho->RowIndex[i + 1]; k++)
		{
			rho_in_d[i * md.size + m->Rho->Col[k]].real = m->Rho->Value[k].re;
			rho_in_d[i * md.size + m->Rho->Col[k]].imag = m->Rho->Value[k].im;
		}
	}

	for (int st_id_1 = 0; st_id_1 < md.size; st_id_1++)
	{
		for (int st_id_2 = 0; st_id_2 < md.size; st_id_2++)
		{
			pd.deep_avg_rho[st_id_1 * md.size + st_id_2].real += rho_in_d[st_id_1 * md.size + st_id_2].real;
			pd.deep_avg_rho[st_id_1 * md.size + st_id_2].imag += rho_in_d[st_id_1 * md.size + st_id_2].imag;
		}
	}

	int info = LAPACKE_zgeev(
		LAPACK_ROW_MAJOR,
		'N',
		'N',
		md.size,
		rho_in_d,
		md.size,
		evals,
		NULL,
		md.size,
		NULL,
		md.size);

	if (info > 0) 
	{
		printf("The algorithm failed to compute eigenvalues.\n");
		exit(1);
	}

	for (int st_id = 0; st_id < md.size; st_id++)
	{
		pd.deep_evals[dump_id * md.size + st_id] = evals[st_id].real;
	}

	delete[] evals;
	delete[] rho_in_d;
}

void f_basis_init(Model* model, RunParam &rp, ConfigParam &cp, MainData &md)
{
	string fn;

	double time = omp_get_wtime();
	double init_time = time;

	time = omp_get_wtime() - init_time;
	cout << "Time of createModel: " << time << endl << endl;

	initFs(model->Fs, model->N);
	time = omp_get_wtime() - init_time;;
	cout << "Time of initFs: " << time << endl << endl;

	init_h_0_vector(model, rp, cp, md);
	time = omp_get_wtime() - init_time;
	cout << "Time of init_h_0_vector: " << time << endl << endl;

	init_h_1_vector(model, rp, cp, md);
	time = omp_get_wtime() - init_time;
	cout << "Time of init_h_1_vector: " << time << endl << endl;

	init_H0(model, rp, cp, md); 
	init_H1(model, rp, cp, md);

	if (rp.issmtx == 1)
	{
		fn = rp.path + "H0" + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, model->H0, 16, false);

		fn = rp.path + "H1" + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, model->H1, 16, false);
	}

	if (cp.dt == 1)
	{
		init_diss_1(model, rp, cp, md);
	}
	else
	{
		init_diss_1(model, rp, cp, md);
	}

	time = omp_get_wtime() - init_time;
	cout << "time of init_diss_" << to_string(cp.dt) << ": " << time << endl << endl;

	init_f_d_valentin(model);
	time = omp_get_wtime() - init_time;
	cout << "time of init_f_d_valentin: " << time << endl << endl;

	transpFs(model);
	time = omp_get_wtime() - init_time;
	cout << "time of transpFs: " << time << endl << endl;

	calc_Q_0(model);
	time = omp_get_wtime() - init_time;
	cout << "time of calc_Q_0: " << time << endl << endl;
	if (rp.debug == 1)
	{
		fn = rp.path + "Q_0" + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, model->Q_0, 16, false);
	}

	calc_Q_1(model);
	time = omp_get_wtime() - init_time;
	cout << "time of calc_Q_1: " << time << endl << endl;
	if (rp.debug == 1)
	{
		fn = rp.path + "Q_1" + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, model->Q_1, 16, false);
	}

	calcKs(model);
	time = omp_get_wtime() - init_time;
	cout << "time of calcKs: " << time << endl << endl;

	if (rp.debug == 1)
	{
		fn = rp.path + "Ks" + file_name_suffix(cp, 4);
		save_complex_data(fn, (MKL_Complex16 *)model->Ks, model->N_mat, 15, false);
	}

	calcRs(model);
	time = omp_get_wtime() - init_time;
	cout << "time of calcRs: " << time << endl << endl;
	if (rp.debug == 1)
	{
		fn = rp.path + "Rs" + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, model->Rs, 16, false);
	}

	calc_G_0_s(model, cp);
	time = omp_get_wtime() - init_time;
	cout << "time of calc_G_0_s: " << time << endl << endl;

	if (rp.debug == 1)
	{
		fn = rp.path + "G_0_s" + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, model->G_0_s, 16, false);
	}

	calc_G_1_s(model, cp);
	time = omp_get_wtime() - init_time;
	cout << "time of calc_G_1_s: " << time << endl << endl;

	if (rp.debug == 1)
	{
		fn = rp.path + "G_1_s" + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, model->G_1_s, 16, false);
	}

	initRhoODE(model, rp, cp, md);
	time = omp_get_wtime() - init_time;
	cout << "time of initRhoODE: " << time << endl << endl;

	if (rp.debug == 1)
	{
		fn = rp.path + "init_rho_f" + file_name_suffix(cp, 4);
		save_complex_data(fn, (MKL_Complex16 *)model->RhoF, model->N_mat, 16, false);
	}
}

