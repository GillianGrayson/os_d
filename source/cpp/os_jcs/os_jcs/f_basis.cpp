#include "f_basis.h"
#include <string.h>
#include <map>

#define IND(i, j, k) ((ulli)i)*((ulli)N * N - 1)*((ulli)N * N - 1) + ((ulli)j)*((ulli)N * N - 1) + ((ulli)k)

#define IndS(i, j) ((int)(((N - 1) + (N - i))/2.0 * (i) + (j - i - 1)) * 2)
#define IndJ(i, j) ((int)(((N - 1) + (N - i))/2.0 * (i) + (j - i - 1)) * 2 + 1)
#define IndD(l)    (N * (N-1) + l - 1)

ulli max_bit(ulli * mas, unsigned int n)
{
	ulli mask = 0;

	for (unsigned int i = 0; i < n; i++)
	{
		mask |= mas[i];
	}

	ulli res = 1;
	mask = mask >> 1;

	while (res < mask)
	{
		res = res << 1;
	}

	return res;
}
inline void swap(Tensor_Coordinates  * mas, unsigned int i, unsigned int j)
{
	ulli tmp_ulli = mas->hash[i];
	mas->hash[i] = mas->hash[j];
	mas->hash[j] = tmp_ulli;

	unsigned int tmp_ui = mas->coord1[i];
	mas->coord1[i] = mas->coord1[j];
	mas->coord1[j] = tmp_ui;

	tmp_ui = mas->coord2[i];
	mas->coord2[i] = mas->coord2[j];
	mas->coord2[j] = tmp_ui;

	tmp_ui = mas->coord3[i];
	mas->coord3[i] = mas->coord3[j];
	mas->coord3[j] = tmp_ui;

	dcomplex tmp_mklc = mas->data[i];
	mas->data[i] = mas->data[j];
	mas->data[j] = tmp_mklc;
}
void msd_sort(Tensor_Coordinates * mas, unsigned int from, unsigned int to, ulli bit, int threads_level)
{
	if (!bit || to < from + 1) return;

	unsigned int left = from, right = to - 1;

	while (true) {

		while (left < right && !(mas->hash[left] & bit)) left++;

		while (left < right && (mas->hash[right] & bit)) right--;

		if (left >= right)
			break;
		else
			swap(mas, left, right);
	}

	if (!(bit & mas->hash[left]) && left < to) left++;

	bit >>= 1;
	if (threads_level == 0)
	{
		msd_sort(mas, from, left, bit, threads_level);
		msd_sort(mas, left, to, bit, threads_level);
	}
	else
	{
		threads_level = threads_level - 1;
#pragma omp parallel
#pragma omp single nowait
		{
#pragma omp task
			{
				msd_sort(mas, from, left, bit, threads_level);
			}
#pragma omp task
			{
				msd_sort(mas, left, to, bit, threads_level);
			}
		}
	}
}
void sort_matrix(Tensor_Coordinates * matrix)
{
	unsigned int threads_level = 0;
	ulli bit = max_bit(matrix->hash, matrix->k);

	msd_sort(matrix, 0, matrix->k, bit, threads_level);
}

Tensor_Coordinates * create_matrix(int NZ)
{
	Tensor_Coordinates * matrix = new Tensor_Coordinates;
	matrix->data = new dcomplex[NZ];
	matrix->coord1 = new unsigned int[NZ];
	matrix->coord2 = new unsigned int[NZ];
	matrix->coord3 = new unsigned int[NZ];
	matrix->hash = new unsigned long long int[NZ];
	matrix->k = NZ;
	return matrix;
}
void free_matrix(Tensor_Coordinates * matrix)
{
	delete[](matrix->data);
	delete[](matrix->coord1);
	delete[](matrix->coord2);
	delete[](matrix->coord3);
	delete[](matrix->hash);
	matrix->k = 0;
}

void fijk_coord(Tensor_Coordinates * f_ijk, int N)
{
	unsigned int size = 5 * N * N * N - 9 * N * N - 2 * N + 6;
	double tmp = 0.0;
	f_ijk->data = new dcomplex[size];
	memset(f_ijk->data, 0, size * sizeof(dcomplex));
	f_ijk->coord1 = new unsigned int[size];
	memset(f_ijk->coord1, 0, size * sizeof(unsigned int));
	f_ijk->coord2 = new unsigned int[size];
	memset(f_ijk->coord2, 0, size * sizeof(unsigned int));
	f_ijk->coord3 = new unsigned int[size];
	memset(f_ijk->coord3, 0, size * sizeof(unsigned int));
	f_ijk->hash = new ulli[size];
	memset(f_ijk->hash, 0, size * sizeof(unsigned int));
	f_ijk->k = 0;

	for (int i = 1; i < N; i++)
	{
		for (int j = i + 1; j < N; j++)
		{
			f_ijk->data[f_ijk->k].re = (i) / sqrt((double)(i) * (i + 1));
			f_ijk->coord1[f_ijk->k] = IndD(i);
			f_ijk->coord2[f_ijk->k] = IndJ(i, j);
			f_ijk->coord3[f_ijk->k] = IndS(i, j);
			f_ijk->hash[f_ijk->k] = IND(IndD(i), IndJ(i, j), IndS(i, j));
			f_ijk->k++;

			f_ijk->data[f_ijk->k].re = -(i) / sqrt((double)(i) * (i + 1));
			f_ijk->coord1[f_ijk->k] = IndD(i);
			f_ijk->coord2[f_ijk->k] = IndS(i, j);
			f_ijk->coord3[f_ijk->k] = IndJ(i, j);
			f_ijk->hash[f_ijk->k] = IND(IndD(i), IndS(i, j), IndJ(i, j));
			f_ijk->k++;

			f_ijk->data[f_ijk->k].re = -(i) / sqrt((double)(i) * (i + 1));
			f_ijk->coord1[f_ijk->k] = IndJ(i, j);
			f_ijk->coord2[f_ijk->k] = IndD(i);
			f_ijk->coord3[f_ijk->k] = IndS(i, j);
			f_ijk->hash[f_ijk->k] = IND(IndJ(i, j), IndD(i), IndS(i, j));
			f_ijk->k++;

			f_ijk->data[f_ijk->k].re = (i) / sqrt((double)(i) * (i + 1));
			f_ijk->coord1[f_ijk->k] = IndS(i, j);
			f_ijk->coord2[f_ijk->k] = IndD(i);
			f_ijk->coord3[f_ijk->k] = IndJ(i, j);
			f_ijk->hash[f_ijk->k] = IND(IndS(i, j), IndD(i), IndJ(i, j));
			f_ijk->k++;

			f_ijk->data[f_ijk->k].re = (i) / sqrt((double)(i) * (i + 1));
			f_ijk->coord1[f_ijk->k] = IndJ(i, j);
			f_ijk->coord2[f_ijk->k] = IndS(i, j);
			f_ijk->coord3[f_ijk->k] = IndD(i);
			f_ijk->hash[f_ijk->k] = IND(IndJ(i, j), IndS(i, j), IndD(i));
			f_ijk->k++;

			f_ijk->data[f_ijk->k].re = -(i) / sqrt((double)(i) * (i + 1));
			f_ijk->coord1[f_ijk->k] = IndS(i, j);
			f_ijk->coord2[f_ijk->k] = IndJ(i, j);
			f_ijk->coord3[f_ijk->k] = IndD(i);
			f_ijk->hash[f_ijk->k] = IND(IndS(i, j), IndJ(i, j), IndD(i));
			f_ijk->k++;
		}
	}

	for (int i = 0; i < N; i++)
	{
		for (int j = i + 1; j < N; j++)
		{
			for (int m = i + 1; m < j; m++)
			{
				f_ijk->data[f_ijk->k].re = -1.0 / sqrt((double)(m) * (m + 1));
				f_ijk->coord1[f_ijk->k] = IndD(m);
				f_ijk->coord2[f_ijk->k] = IndJ(i, j);
				f_ijk->coord3[f_ijk->k] = IndS(i, j);
				f_ijk->hash[f_ijk->k] = IND(IndD(m), IndJ(i, j), IndS(i, j));
				f_ijk->k++;

				f_ijk->data[f_ijk->k].re = 1.0 / sqrt((double)(m) * (m + 1));
				f_ijk->coord1[f_ijk->k] = IndD(m);
				f_ijk->coord2[f_ijk->k] = IndS(i, j);
				f_ijk->coord3[f_ijk->k] = IndJ(i, j);
				f_ijk->hash[f_ijk->k] = IND(IndD(m), IndS(i, j), IndJ(i, j));
				f_ijk->k++;

				f_ijk->data[f_ijk->k].re = 1.0 / sqrt((double)(m) * (m + 1));
				f_ijk->coord1[f_ijk->k] = IndJ(i, j);
				f_ijk->coord2[f_ijk->k] = IndD(m);
				f_ijk->coord3[f_ijk->k] = IndS(i, j);
				f_ijk->hash[f_ijk->k] = IND(IndJ(i, j), IndD(m), IndS(i, j));
				f_ijk->k++;

				f_ijk->data[f_ijk->k].re = -1.0 / sqrt((double)(m) * (m + 1));
				f_ijk->coord1[f_ijk->k] = IndS(i, j);
				f_ijk->coord2[f_ijk->k] = IndD(m);
				f_ijk->coord3[f_ijk->k] = IndJ(i, j);
				f_ijk->hash[f_ijk->k] = IND(IndS(i, j), IndD(m), IndJ(i, j));
				f_ijk->k++;

				f_ijk->data[f_ijk->k].re = -1.0 / sqrt((double)(m) * (m + 1));
				f_ijk->coord1[f_ijk->k] = IndJ(i, j);
				f_ijk->coord2[f_ijk->k] = IndS(i, j);
				f_ijk->coord3[f_ijk->k] = IndD(m);
				f_ijk->hash[f_ijk->k] = IND(IndJ(i, j), IndS(i, j), IndD(m));
				f_ijk->k++;

				f_ijk->data[f_ijk->k].re = 1.0 / sqrt((double)(m) * (m + 1));
				f_ijk->coord1[f_ijk->k] = IndS(i, j);
				f_ijk->coord2[f_ijk->k] = IndJ(i, j);
				f_ijk->coord3[f_ijk->k] = IndD(m);
				f_ijk->hash[f_ijk->k] = IND(IndS(i, j), IndJ(i, j), IndD(m));
				f_ijk->k++;
			}
		}
	}

	for (int i = 0; i < N; i++)
	{
		for (int j = i + 1; j < N; j++)
		{
			f_ijk->data[f_ijk->k].re = -(1 + j) / sqrt((double)(j) * (j + 1));
			f_ijk->coord1[f_ijk->k] = IndD(j);
			f_ijk->coord2[f_ijk->k] = IndJ(i, j);
			f_ijk->coord3[f_ijk->k] = IndS(i, j);
			f_ijk->hash[f_ijk->k] = IND(IndD(j), IndJ(i, j), IndS(i, j));
			f_ijk->k++;

			f_ijk->data[f_ijk->k].re = (1 + j) / sqrt((double)(j) * (j + 1));
			f_ijk->coord1[f_ijk->k] = IndD(j);
			f_ijk->coord2[f_ijk->k] = IndS(i, j);
			f_ijk->coord3[f_ijk->k] = IndJ(i, j);
			f_ijk->hash[f_ijk->k] = IND(IndD(j), IndS(i, j), IndJ(i, j));
			f_ijk->k++;

			f_ijk->data[f_ijk->k].re = (1 + j) / sqrt((double)(j) * (j + 1));
			f_ijk->coord1[f_ijk->k] = IndJ(i, j);
			f_ijk->coord2[f_ijk->k] = IndD(j);
			f_ijk->coord3[f_ijk->k] = IndS(i, j);
			f_ijk->hash[f_ijk->k] = IND(IndJ(i, j), IndD(j), IndS(i, j));
			f_ijk->k++;

			f_ijk->data[f_ijk->k].re = -(1 + j) / sqrt((double)(j) * (j + 1));
			f_ijk->coord1[f_ijk->k] = IndS(i, j);
			f_ijk->coord2[f_ijk->k] = IndD(j);
			f_ijk->coord3[f_ijk->k] = IndJ(i, j);
			f_ijk->hash[f_ijk->k] = IND(IndS(i, j), IndD(j), IndJ(i, j));
			f_ijk->k++;

			f_ijk->data[f_ijk->k].re = -(1 + j) / sqrt((double)(j) * (j + 1));
			f_ijk->coord1[f_ijk->k] = IndJ(i, j);
			f_ijk->coord2[f_ijk->k] = IndS(i, j);
			f_ijk->coord3[f_ijk->k] = IndD(j);
			f_ijk->hash[f_ijk->k] = IND(IndJ(i, j), IndS(i, j), IndD(j));
			f_ijk->k++;

			f_ijk->data[f_ijk->k].re = (1 + j) / sqrt((double)(j) * (j + 1));
			f_ijk->coord1[f_ijk->k] = IndS(i, j);
			f_ijk->coord2[f_ijk->k] = IndJ(i, j);
			f_ijk->coord3[f_ijk->k] = IndD(j);
			f_ijk->hash[f_ijk->k] = IND(IndS(i, j), IndJ(i, j), IndD(j));
			f_ijk->k++;
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
						f_ijk->data[f_ijk->k].re = 1.0 / sqrt(2.0); //-i*Jjk*[Sij,Sik]  //and symmetric Ski, Jkj
						f_ijk->coord1[f_ijk->k] = IndJ(j, k);
						f_ijk->coord2[f_ijk->k] = IndS(i, j);
						f_ijk->coord3[f_ijk->k] = IndS(i, k);
						f_ijk->hash[f_ijk->k] = IND(IndJ(j, k), IndS(i, j), IndS(i, k));
						f_ijk->k++;

						f_ijk->data[f_ijk->k].re = 1.0 / sqrt(2.0); //-i*Jik*[Sij,Sjk]  //and symmetric Skj, Jki
						f_ijk->coord1[f_ijk->k] = IndJ(i, k);
						f_ijk->coord2[f_ijk->k] = IndS(i, j);
						f_ijk->coord3[f_ijk->k] = IndS(j, k);
						f_ijk->hash[f_ijk->k] = IND(IndJ(i, k), IndS(i, j), IndS(j, k));
						f_ijk->k++;

						f_ijk->data[f_ijk->k].re = -1.0 / sqrt(2.0); //-i*Sjk*[Sij,Jik]  //and symmetric Skj, Jki
						f_ijk->coord1[f_ijk->k] = IndS(j, k);
						f_ijk->coord2[f_ijk->k] = IndS(i, j);
						f_ijk->coord3[f_ijk->k] = IndJ(i, k);
						f_ijk->hash[f_ijk->k] = IND(IndS(j, k), IndS(i, j), IndJ(i, k));
						f_ijk->k++;

						f_ijk->data[f_ijk->k].re = -1.0 / sqrt(2.0); //-i*Sik*[Sij,Jjk]  //and symmetric Ski, Jkj
						f_ijk->coord1[f_ijk->k] = IndS(i, k);
						f_ijk->coord2[f_ijk->k] = IndS(i, j);
						f_ijk->coord3[f_ijk->k] = IndJ(j, k);
						f_ijk->hash[f_ijk->k] = IND(IndS(i, k), IndS(i, j), IndJ(j, k));
						f_ijk->k++;

						f_ijk->data[f_ijk->k].re = 1.0 / sqrt(2.0); //-i*Sjk*[Jij,Sik]  //and symmetric Ski, Skj
						f_ijk->coord1[f_ijk->k] = IndS(j, k);
						f_ijk->coord2[f_ijk->k] = IndJ(i, j);
						f_ijk->coord3[f_ijk->k] = IndS(i, k);
						f_ijk->hash[f_ijk->k] = IND(IndS(j, k), IndJ(i, j), IndS(i, k));
						f_ijk->k++;

						f_ijk->data[f_ijk->k].re = -1.0 / sqrt(2.0); //-i*Sik*[Jij,Sjk]  //and symmetric Skj, Ski
						f_ijk->coord1[f_ijk->k] = IndS(i, k);
						f_ijk->coord2[f_ijk->k] = IndJ(i, j);
						f_ijk->coord3[f_ijk->k] = IndS(j, k);
						f_ijk->hash[f_ijk->k] = IND(IndS(i, k), IndJ(i, j), IndS(j, k));
						f_ijk->k++;

						f_ijk->data[f_ijk->k].re = 1.0 / sqrt(2.0); //-i*Jjk*[Jij,Jik]  //and symmetric Jkj, Jki
						f_ijk->coord1[f_ijk->k] = IndJ(j, k);
						f_ijk->coord2[f_ijk->k] = IndJ(i, j);
						f_ijk->coord3[f_ijk->k] = IndJ(i, k);
						f_ijk->hash[f_ijk->k] = IND(IndJ(j, k), IndJ(i, j), IndJ(i, k));
						f_ijk->k++;

						f_ijk->data[f_ijk->k].re = -1.0 / sqrt(2.0); //-i*Jik*[Jij,Jjk]  //and symmetric Jki, Jkj
						f_ijk->coord1[f_ijk->k] = IndJ(i, k);
						f_ijk->coord2[f_ijk->k] = IndJ(i, j);
						f_ijk->coord3[f_ijk->k] = IndJ(j, k);
						f_ijk->hash[f_ijk->k] = IND(IndJ(i, k), IndJ(i, j), IndJ(j, k));
						f_ijk->k++;
					}
				}
				if (k < j)
				{
					if (k > i)
					{
						f_ijk->data[f_ijk->k].re = -1.0 / sqrt(2.0);
						f_ijk->coord1[f_ijk->k] = IndJ(k, j);
						f_ijk->coord2[f_ijk->k] = IndS(i, j);
						f_ijk->coord3[f_ijk->k] = IndS(i, k);
						f_ijk->hash[f_ijk->k] = IND(IndJ(k, j), IndS(i, j), IndS(i, k));
						f_ijk->k++;

						f_ijk->data[f_ijk->k].re = 1.0 / sqrt(2.0);
						f_ijk->coord1[f_ijk->k] = IndJ(i, k);
						f_ijk->coord2[f_ijk->k] = IndS(i, j);
						f_ijk->coord3[f_ijk->k] = IndS(k, j);
						f_ijk->hash[f_ijk->k] = IND(IndJ(i, k), IndS(i, j), IndS(k, j));
						f_ijk->k++;

						f_ijk->data[f_ijk->k].re = -1.0 / sqrt(2.0);
						f_ijk->coord1[f_ijk->k] = IndS(k, j);
						f_ijk->coord2[f_ijk->k] = IndS(i, j);
						f_ijk->coord3[f_ijk->k] = IndJ(i, k);
						f_ijk->hash[f_ijk->k] = IND(IndS(k, j), IndS(i, j), IndJ(i, k));
						f_ijk->k++;

						f_ijk->data[f_ijk->k].re = 1.0 / sqrt(2.0);
						f_ijk->coord1[f_ijk->k] = IndS(i, k);
						f_ijk->coord2[f_ijk->k] = IndS(i, j);
						f_ijk->coord3[f_ijk->k] = IndJ(k, j);
						f_ijk->hash[f_ijk->k] = IND(IndS(i, k), IndS(i, j), IndJ(k, j));
						f_ijk->k++;

						f_ijk->data[f_ijk->k].re = 1.0 / sqrt(2.0);
						f_ijk->coord1[f_ijk->k] = IndS(k, j);
						f_ijk->coord2[f_ijk->k] = IndJ(i, j);
						f_ijk->coord3[f_ijk->k] = IndS(i, k);
						f_ijk->hash[f_ijk->k] = IND(IndS(k, j), IndJ(i, j), IndS(i, k));
						f_ijk->k++;

						f_ijk->data[f_ijk->k].re = -1.0 / sqrt(2.0);
						f_ijk->coord1[f_ijk->k] = IndS(i, k);
						f_ijk->coord2[f_ijk->k] = IndJ(i, j);
						f_ijk->coord3[f_ijk->k] = IndS(k, j);
						f_ijk->hash[f_ijk->k] = IND(IndS(i, k), IndJ(i, j), IndS(k, j));
						f_ijk->k++;

						f_ijk->data[f_ijk->k].re = -1.0 / sqrt(2.0);
						f_ijk->coord1[f_ijk->k] = IndJ(k, j);
						f_ijk->coord2[f_ijk->k] = IndJ(i, j);
						f_ijk->coord3[f_ijk->k] = IndJ(i, k);
						f_ijk->hash[f_ijk->k] = IND(IndJ(k, j), IndJ(i, j), IndJ(i, k));
						f_ijk->k++;

						f_ijk->data[f_ijk->k].re = 1.0 / sqrt(2.0);
						f_ijk->coord1[f_ijk->k] = IndJ(i, k);
						f_ijk->coord2[f_ijk->k] = IndJ(i, j);
						f_ijk->coord3[f_ijk->k] = IndJ(k, j);
						f_ijk->hash[f_ijk->k] = IND(IndJ(i, k), IndJ(i, j), IndJ(k, j));
						f_ijk->k++;
					}
					if (k < i)
					{
						f_ijk->data[f_ijk->k].re = -1.0 / sqrt(2.0);
						f_ijk->coord1[f_ijk->k] = IndJ(k, j);
						f_ijk->coord2[f_ijk->k] = IndS(i, j);
						f_ijk->coord3[f_ijk->k] = IndS(k, i);
						f_ijk->hash[f_ijk->k] = IND(IndJ(k, j), IndS(i, j), IndS(k, i));
						f_ijk->k++;

						f_ijk->data[f_ijk->k].re = -1.0 / sqrt(2.0);
						f_ijk->coord1[f_ijk->k] = IndJ(k, i);
						f_ijk->coord2[f_ijk->k] = IndS(i, j);
						f_ijk->coord3[f_ijk->k] = IndS(k, j);
						f_ijk->hash[f_ijk->k] = IND(IndJ(k, i), IndS(i, j), IndS(k, j));
						f_ijk->k++;

						f_ijk->data[f_ijk->k].re = 1.0 / sqrt(2.0);
						f_ijk->coord1[f_ijk->k] = IndS(k, j);
						f_ijk->coord2[f_ijk->k] = IndS(i, j);
						f_ijk->coord3[f_ijk->k] = IndJ(k, i);
						f_ijk->hash[f_ijk->k] = IND(IndS(k, j), IndS(i, j), IndJ(k, i));
						f_ijk->k++;

						f_ijk->data[f_ijk->k].re = 1.0 / sqrt(2.0);
						f_ijk->coord1[f_ijk->k] = IndS(k, i);
						f_ijk->coord2[f_ijk->k] = IndS(i, j);
						f_ijk->coord3[f_ijk->k] = IndJ(k, j);
						f_ijk->hash[f_ijk->k] = IND(IndS(k, i), IndS(i, j), IndJ(k, j));
						f_ijk->k++;

						f_ijk->data[f_ijk->k].re = 1.0 / sqrt(2.0);
						f_ijk->coord1[f_ijk->k] = IndS(k, j);
						f_ijk->coord2[f_ijk->k] = IndJ(i, j);
						f_ijk->coord3[f_ijk->k] = IndS(k, i);
						f_ijk->hash[f_ijk->k] = IND(IndS(k, j), IndJ(i, j), IndS(k, i));
						f_ijk->k++;

						f_ijk->data[f_ijk->k].re = -1.0 / sqrt(2.0);
						f_ijk->coord1[f_ijk->k] = IndS(k, i);
						f_ijk->coord2[f_ijk->k] = IndJ(i, j);
						f_ijk->coord3[f_ijk->k] = IndS(k, j);
						f_ijk->hash[f_ijk->k] = IND(IndS(k, i), IndJ(i, j), IndS(k, j));
						f_ijk->k++;

						f_ijk->data[f_ijk->k].re = 1.0 / sqrt(2.0);
						f_ijk->coord1[f_ijk->k] = IndJ(k, j);
						f_ijk->coord2[f_ijk->k] = IndJ(i, j);
						f_ijk->coord3[f_ijk->k] = IndJ(k, i);
						f_ijk->hash[f_ijk->k] = IND(IndJ(k, j), IndJ(i, j), IndJ(k, i));
						f_ijk->k++;

						f_ijk->data[f_ijk->k].re = -1.0 / sqrt(2.0);
						f_ijk->coord1[f_ijk->k] = IndJ(k, i);
						f_ijk->coord2[f_ijk->k] = IndJ(i, j);
						f_ijk->coord3[f_ijk->k] = IndJ(k, j);
						f_ijk->hash[f_ijk->k] = IND(IndJ(k, i), IndJ(i, j), IndJ(k, j));
						f_ijk->k++;
					}
				}
			}
		}
	}
}
void dijk_coord(Tensor_Coordinates * d_ijk, int N)
{
	unsigned int size = 6 * N * N * N - (N * (21 * N + 7)) / 2 + 1;
	d_ijk->data = new dcomplex[size];
	memset(d_ijk->data, 0, size * sizeof(dcomplex));
	d_ijk->coord1 = new unsigned int[size];
	memset(d_ijk->coord1, 0, size * sizeof(unsigned int));
	d_ijk->coord2 = new unsigned int[size];
	memset(d_ijk->coord2, 0, size * sizeof(unsigned int));
	d_ijk->coord3 = new unsigned int[size];
	memset(d_ijk->coord3, 0, size * sizeof(unsigned int));
	d_ijk->hash = new ulli[size];
	memset(d_ijk->hash, 0, size * sizeof(unsigned int));
	d_ijk->k = 0;

	for (int i = 1; i < N; i++)
	{
		for (int j = i + 1; j < N; j++)
		{
			d_ijk->data[d_ijk->k].re = -(i) / sqrt((double)(i) * (i + 1));
			d_ijk->coord1[d_ijk->k] = IndD(i);
			d_ijk->coord2[d_ijk->k] = IndS(i, j);
			d_ijk->coord3[d_ijk->k] = IndS(i, j);
			d_ijk->hash[d_ijk->k] = IND(IndD(i), IndS(i, j), IndS(i, j));
			d_ijk->k++;

			d_ijk->data[d_ijk->k].re = -(i) / sqrt((double)(i) * (i + 1));
			d_ijk->coord1[d_ijk->k] = IndS(i, j);
			d_ijk->coord2[d_ijk->k] = IndD(i);
			d_ijk->coord3[d_ijk->k] = IndS(i, j);
			d_ijk->hash[d_ijk->k] = IND(IndS(i, j), IndD(i), IndS(i, j));
			d_ijk->k++;

			d_ijk->data[d_ijk->k].re = -(i) / sqrt((double)(i) * (i + 1));
			d_ijk->coord1[d_ijk->k] = IndD(i);
			d_ijk->coord2[d_ijk->k] = IndJ(i, j);
			d_ijk->coord3[d_ijk->k] = IndJ(i, j);
			d_ijk->hash[d_ijk->k] = IND(IndD(i), IndJ(i, j), IndJ(i, j));
			d_ijk->k++;

			d_ijk->data[d_ijk->k].re = -(i) / sqrt((double)(i) * (i + 1));
			d_ijk->coord1[d_ijk->k] = IndJ(i, j);
			d_ijk->coord2[d_ijk->k] = IndD(i);
			d_ijk->coord3[d_ijk->k] = IndJ(i, j);
			d_ijk->hash[d_ijk->k] = IND(IndJ(i, j), IndD(i), IndJ(i, j));
			d_ijk->k++;

			d_ijk->data[d_ijk->k].re = -(i) / sqrt((double)(i) * (i + 1));
			d_ijk->coord1[d_ijk->k] = IndS(i, j);
			d_ijk->coord2[d_ijk->k] = IndS(i, j);
			d_ijk->coord3[d_ijk->k] = IndD(i);
			d_ijk->hash[d_ijk->k] = IND(IndS(i, j), IndS(i, j), IndD(i));
			d_ijk->k++;

			d_ijk->data[d_ijk->k].re = -(i) / sqrt((double)(i) * (i + 1));
			d_ijk->coord1[d_ijk->k] = IndJ(i, j);
			d_ijk->coord2[d_ijk->k] = IndJ(i, j);
			d_ijk->coord3[d_ijk->k] = IndD(i);
			d_ijk->hash[d_ijk->k] = IND(IndJ(i, j), IndJ(i, j), IndD(i));
			d_ijk->k++;
		}
	}

	for (int i = 0; i < N; i++)
	{
		for (int j = i + 1; j < N; j++)
		{
			for (int m = i + 1; m < j; m++)
			{
				d_ijk->data[d_ijk->k].re = 1.0 / sqrt((double)(m) * (m + 1));
				d_ijk->coord1[d_ijk->k] = IndD(m);
				d_ijk->coord2[d_ijk->k] = IndS(i, j);
				d_ijk->coord3[d_ijk->k] = IndS(i, j);
				d_ijk->hash[d_ijk->k] = IND(IndD(m), IndS(i, j), IndS(i, j));
				d_ijk->k++;

				d_ijk->data[d_ijk->k].re = 1.0 / sqrt((double)(m) * (m + 1));
				d_ijk->coord1[d_ijk->k] = IndS(i, j);
				d_ijk->coord2[d_ijk->k] = IndD(m);
				d_ijk->coord3[d_ijk->k] = IndS(i, j);
				d_ijk->hash[d_ijk->k] = IND(IndS(i, j), IndD(m), IndS(i, j));
				d_ijk->k++;

				d_ijk->data[d_ijk->k].re = 1.0 / sqrt((double)(m) * (m + 1));
				d_ijk->coord1[d_ijk->k] = IndD(m);
				d_ijk->coord2[d_ijk->k] = IndJ(i, j);
				d_ijk->coord3[d_ijk->k] = IndJ(i, j);
				d_ijk->hash[d_ijk->k] = IND(IndD(m), IndJ(i, j), IndJ(i, j));
				d_ijk->k++;

				d_ijk->data[d_ijk->k].re = 1.0 / sqrt((double)(m) * (m + 1));
				d_ijk->coord1[d_ijk->k] = IndJ(i, j);
				d_ijk->coord2[d_ijk->k] = IndD(m);
				d_ijk->coord3[d_ijk->k] = IndJ(i, j);
				d_ijk->hash[d_ijk->k] = IND(IndJ(i, j), IndD(m), IndJ(i, j));
				d_ijk->k++;

				d_ijk->data[d_ijk->k].re = 1.0 / sqrt((double)(m) * (m + 1));
				d_ijk->coord1[d_ijk->k] = IndS(i, j);
				d_ijk->coord2[d_ijk->k] = IndS(i, j);
				d_ijk->coord3[d_ijk->k] = IndD(m);
				d_ijk->hash[d_ijk->k] = IND(IndS(i, j), IndS(i, j), IndD(m));
				d_ijk->k++;

				d_ijk->data[d_ijk->k].re = 1.0 / sqrt((double)(m) * (m + 1));
				d_ijk->coord1[d_ijk->k] = IndJ(i, j);
				d_ijk->coord2[d_ijk->k] = IndJ(i, j);
				d_ijk->coord3[d_ijk->k] = IndD(m);
				d_ijk->hash[d_ijk->k] = IND(IndJ(i, j), IndJ(i, j), IndD(m));
				d_ijk->k++;
			}
		}
	}

	for (int j = 2; j < N; j++)
	{
		int i = 0;
		d_ijk->data[d_ijk->k].re = (1 - j) / sqrt((double)(j) * (j + 1));
		d_ijk->coord1[d_ijk->k] = IndD(j);
		d_ijk->coord2[d_ijk->k] = IndS(i, j);
		d_ijk->coord3[d_ijk->k] = IndS(i, j);
		d_ijk->hash[d_ijk->k] = IND(IndD(j), IndS(i, j), IndS(i, j));
		d_ijk->k++;

		d_ijk->data[d_ijk->k].re = (1 - j) / sqrt((double)(j) * (j + 1));
		d_ijk->coord1[d_ijk->k] = IndS(i, j);
		d_ijk->coord2[d_ijk->k] = IndD(j);
		d_ijk->coord3[d_ijk->k] = IndS(i, j);
		d_ijk->hash[d_ijk->k] = IND(IndS(i, j), IndD(j), IndS(i, j));
		d_ijk->k++;

		d_ijk->data[d_ijk->k].re = (1 - j) / sqrt((double)(j) * (j + 1));
		d_ijk->coord1[d_ijk->k] = IndD(j);
		d_ijk->coord2[d_ijk->k] = IndJ(i, j);
		d_ijk->coord3[d_ijk->k] = IndJ(i, j);
		d_ijk->hash[d_ijk->k] = IND(IndD(j), IndJ(i, j), IndJ(i, j));
		d_ijk->k++;

		d_ijk->data[d_ijk->k].re = (1 - j) / sqrt((double)(j) * (j + 1));
		d_ijk->coord1[d_ijk->k] = IndJ(i, j);
		d_ijk->coord2[d_ijk->k] = IndD(j);
		d_ijk->coord3[d_ijk->k] = IndJ(i, j);
		d_ijk->hash[d_ijk->k] = IND(IndJ(i, j), IndD(j), IndJ(i, j));
		d_ijk->k++;

		d_ijk->data[d_ijk->k].re = (1 - j) / sqrt((double)(j) * (j + 1));
		d_ijk->coord1[d_ijk->k] = IndS(i, j);
		d_ijk->coord2[d_ijk->k] = IndS(i, j);
		d_ijk->coord3[d_ijk->k] = IndD(j);
		d_ijk->hash[d_ijk->k] = IND(IndS(i, j), IndS(i, j), IndD(j));
		d_ijk->k++;

		d_ijk->data[d_ijk->k].re = (1 - j) / sqrt((double)(j) * (j + 1));
		d_ijk->coord1[d_ijk->k] = IndJ(i, j);
		d_ijk->coord2[d_ijk->k] = IndJ(i, j);
		d_ijk->coord3[d_ijk->k] = IndD(j);
		d_ijk->hash[d_ijk->k] = IND(IndJ(i, j), IndJ(i, j), IndD(j));
		d_ijk->k++;
	}

	for (int i = 1; i < N; i++)
	{
		for (int j = i + 1; j < N; j++)
		{
			d_ijk->data[d_ijk->k].re = (1 - j) / sqrt((double)(j) * (j + 1));
			d_ijk->coord1[d_ijk->k] = IndD(j);
			d_ijk->coord2[d_ijk->k] = IndS(i, j);
			d_ijk->coord3[d_ijk->k] = IndS(i, j);
			d_ijk->hash[d_ijk->k] = IND(IndD(j), IndS(i, j), IndS(i, j));
			d_ijk->k++;

			d_ijk->data[d_ijk->k].re = (1 - j) / sqrt((double)(j) * (j + 1));
			d_ijk->coord1[d_ijk->k] = IndS(i, j);
			d_ijk->coord2[d_ijk->k] = IndD(j);
			d_ijk->coord3[d_ijk->k] = IndS(i, j);
			d_ijk->hash[d_ijk->k] = IND(IndS(i, j), IndD(j), IndS(i, j));
			d_ijk->k++;

			d_ijk->data[d_ijk->k].re = (1 - j) / sqrt((double)(j) * (j + 1));
			d_ijk->coord1[d_ijk->k] = IndD(j);
			d_ijk->coord2[d_ijk->k] = IndJ(i, j);
			d_ijk->coord3[d_ijk->k] = IndJ(i, j);
			d_ijk->hash[d_ijk->k] = IND(IndD(j), IndJ(i, j), IndJ(i, j));
			d_ijk->k++;

			d_ijk->data[d_ijk->k].re = (1 - j) / sqrt((double)(j) * (j + 1));
			d_ijk->coord1[d_ijk->k] = IndJ(i, j);
			d_ijk->coord2[d_ijk->k] = IndD(j);
			d_ijk->coord3[d_ijk->k] = IndJ(i, j);
			d_ijk->hash[d_ijk->k] = IND(IndJ(i, j), IndD(j), IndJ(i, j));
			d_ijk->k++;

			d_ijk->data[d_ijk->k].re = (1 - j) / sqrt((double)(j) * (j + 1));
			d_ijk->coord1[d_ijk->k] = IndS(i, j);
			d_ijk->coord2[d_ijk->k] = IndS(i, j);
			d_ijk->coord3[d_ijk->k] = IndD(j);
			d_ijk->hash[d_ijk->k] = IND(IndS(i, j), IndS(i, j), IndD(j));
			d_ijk->k++;

			d_ijk->data[d_ijk->k].re = (1 - j) / sqrt((double)(j) * (j + 1));
			d_ijk->coord1[d_ijk->k] = IndJ(i, j);
			d_ijk->coord2[d_ijk->k] = IndJ(i, j);
			d_ijk->coord3[d_ijk->k] = IndD(j);
			d_ijk->hash[d_ijk->k] = IND(IndJ(i, j), IndJ(i, j), IndD(j));
			d_ijk->k++;
		}
	}

	for (int i = 0; i < N; i++)
	{
		for (int j = i + 1; j < N; j++)
		{
			for (int z = j + 1; z < N; z++)
			{
				d_ijk->data[d_ijk->k].re = 2.0 / sqrt((double)(z) * (z + 1));
				d_ijk->coord1[d_ijk->k] = IndD(z);
				d_ijk->coord2[d_ijk->k] = IndS(i, j);
				d_ijk->coord3[d_ijk->k] = IndS(i, j);
				d_ijk->hash[d_ijk->k] = IND(IndD(z), IndS(i, j), IndS(i, j));
				d_ijk->k++;

				d_ijk->data[d_ijk->k].re = 2.0 / sqrt((double)(z) * (z + 1));
				d_ijk->coord1[d_ijk->k] = IndS(i, j);
				d_ijk->coord2[d_ijk->k] = IndD(z);
				d_ijk->coord3[d_ijk->k] = IndS(i, j);
				d_ijk->hash[d_ijk->k] = IND(IndS(i, j), IndD(z), IndS(i, j));
				d_ijk->k++;

				d_ijk->data[d_ijk->k].re = 2.0 / sqrt((double)(z) * (z + 1));
				d_ijk->coord1[d_ijk->k] = IndD(z);
				d_ijk->coord2[d_ijk->k] = IndJ(i, j);
				d_ijk->coord3[d_ijk->k] = IndJ(i, j);
				d_ijk->hash[d_ijk->k] = IND(IndD(z), IndJ(i, j), IndJ(i, j));
				d_ijk->k++;

				d_ijk->data[d_ijk->k].re = 2.0 / sqrt((double)(z) * (z + 1));
				d_ijk->coord1[d_ijk->k] = IndJ(i, j);
				d_ijk->coord2[d_ijk->k] = IndD(z);
				d_ijk->coord3[d_ijk->k] = IndJ(i, j);
				d_ijk->hash[d_ijk->k] = IND(IndJ(i, j), IndD(z), IndJ(i, j));
				d_ijk->k++;

				d_ijk->data[d_ijk->k].re = 2.0 / sqrt((double)(z) * (z + 1));
				d_ijk->coord1[d_ijk->k] = IndS(i, j);
				d_ijk->coord2[d_ijk->k] = IndS(i, j);
				d_ijk->coord3[d_ijk->k] = IndD(z);
				d_ijk->hash[d_ijk->k] = IND(IndS(i, j), IndS(i, j), IndD(z));
				d_ijk->k++;

				d_ijk->data[d_ijk->k].re = 2.0 / sqrt((double)(z) * (z + 1));
				d_ijk->coord1[d_ijk->k] = IndJ(i, j);
				d_ijk->coord2[d_ijk->k] = IndJ(i, j);
				d_ijk->coord3[d_ijk->k] = IndD(z);
				d_ijk->hash[d_ijk->k] = IND(IndJ(i, j), IndJ(i, j), IndD(z));
				d_ijk->k++;
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
						d_ijk->data[d_ijk->k].re = 1.0 / sqrt(2.0);  //Sjk*{Sij,Sik}  //and symmetric Ski, Skj
						d_ijk->coord1[d_ijk->k] = IndS(j, k);
						d_ijk->coord2[d_ijk->k] = IndS(i, j);
						d_ijk->coord3[d_ijk->k] = IndS(i, k);
						d_ijk->hash[d_ijk->k] = IND(IndS(j, k), IndS(i, j), IndS(i, k));
						d_ijk->k++;

						d_ijk->data[d_ijk->k].re = 1.0 / sqrt(2.0);  //Sik*{Sij,Sjk}  //and symmetric Skj, Ski
						d_ijk->coord1[d_ijk->k] = IndS(i, k);
						d_ijk->coord2[d_ijk->k] = IndS(i, j);
						d_ijk->coord3[d_ijk->k] = IndS(j, k);
						d_ijk->hash[d_ijk->k] = IND(IndS(i, k), IndS(i, j), IndS(j, k));
						d_ijk->k++;

						d_ijk->data[d_ijk->k].re = 1.0 / sqrt(2.0);  //Jjk*{Sij,Jik}  //and symmetric Jkj, Jki
						d_ijk->coord1[d_ijk->k] = IndJ(j, k);
						d_ijk->coord2[d_ijk->k] = IndS(i, j);
						d_ijk->coord3[d_ijk->k] = IndJ(i, k);
						d_ijk->hash[d_ijk->k] = IND(IndJ(j, k), IndS(i, j), IndJ(i, k));
						d_ijk->k++;

						d_ijk->data[d_ijk->k].re = 1.0 / sqrt(2.0);  //Jik*{Sij,Jjk}  //and symmetric Jki, Jkj
						d_ijk->coord1[d_ijk->k] = IndJ(i, k);
						d_ijk->coord2[d_ijk->k] = IndS(i, j);
						d_ijk->coord3[d_ijk->k] = IndJ(j, k);
						d_ijk->hash[d_ijk->k] = IND(IndJ(i, k), IndS(i, j), IndJ(j, k));
						d_ijk->k++;

						d_ijk->data[d_ijk->k].re = -1.0 / sqrt(2.0);  //Jjk*{Jij,Sik}  //and symmetric Ski, Jkj
						d_ijk->coord1[d_ijk->k] = IndJ(j, k);
						d_ijk->coord2[d_ijk->k] = IndJ(i, j);
						d_ijk->coord3[d_ijk->k] = IndS(i, k);
						d_ijk->hash[d_ijk->k] = IND(IndJ(j, k), IndJ(i, j), IndS(i, k));
						d_ijk->k++;

						d_ijk->data[d_ijk->k].re = 1.0 / sqrt(2.0);  //Jik*{Jij,Sjk}  //and symmetric Skj, Jki
						d_ijk->coord1[d_ijk->k] = IndJ(i, k);
						d_ijk->coord2[d_ijk->k] = IndJ(i, j);
						d_ijk->coord3[d_ijk->k] = IndS(j, k);
						d_ijk->hash[d_ijk->k] = IND(IndJ(i, k), IndJ(i, j), IndS(j, k));
						d_ijk->k++;

						d_ijk->data[d_ijk->k].re = 1.0 / sqrt(2.0);  //Sjk*{Jij,Jik}  //and symmetric Skj, Jki
						d_ijk->coord1[d_ijk->k] = IndS(j, k);
						d_ijk->coord2[d_ijk->k] = IndJ(i, j);
						d_ijk->coord3[d_ijk->k] = IndJ(i, k);
						d_ijk->hash[d_ijk->k] = IND(IndS(j, k), IndJ(i, j), IndJ(i, k));
						d_ijk->k++;

						d_ijk->data[d_ijk->k].re = -1.0 / sqrt(2.0);  //Sik*{Jij,Jjk}  //and symmetric Ski, Jkj
						d_ijk->coord1[d_ijk->k] = IndS(i, k);
						d_ijk->coord2[d_ijk->k] = IndJ(i, j);
						d_ijk->coord3[d_ijk->k] = IndJ(j, k);
						d_ijk->hash[d_ijk->k] = IND(IndS(i, k), IndJ(i, j), IndJ(j, k));
						d_ijk->k++;
					}
				}
				if (k < j)
				{
					if (k > i)
					{
						d_ijk->data[d_ijk->k].re = 1.0 / sqrt(2.0);
						d_ijk->coord1[d_ijk->k] = IndS(k, j);
						d_ijk->coord2[d_ijk->k] = IndS(i, j);
						d_ijk->coord3[d_ijk->k] = IndS(i, k);
						d_ijk->hash[d_ijk->k] = IND(IndS(k, j), IndS(i, j), IndS(i, k));
						d_ijk->k++;

						d_ijk->data[d_ijk->k].re = 1.0 / sqrt(2.0);
						d_ijk->coord1[d_ijk->k] = IndS(i, k);
						d_ijk->coord2[d_ijk->k] = IndS(i, j);
						d_ijk->coord3[d_ijk->k] = IndS(k, j);
						d_ijk->hash[d_ijk->k] = IND(IndS(i, k), IndS(i, j), IndS(k, j));
						d_ijk->k++;

						d_ijk->data[d_ijk->k].re = -1.0 / sqrt(2.0);
						d_ijk->coord1[d_ijk->k] = IndJ(k, j);
						d_ijk->coord2[d_ijk->k] = IndS(i, j);
						d_ijk->coord3[d_ijk->k] = IndJ(i, k);
						d_ijk->hash[d_ijk->k] = IND(IndJ(k, j), IndS(i, j), IndJ(i, k));
						d_ijk->k++;

						d_ijk->data[d_ijk->k].re = -1.0 / sqrt(2.0);
						d_ijk->coord1[d_ijk->k] = IndJ(i, k);
						d_ijk->coord2[d_ijk->k] = IndS(i, j);
						d_ijk->coord3[d_ijk->k] = IndJ(k, j);
						d_ijk->hash[d_ijk->k] = IND(IndJ(i, k), IndS(i, j), IndJ(k, j));
						d_ijk->k++;

						d_ijk->data[d_ijk->k].re = 1.0 / sqrt(2.0);
						d_ijk->coord1[d_ijk->k] = IndJ(k, j);
						d_ijk->coord2[d_ijk->k] = IndJ(i, j);
						d_ijk->coord3[d_ijk->k] = IndS(i, k);
						d_ijk->hash[d_ijk->k] = IND(IndJ(k, j), IndJ(i, j), IndS(i, k));
						d_ijk->k++;

						d_ijk->data[d_ijk->k].re = 1.0 / sqrt(2.0);
						d_ijk->coord1[d_ijk->k] = IndJ(i, k);
						d_ijk->coord2[d_ijk->k] = IndJ(i, j);
						d_ijk->coord3[d_ijk->k] = IndS(k, j);
						d_ijk->hash[d_ijk->k] = IND(IndJ(i, k), IndJ(i, j), IndS(k, j));
						d_ijk->k++;

						d_ijk->data[d_ijk->k].re = 1.0 / sqrt(2.0);
						d_ijk->coord1[d_ijk->k] = IndS(k, j);
						d_ijk->coord2[d_ijk->k] = IndJ(i, j);
						d_ijk->coord3[d_ijk->k] = IndJ(i, k);
						d_ijk->hash[d_ijk->k] = IND(IndS(k, j), IndJ(i, j), IndJ(i, k));
						d_ijk->k++;

						d_ijk->data[d_ijk->k].re = 1.0 / sqrt(2.0);
						d_ijk->coord1[d_ijk->k] = IndS(i, k);
						d_ijk->coord2[d_ijk->k] = IndJ(i, j);
						d_ijk->coord3[d_ijk->k] = IndJ(k, j);
						d_ijk->hash[d_ijk->k] = IND(IndS(i, k), IndJ(i, j), IndJ(k, j));
						d_ijk->k++;
					}
					if (k < i)
					{
						d_ijk->data[d_ijk->k].re = 1.0 / sqrt(2.0);
						d_ijk->coord1[d_ijk->k] = IndS(k, j);
						d_ijk->coord2[d_ijk->k] = IndS(i, j);
						d_ijk->coord3[d_ijk->k] = IndS(k, i);
						d_ijk->hash[d_ijk->k] = IND(IndS(k, j), IndS(i, j), IndS(k, i));
						d_ijk->k++;

						d_ijk->data[d_ijk->k].re = 1.0 / sqrt(2.0);
						d_ijk->coord1[d_ijk->k] = IndS(k, i);
						d_ijk->coord2[d_ijk->k] = IndS(i, j);
						d_ijk->coord3[d_ijk->k] = IndS(k, j);
						d_ijk->hash[d_ijk->k] = IND(IndS(k, i), IndS(i, j), IndS(k, j));
						d_ijk->k++;

						d_ijk->data[d_ijk->k].re = 1.0 / sqrt(2.0);
						d_ijk->coord1[d_ijk->k] = IndJ(k, j);
						d_ijk->coord2[d_ijk->k] = IndS(i, j);
						d_ijk->coord3[d_ijk->k] = IndJ(k, i);
						d_ijk->hash[d_ijk->k] = IND(IndJ(k, j), IndS(i, j), IndJ(k, i));
						d_ijk->k++;

						d_ijk->data[d_ijk->k].re = 1.0 / sqrt(2.0);
						d_ijk->coord1[d_ijk->k] = IndJ(k, i);
						d_ijk->coord2[d_ijk->k] = IndS(i, j);
						d_ijk->coord3[d_ijk->k] = IndJ(k, j);
						d_ijk->hash[d_ijk->k] = IND(IndJ(k, i), IndS(i, j), IndJ(k, j));
						d_ijk->k++;

						d_ijk->data[d_ijk->k].re = 1.0 / sqrt(2.0);
						d_ijk->coord1[d_ijk->k] = IndJ(k, j);
						d_ijk->coord2[d_ijk->k] = IndJ(i, j);
						d_ijk->coord3[d_ijk->k] = IndS(k, i);
						d_ijk->hash[d_ijk->k] = IND(IndJ(k, j), IndJ(i, j), IndS(k, i));
						d_ijk->k++;

						d_ijk->data[d_ijk->k].re = -1.0 / sqrt(2.0);
						d_ijk->coord1[d_ijk->k] = IndJ(k, i);
						d_ijk->coord2[d_ijk->k] = IndJ(i, j);
						d_ijk->coord3[d_ijk->k] = IndS(k, j);
						d_ijk->hash[d_ijk->k] = IND(IndJ(k, i), IndJ(i, j), IndS(k, j));
						d_ijk->k++;

						d_ijk->data[d_ijk->k].re = -1.0 / sqrt(2.0);
						d_ijk->coord1[d_ijk->k] = IndS(k, j);
						d_ijk->coord2[d_ijk->k] = IndJ(i, j);
						d_ijk->coord3[d_ijk->k] = IndJ(k, i);
						d_ijk->hash[d_ijk->k] = IND(IndS(k, j), IndJ(i, j), IndJ(k, i));
						d_ijk->k++;

						d_ijk->data[d_ijk->k].re = 1.0 / sqrt(2.0);
						d_ijk->coord1[d_ijk->k] = IndS(k, i);
						d_ijk->coord2[d_ijk->k] = IndJ(i, j);
						d_ijk->coord3[d_ijk->k] = IndJ(k, j);
						d_ijk->hash[d_ijk->k] = IND(IndS(k, i), IndJ(i, j), IndJ(k, j));
						d_ijk->k++;
					}
				}
			}
		}
	}


	for (int j = 2; j < N; j++)
	{
		int i = 1;
		d_ijk->data[d_ijk->k].re = 2.0 / sqrt((double)(j) * (j + 1));
		d_ijk->coord1[d_ijk->k] = IndD(j);
		d_ijk->coord2[d_ijk->k] = IndD(i);
		d_ijk->coord3[d_ijk->k] = IndD(i);
		d_ijk->hash[d_ijk->k] = IND(IndD(j), IndD(i), IndD(i));
		d_ijk->k++;

		d_ijk->data[d_ijk->k].re = 2.0 / sqrt((double)(j) * (j + 1));
		d_ijk->coord1[d_ijk->k] = IndD(i);
		d_ijk->coord2[d_ijk->k] = IndD(i);
		d_ijk->coord3[d_ijk->k] = IndD(j);
		d_ijk->hash[d_ijk->k] = IND(IndD(i), IndD(i), IndD(j));
		d_ijk->k++;

		d_ijk->data[d_ijk->k].re = 2.0 / sqrt((double)(j) * (j + 1));
		d_ijk->coord1[d_ijk->k] = IndD(i);
		d_ijk->coord2[d_ijk->k] = IndD(j);
		d_ijk->coord3[d_ijk->k] = IndD(i);
		d_ijk->hash[d_ijk->k] = IND(IndD(i), IndD(j), IndD(i));
		d_ijk->k++;
	}

	for (int i = 2; i < N; i++)
	{
		d_ijk->data[d_ijk->k].re = 2.0 * (1 - i) / sqrt((double)(i) * (i + 1));
		d_ijk->coord1[d_ijk->k] = IndD(i);
		d_ijk->coord2[d_ijk->k] = IndD(i);
		d_ijk->coord3[d_ijk->k] = IndD(i);
		d_ijk->hash[d_ijk->k] = IND(IndD(i), IndD(i), IndD(i));
		d_ijk->k++;

		for (int j = i + 1; j < N; j++)
		{
			d_ijk->data[d_ijk->k].re = 2.0 / sqrt((double)(j) * (j + 1));
			d_ijk->coord1[d_ijk->k] = IndD(j);
			d_ijk->coord2[d_ijk->k] = IndD(i);
			d_ijk->coord3[d_ijk->k] = IndD(i);
			d_ijk->hash[d_ijk->k] = IND(IndD(j), IndD(i), IndD(i));
			d_ijk->k++;

			d_ijk->data[d_ijk->k].re = 2.0 / sqrt((double)(j) * (j + 1));
			d_ijk->coord1[d_ijk->k] = IndD(i);
			d_ijk->coord2[d_ijk->k] = IndD(i);
			d_ijk->coord3[d_ijk->k] = IndD(j);
			d_ijk->hash[d_ijk->k] = IND(IndD(i), IndD(i), IndD(j));
			d_ijk->k++;

			d_ijk->data[d_ijk->k].re = 2.0 / sqrt((double)(j) * (j + 1));
			d_ijk->coord1[d_ijk->k] = IndD(i);
			d_ijk->coord2[d_ijk->k] = IndD(j);
			d_ijk->coord3[d_ijk->k] = IndD(i);
			d_ijk->hash[d_ijk->k] = IND(IndD(i), IndD(j), IndD(i));
			d_ijk->k++;
		}
	}
}

ulli fijk_coord_sym(crsMatrix *sel, int N)
{
	ulli cnt = 0, ind;

	for (int i = 1; i < N; i++)
	{
		for (int j = i + 1; j < N; j++)
		{
			ind = IndD(i);
			cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

			ind = IndD(i);
			cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

			ind = IndJ(i, j);
			cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

			ind = IndS(i, j);
			cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

			ind = IndJ(i, j);
			cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

			ind = IndS(i, j);
			cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];
		}
	}

	for (int i = 0; i < N; i++)
	{
		for (int j = i + 1; j < N; j++)
		{
			for (int m = i + 1; m < j; m++)
			{
				ind = IndD(m);
				cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

				ind = IndD(m);
				cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

				ind = IndJ(i, j);
				cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

				ind = IndS(i, j);
				cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

				ind = IndJ(i, j);
				cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

				ind = IndS(i, j);
				cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];
			}
		}
	}

	for (int i = 0; i < N; i++)
	{
		for (int j = i + 1; j < N; j++)
		{
			ind = IndD(j);
			cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

			ind = IndD(j);
			cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

			ind = IndJ(i, j);
			cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

			ind = IndS(i, j);
			cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

			ind = IndJ(i, j);
			cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

			ind = IndS(i, j);
			cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];
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
						ind = IndJ(j, k);
						cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						ind = IndJ(i, k);
						cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						ind = IndS(j, k);
						cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						ind = IndS(i, k);
						cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						ind = IndS(j, k);
						cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						ind = IndS(i, k);
						cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						ind = IndJ(j, k);
						cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						ind = IndJ(i, k);
						cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];
					}
				}
				if (k < j)
				{
					if (k > i)
					{
						ind = IndJ(k, j);
						cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						ind = IndJ(i, k);
						cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						ind = IndS(k, j);
						cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						ind = IndS(i, k);
						cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						ind = IndS(k, j);
						cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						ind = IndS(i, k);
						cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						ind = IndJ(k, j);
						cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						ind = IndJ(i, k);
						cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];
					}
					if (k < i)
					{
						ind = IndJ(k, j);
						cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						ind = IndJ(k, i);
						cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						ind = IndS(k, j);
						cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						ind = IndS(k, i);
						cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						ind = IndS(k, j);
						cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						ind = IndS(k, i);
						cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						ind = IndJ(k, j);
						cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						ind = IndJ(k, i);
						cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];
					}
				}
			}
		}
	}
	return cnt;
}
ulli dijk_coord_sym(crsMatrix *sel, int N)
{
	ulli cnt = 0, ind;

	for (int i = 1; i < N; i++)
	{
		for (int j = i + 1; j < N; j++)
		{
			ind = IndD(i);
			cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

			ind = IndS(i, j);
			cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

			ind = IndD(i);
			cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

			ind = IndJ(i, j);
			cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

			ind = IndS(i, j);
			cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

			ind = IndJ(i, j);
			cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];
		}
	}

	for (int i = 0; i < N; i++)
	{
		for (int j = i + 1; j < N; j++)
		{
			for (int m = i + 1; m < j; m++)
			{
				ind = IndD(m);
				cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

				ind = IndS(i, j);
				cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

				ind = IndD(m);
				cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

				ind = IndJ(i, j);
				cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

				ind = IndS(i, j);
				cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

				ind = IndJ(i, j);
				cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];
			}
		}
	}

	for (int j = 2; j < N; j++)
	{
		int i = 0;
		ind = IndD(j);
		cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

		ind = IndS(i, j);
		cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

		ind = IndD(j);
		cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

		ind = IndJ(i, j);
		cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

		ind = IndS(i, j);
		cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

		ind = IndJ(i, j);
		cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];
	}

	for (int i = 1; i < N; i++)
	{
		for (int j = i + 1; j < N; j++)
		{
			ind = IndD(j);
			cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

			ind = IndS(i, j);
			cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

			ind = IndD(j);
			cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

			ind = IndJ(i, j);
			cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

			ind = IndS(i, j);
			cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

			ind = IndJ(i, j);
			cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];
		}
	}

	for (int i = 0; i < N; i++)
	{
		for (int j = i + 1; j < N; j++)
		{
			for (int z = j + 1; z < N; z++)
			{
				ind = IndD(z);
				cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

				ind = IndS(i, j);
				cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

				ind = IndD(z);
				cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

				ind = IndJ(i, j);
				cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

				ind = IndS(i, j);
				cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

				ind = IndJ(i, j);
				cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];
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
						ind = IndS(j, k);
						cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						ind = IndS(i, k);
						cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						ind = IndJ(j, k);
						cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						ind = IndJ(i, k);
						cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						ind = IndJ(j, k);
						cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						ind = IndJ(i, k);
						cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						ind = IndS(j, k);
						cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						ind = IndS(i, k);
						cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];
					}
				}
				if (k < j)
				{
					if (k > i)
					{
						ind = IndS(k, j);
						cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						ind = IndS(i, k);
						cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						ind = IndJ(k, j);
						cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						ind = IndJ(i, k);
						cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						ind = IndJ(k, j);
						cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						ind = IndJ(i, k);
						cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						ind = IndS(k, j);
						cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						ind = IndS(i, k);
						cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];
					}
					if (k < i)
					{
						ind = IndS(k, j);
						cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						ind = IndS(k, i);
						cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						ind = IndJ(k, j);
						cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						ind = IndJ(k, i);
						cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						ind = IndJ(k, j);
						cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						ind = IndJ(k, i);
						cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						ind = IndS(k, j);
						cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						ind = IndS(k, i);
						cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];
					}
				}
			}
		}
	}


	for (int j = 2; j < N; j++)
	{
		int i = 1;
		ind = IndD(j);
		cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

		ind = IndD(i);
		cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

		ind = IndD(i);
		cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];
	}

	for (int i = 2; i < N; i++)
	{
		ind = IndD(i);
		cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

		for (int j = i + 1; j < N; j++)
		{
			ind = IndD(j);
			cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

			ind = IndD(i);
			cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

			ind = IndD(i);
			cnt += sel->RowIndex[ind + 1] - sel->RowIndex[ind];
		}
	}
	return cnt;
}

void fijk_coord_ch(Tensor_Coordinates * f_ijk, crsMatrix *sel, ulli NZ, int N)
{
	ulli ind;
	ulli size = NZ + 1;
	ulli step1 = N*N;
	ulli step2 = step1 * step1;

	double tmp = 0.0;
	f_ijk->data = new dcomplex[size];
	memset(f_ijk->data, 0, size * sizeof(dcomplex));
	f_ijk->coord1 = new unsigned int[size];
	memset(f_ijk->coord1, 0, size * sizeof(unsigned int));
	f_ijk->coord2 = new unsigned int[size];
	memset(f_ijk->coord2, 0, size * sizeof(unsigned int));
	f_ijk->coord3 = new unsigned int[size];
	memset(f_ijk->coord3, 0, size * sizeof(unsigned int));
	f_ijk->hash = new ulli[size];
	memset(f_ijk->hash, 0, size * sizeof(unsigned int));
	f_ijk->k = 0;

	for (int i = 1; i < N; i++)
	{
		for (int j = i + 1; j < N; j++)
		{
			f_ijk->data[f_ijk->k].re = (i) / sqrt((double)(i)* (i + 1));
			ind = f_ijk->coord1[f_ijk->k] = IndD(i);
			f_ijk->coord2[f_ijk->k] = IndJ(i, j);
			f_ijk->coord3[f_ijk->k] = IndS(i, j);
			f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
			f_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

			f_ijk->data[f_ijk->k].re = -(i) / sqrt((double)(i)* (i + 1));
			ind = f_ijk->coord1[f_ijk->k] = IndD(i);
			f_ijk->coord2[f_ijk->k] = IndS(i, j);
			f_ijk->coord3[f_ijk->k] = IndJ(i, j);
			f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
			f_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

			f_ijk->data[f_ijk->k].re = -(i) / sqrt((double)(i)* (i + 1));
			ind = f_ijk->coord1[f_ijk->k] = IndJ(i, j);
			f_ijk->coord2[f_ijk->k] = IndD(i);
			f_ijk->coord3[f_ijk->k] = IndS(i, j);
			f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
			f_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

			f_ijk->data[f_ijk->k].re = (i) / sqrt((double)(i)* (i + 1));
			ind = f_ijk->coord1[f_ijk->k] = IndS(i, j);
			f_ijk->coord2[f_ijk->k] = IndD(i);
			f_ijk->coord3[f_ijk->k] = IndJ(i, j);
			f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
			f_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

			f_ijk->data[f_ijk->k].re = (i) / sqrt((double)(i)* (i + 1));
			ind = f_ijk->coord1[f_ijk->k] = IndJ(i, j);
			f_ijk->coord2[f_ijk->k] = IndS(i, j);
			f_ijk->coord3[f_ijk->k] = IndD(i);
			f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
			f_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

			f_ijk->data[f_ijk->k].re = -(i) / sqrt((double)(i)* (i + 1));
			ind = f_ijk->coord1[f_ijk->k] = IndS(i, j);
			f_ijk->coord2[f_ijk->k] = IndJ(i, j);
			f_ijk->coord3[f_ijk->k] = IndD(i);
			f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
			f_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];
		}
	}

	for (int i = 0; i < N; i++)
	{
		for (int j = i + 1; j < N; j++)
		{
			for (int m = i + 1; m < j; m++)
			{
				f_ijk->data[f_ijk->k].re = -1.0 / sqrt((double)(m)* (m + 1));
				ind = f_ijk->coord1[f_ijk->k] = IndD(m);
				f_ijk->coord2[f_ijk->k] = IndJ(i, j);
				f_ijk->coord3[f_ijk->k] = IndS(i, j);
				f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
				f_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

				f_ijk->data[f_ijk->k].re = 1.0 / sqrt((double)(m)* (m + 1));
				ind = f_ijk->coord1[f_ijk->k] = IndD(m);
				f_ijk->coord2[f_ijk->k] = IndS(i, j);
				f_ijk->coord3[f_ijk->k] = IndJ(i, j);
				f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
				f_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

				f_ijk->data[f_ijk->k].re = 1.0 / sqrt((double)(m)* (m + 1));
				ind = f_ijk->coord1[f_ijk->k] = IndJ(i, j);
				f_ijk->coord2[f_ijk->k] = IndD(m);
				f_ijk->coord3[f_ijk->k] = IndS(i, j);
				f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
				f_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

				f_ijk->data[f_ijk->k].re = -1.0 / sqrt((double)(m)* (m + 1));
				ind = f_ijk->coord1[f_ijk->k] = IndS(i, j);
				f_ijk->coord2[f_ijk->k] = IndD(m);
				f_ijk->coord3[f_ijk->k] = IndJ(i, j);
				f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
				f_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

				f_ijk->data[f_ijk->k].re = -1.0 / sqrt((double)(m)* (m + 1));
				ind = f_ijk->coord1[f_ijk->k] = IndJ(i, j);
				f_ijk->coord2[f_ijk->k] = IndS(i, j);
				f_ijk->coord3[f_ijk->k] = IndD(m);
				f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
				f_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

				f_ijk->data[f_ijk->k].re = 1.0 / sqrt((double)(m)* (m + 1));
				ind = f_ijk->coord1[f_ijk->k] = IndS(i, j);
				f_ijk->coord2[f_ijk->k] = IndJ(i, j);
				f_ijk->coord3[f_ijk->k] = IndD(m);
				f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
				f_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];
			}
		}
	}

	for (int i = 0; i < N; i++)
	{
		for (int j = i + 1; j < N; j++)
		{
			f_ijk->data[f_ijk->k].re = -(1 + j) / sqrt((double)(j)* (j + 1));
			ind = f_ijk->coord1[f_ijk->k] = IndD(j);
			f_ijk->coord2[f_ijk->k] = IndJ(i, j);
			f_ijk->coord3[f_ijk->k] = IndS(i, j);
			f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
			f_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

			f_ijk->data[f_ijk->k].re = (1 + j) / sqrt((double)(j)* (j + 1));
			ind = f_ijk->coord1[f_ijk->k] = IndD(j);
			f_ijk->coord2[f_ijk->k] = IndS(i, j);
			f_ijk->coord3[f_ijk->k] = IndJ(i, j);
			f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
			f_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

			f_ijk->data[f_ijk->k].re = (1 + j) / sqrt((double)(j)* (j + 1));
			ind = f_ijk->coord1[f_ijk->k] = IndJ(i, j);
			f_ijk->coord2[f_ijk->k] = IndD(j);
			f_ijk->coord3[f_ijk->k] = IndS(i, j);
			f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
			f_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

			f_ijk->data[f_ijk->k].re = -(1 + j) / sqrt((double)(j)* (j + 1));
			ind = f_ijk->coord1[f_ijk->k] = IndS(i, j);
			f_ijk->coord2[f_ijk->k] = IndD(j);
			f_ijk->coord3[f_ijk->k] = IndJ(i, j);
			f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
			f_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

			f_ijk->data[f_ijk->k].re = -(1 + j) / sqrt((double)(j)* (j + 1));
			ind = f_ijk->coord1[f_ijk->k] = IndJ(i, j);
			f_ijk->coord2[f_ijk->k] = IndS(i, j);
			f_ijk->coord3[f_ijk->k] = IndD(j);
			f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
			f_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

			f_ijk->data[f_ijk->k].re = (1 + j) / sqrt((double)(j)* (j + 1));
			ind = f_ijk->coord1[f_ijk->k] = IndS(i, j);
			f_ijk->coord2[f_ijk->k] = IndJ(i, j);
			f_ijk->coord3[f_ijk->k] = IndD(j);
			f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
			f_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];
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
						f_ijk->data[f_ijk->k].re = 1.0 / sqrt(2.0); //-i*Jjk*[Sij,Sik]  //and symmetric Ski, Jkj
						ind = f_ijk->coord1[f_ijk->k] = IndJ(j, k);
						f_ijk->coord2[f_ijk->k] = IndS(i, j);
						f_ijk->coord3[f_ijk->k] = IndS(i, k);
						f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
						f_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						f_ijk->data[f_ijk->k].re = 1.0 / sqrt(2.0); //-i*Jik*[Sij,Sjk]  //and symmetric Skj, Jki
						ind = f_ijk->coord1[f_ijk->k] = IndJ(i, k);
						f_ijk->coord2[f_ijk->k] = IndS(i, j);
						f_ijk->coord3[f_ijk->k] = IndS(j, k);
						f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
						f_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						f_ijk->data[f_ijk->k].re = -1.0 / sqrt(2.0); //-i*Sjk*[Sij,Jik]  //and symmetric Skj, Jki
						ind = f_ijk->coord1[f_ijk->k] = IndS(j, k);
						f_ijk->coord2[f_ijk->k] = IndS(i, j);
						f_ijk->coord3[f_ijk->k] = IndJ(i, k);
						f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
						f_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						f_ijk->data[f_ijk->k].re = -1.0 / sqrt(2.0); //-i*Sik*[Sij,Jjk]  //and symmetric Ski, Jkj
						ind = f_ijk->coord1[f_ijk->k] = IndS(i, k);
						f_ijk->coord2[f_ijk->k] = IndS(i, j);
						f_ijk->coord3[f_ijk->k] = IndJ(j, k);
						f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
						f_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						f_ijk->data[f_ijk->k].re = 1.0 / sqrt(2.0); //-i*Sjk*[Jij,Sik]  //and symmetric Ski, Skj
						ind = f_ijk->coord1[f_ijk->k] = IndS(j, k);
						f_ijk->coord2[f_ijk->k] = IndJ(i, j);
						f_ijk->coord3[f_ijk->k] = IndS(i, k);
						f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
						f_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						f_ijk->data[f_ijk->k].re = -1.0 / sqrt(2.0); //-i*Sik*[Jij,Sjk]  //and symmetric Skj, Ski
						ind = f_ijk->coord1[f_ijk->k] = IndS(i, k);
						f_ijk->coord2[f_ijk->k] = IndJ(i, j);
						f_ijk->coord3[f_ijk->k] = IndS(j, k);
						f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
						f_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						f_ijk->data[f_ijk->k].re = 1.0 / sqrt(2.0); //-i*Jjk*[Jij,Jik]  //and symmetric Jkj, Jki
						ind = f_ijk->coord1[f_ijk->k] = IndJ(j, k);
						f_ijk->coord2[f_ijk->k] = IndJ(i, j);
						f_ijk->coord3[f_ijk->k] = IndJ(i, k);
						f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
						f_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						f_ijk->data[f_ijk->k].re = -1.0 / sqrt(2.0); //-i*Jik*[Jij,Jjk]  //and symmetric Jki, Jkj
						ind = f_ijk->coord1[f_ijk->k] = IndJ(i, k);
						f_ijk->coord2[f_ijk->k] = IndJ(i, j);
						f_ijk->coord3[f_ijk->k] = IndJ(j, k);
						f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
						f_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];
					}
				}
				if (k < j)
				{
					if (k > i)
					{
						f_ijk->data[f_ijk->k].re = -1.0 / sqrt(2.0);
						ind = f_ijk->coord1[f_ijk->k] = IndJ(k, j);
						f_ijk->coord2[f_ijk->k] = IndS(i, j);
						f_ijk->coord3[f_ijk->k] = IndS(i, k);
						f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
						f_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						f_ijk->data[f_ijk->k].re = 1.0 / sqrt(2.0);
						ind = f_ijk->coord1[f_ijk->k] = IndJ(i, k);
						f_ijk->coord2[f_ijk->k] = IndS(i, j);
						f_ijk->coord3[f_ijk->k] = IndS(k, j);
						f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
						f_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						f_ijk->data[f_ijk->k].re = -1.0 / sqrt(2.0);
						ind = f_ijk->coord1[f_ijk->k] = IndS(k, j);
						f_ijk->coord2[f_ijk->k] = IndS(i, j);
						f_ijk->coord3[f_ijk->k] = IndJ(i, k);
						f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
						f_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						f_ijk->data[f_ijk->k].re = 1.0 / sqrt(2.0);
						ind = f_ijk->coord1[f_ijk->k] = IndS(i, k);
						f_ijk->coord2[f_ijk->k] = IndS(i, j);
						f_ijk->coord3[f_ijk->k] = IndJ(k, j);
						f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
						f_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						f_ijk->data[f_ijk->k].re = 1.0 / sqrt(2.0);
						ind = f_ijk->coord1[f_ijk->k] = IndS(k, j);
						f_ijk->coord2[f_ijk->k] = IndJ(i, j);
						f_ijk->coord3[f_ijk->k] = IndS(i, k);
						f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
						f_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						f_ijk->data[f_ijk->k].re = -1.0 / sqrt(2.0);
						ind = f_ijk->coord1[f_ijk->k] = IndS(i, k);
						f_ijk->coord2[f_ijk->k] = IndJ(i, j);
						f_ijk->coord3[f_ijk->k] = IndS(k, j);
						f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
						f_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						f_ijk->data[f_ijk->k].re = -1.0 / sqrt(2.0);
						ind = f_ijk->coord1[f_ijk->k] = IndJ(k, j);
						f_ijk->coord2[f_ijk->k] = IndJ(i, j);
						f_ijk->coord3[f_ijk->k] = IndJ(i, k);
						f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
						f_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						f_ijk->data[f_ijk->k].re = 1.0 / sqrt(2.0);
						ind = f_ijk->coord1[f_ijk->k] = IndJ(i, k);
						f_ijk->coord2[f_ijk->k] = IndJ(i, j);
						f_ijk->coord3[f_ijk->k] = IndJ(k, j);
						f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
						f_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];
					}
					if (k < i)
					{
						f_ijk->data[f_ijk->k].re = -1.0 / sqrt(2.0);
						ind = f_ijk->coord1[f_ijk->k] = IndJ(k, j);
						f_ijk->coord2[f_ijk->k] = IndS(i, j);
						f_ijk->coord3[f_ijk->k] = IndS(k, i);
						f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
						f_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						f_ijk->data[f_ijk->k].re = -1.0 / sqrt(2.0);
						ind = f_ijk->coord1[f_ijk->k] = IndJ(k, i);
						f_ijk->coord2[f_ijk->k] = IndS(i, j);
						f_ijk->coord3[f_ijk->k] = IndS(k, j);
						f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
						f_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						f_ijk->data[f_ijk->k].re = 1.0 / sqrt(2.0);
						ind = f_ijk->coord1[f_ijk->k] = IndS(k, j);
						f_ijk->coord2[f_ijk->k] = IndS(i, j);
						f_ijk->coord3[f_ijk->k] = IndJ(k, i);
						f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
						f_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						f_ijk->data[f_ijk->k].re = 1.0 / sqrt(2.0);
						ind = f_ijk->coord1[f_ijk->k] = IndS(k, i);
						f_ijk->coord2[f_ijk->k] = IndS(i, j);
						f_ijk->coord3[f_ijk->k] = IndJ(k, j);
						f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
						f_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						f_ijk->data[f_ijk->k].re = 1.0 / sqrt(2.0);
						ind = f_ijk->coord1[f_ijk->k] = IndS(k, j);
						f_ijk->coord2[f_ijk->k] = IndJ(i, j);
						f_ijk->coord3[f_ijk->k] = IndS(k, i);
						f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
						f_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						f_ijk->data[f_ijk->k].re = -1.0 / sqrt(2.0);
						ind = f_ijk->coord1[f_ijk->k] = IndS(k, i);
						f_ijk->coord2[f_ijk->k] = IndJ(i, j);
						f_ijk->coord3[f_ijk->k] = IndS(k, j);
						f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
						f_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						f_ijk->data[f_ijk->k].re = 1.0 / sqrt(2.0);
						ind = f_ijk->coord1[f_ijk->k] = IndJ(k, j);
						f_ijk->coord2[f_ijk->k] = IndJ(i, j);
						f_ijk->coord3[f_ijk->k] = IndJ(k, i);
						f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
						f_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						f_ijk->data[f_ijk->k].re = -1.0 / sqrt(2.0);
						ind = f_ijk->coord1[f_ijk->k] = IndJ(k, i);
						f_ijk->coord2[f_ijk->k] = IndJ(i, j);
						f_ijk->coord3[f_ijk->k] = IndJ(k, j);
						f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
						f_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];
					}
				}
			}
		}
	}
}
void dijk_coord_ch(Tensor_Coordinates * d_ijk, crsMatrix *sel, ulli NZ, int N)
{
	ulli ind;
	ulli size = NZ + 1;
	ulli step1 = N*N;
	ulli step2 = step1 * step1;

	d_ijk->data = new dcomplex[size];
	memset(d_ijk->data, 0, size * sizeof(dcomplex));
	d_ijk->coord1 = new unsigned int[size];
	memset(d_ijk->coord1, 0, size * sizeof(unsigned int));
	d_ijk->coord2 = new unsigned int[size];
	memset(d_ijk->coord2, 0, size * sizeof(unsigned int));
	d_ijk->coord3 = new unsigned int[size];
	memset(d_ijk->coord3, 0, size * sizeof(unsigned int));
	d_ijk->hash = new ulli[size];
	memset(d_ijk->hash, 0, size * sizeof(unsigned int));
	d_ijk->k = 0;

	for (int i = 1; i < N; i++)
	{
		for (int j = i + 1; j < N; j++)
		{
			d_ijk->data[d_ijk->k].re = -(i) / sqrt((double)(i)* (i + 1));
			ind = d_ijk->coord1[d_ijk->k] = IndD(i);
			d_ijk->coord2[d_ijk->k] = IndS(i, j);
			d_ijk->coord3[d_ijk->k] = IndS(i, j);
			d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
			d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

			d_ijk->data[d_ijk->k].re = -(i) / sqrt((double)(i)* (i + 1));
			ind = d_ijk->coord1[d_ijk->k] = IndS(i, j);
			d_ijk->coord2[d_ijk->k] = IndD(i);
			d_ijk->coord3[d_ijk->k] = IndS(i, j);
			d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
			d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

			d_ijk->data[d_ijk->k].re = -(i) / sqrt((double)(i)* (i + 1));
			ind = d_ijk->coord1[d_ijk->k] = IndD(i);
			d_ijk->coord2[d_ijk->k] = IndJ(i, j);
			d_ijk->coord3[d_ijk->k] = IndJ(i, j);
			d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
			d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

			d_ijk->data[d_ijk->k].re = -(i) / sqrt((double)(i)* (i + 1));
			ind = d_ijk->coord1[d_ijk->k] = IndJ(i, j);
			d_ijk->coord2[d_ijk->k] = IndD(i);
			d_ijk->coord3[d_ijk->k] = IndJ(i, j);
			d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
			d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

			d_ijk->data[d_ijk->k].re = -(i) / sqrt((double)(i)* (i + 1));
			ind = d_ijk->coord1[d_ijk->k] = IndS(i, j);
			d_ijk->coord2[d_ijk->k] = IndS(i, j);
			d_ijk->coord3[d_ijk->k] = IndD(i);
			d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
			d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

			d_ijk->data[d_ijk->k].re = -(i) / sqrt((double)(i)* (i + 1));
			ind = d_ijk->coord1[d_ijk->k] = IndJ(i, j);
			d_ijk->coord2[d_ijk->k] = IndJ(i, j);
			d_ijk->coord3[d_ijk->k] = IndD(i);
			d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
			d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];
		}
	}

	for (int i = 0; i < N; i++)
	{
		for (int j = i + 1; j < N; j++)
		{
			for (int m = i + 1; m < j; m++)
			{
				d_ijk->data[d_ijk->k].re = 1.0 / sqrt((double)(m)* (m + 1));
				ind = d_ijk->coord1[d_ijk->k] = IndD(m);
				d_ijk->coord2[d_ijk->k] = IndS(i, j);
				d_ijk->coord3[d_ijk->k] = IndS(i, j);
				d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
				d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

				d_ijk->data[d_ijk->k].re = 1.0 / sqrt((double)(m)* (m + 1));
				ind = d_ijk->coord1[d_ijk->k] = IndS(i, j);
				d_ijk->coord2[d_ijk->k] = IndD(m);
				d_ijk->coord3[d_ijk->k] = IndS(i, j);
				d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
				d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

				d_ijk->data[d_ijk->k].re = 1.0 / sqrt((double)(m)* (m + 1));
				ind = d_ijk->coord1[d_ijk->k] = IndD(m);
				d_ijk->coord2[d_ijk->k] = IndJ(i, j);
				d_ijk->coord3[d_ijk->k] = IndJ(i, j);
				d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
				d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

				d_ijk->data[d_ijk->k].re = 1.0 / sqrt((double)(m)* (m + 1));
				ind = d_ijk->coord1[d_ijk->k] = IndJ(i, j);
				d_ijk->coord2[d_ijk->k] = IndD(m);
				d_ijk->coord3[d_ijk->k] = IndJ(i, j);
				d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
				d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

				d_ijk->data[d_ijk->k].re = 1.0 / sqrt((double)(m)* (m + 1));
				ind = d_ijk->coord1[d_ijk->k] = IndS(i, j);
				d_ijk->coord2[d_ijk->k] = IndS(i, j);
				d_ijk->coord3[d_ijk->k] = IndD(m);
				d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
				d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

				d_ijk->data[d_ijk->k].re = 1.0 / sqrt((double)(m)* (m + 1));
				ind = d_ijk->coord1[d_ijk->k] = IndJ(i, j);
				d_ijk->coord2[d_ijk->k] = IndJ(i, j);
				d_ijk->coord3[d_ijk->k] = IndD(m);
				d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
				d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];
			}
		}
	}

	for (int j = 2; j < N; j++)
	{
		int i = 0;
		d_ijk->data[d_ijk->k].re = (1 - j) / sqrt((double)(j)* (j + 1));
		ind = d_ijk->coord1[d_ijk->k] = IndD(j);
		d_ijk->coord2[d_ijk->k] = IndS(i, j);
		d_ijk->coord3[d_ijk->k] = IndS(i, j);
		d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
		d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

		d_ijk->data[d_ijk->k].re = (1 - j) / sqrt((double)(j)* (j + 1));
		ind = d_ijk->coord1[d_ijk->k] = IndS(i, j);
		d_ijk->coord2[d_ijk->k] = IndD(j);
		d_ijk->coord3[d_ijk->k] = IndS(i, j);
		d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
		d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

		d_ijk->data[d_ijk->k].re = (1 - j) / sqrt((double)(j)* (j + 1));
		ind = d_ijk->coord1[d_ijk->k] = IndD(j);
		d_ijk->coord2[d_ijk->k] = IndJ(i, j);
		d_ijk->coord3[d_ijk->k] = IndJ(i, j);
		d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
		d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

		d_ijk->data[d_ijk->k].re = (1 - j) / sqrt((double)(j)* (j + 1));
		ind = d_ijk->coord1[d_ijk->k] = IndJ(i, j);
		d_ijk->coord2[d_ijk->k] = IndD(j);
		d_ijk->coord3[d_ijk->k] = IndJ(i, j);
		d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
		d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

		d_ijk->data[d_ijk->k].re = (1 - j) / sqrt((double)(j)* (j + 1));
		ind = d_ijk->coord1[d_ijk->k] = IndS(i, j);
		d_ijk->coord2[d_ijk->k] = IndS(i, j);
		d_ijk->coord3[d_ijk->k] = IndD(j);
		d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
		d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

		d_ijk->data[d_ijk->k].re = (1 - j) / sqrt((double)(j)* (j + 1));
		ind = d_ijk->coord1[d_ijk->k] = IndJ(i, j);
		d_ijk->coord2[d_ijk->k] = IndJ(i, j);
		d_ijk->coord3[d_ijk->k] = IndD(j);
		d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
		d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];
	}

	for (int i = 1; i < N; i++)
	{
		for (int j = i + 1; j < N; j++)
		{
			d_ijk->data[d_ijk->k].re = (1 - j) / sqrt((double)(j)* (j + 1));
			ind = d_ijk->coord1[d_ijk->k] = IndD(j);
			d_ijk->coord2[d_ijk->k] = IndS(i, j);
			d_ijk->coord3[d_ijk->k] = IndS(i, j);
			d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
			d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

			d_ijk->data[d_ijk->k].re = (1 - j) / sqrt((double)(j)* (j + 1));
			ind = d_ijk->coord1[d_ijk->k] = IndS(i, j);
			d_ijk->coord2[d_ijk->k] = IndD(j);
			d_ijk->coord3[d_ijk->k] = IndS(i, j);
			d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
			d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

			d_ijk->data[d_ijk->k].re = (1 - j) / sqrt((double)(j)* (j + 1));
			ind = d_ijk->coord1[d_ijk->k] = IndD(j);
			d_ijk->coord2[d_ijk->k] = IndJ(i, j);
			d_ijk->coord3[d_ijk->k] = IndJ(i, j);
			d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
			d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

			d_ijk->data[d_ijk->k].re = (1 - j) / sqrt((double)(j)* (j + 1));
			ind = d_ijk->coord1[d_ijk->k] = IndJ(i, j);
			d_ijk->coord2[d_ijk->k] = IndD(j);
			d_ijk->coord3[d_ijk->k] = IndJ(i, j);
			d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
			d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

			d_ijk->data[d_ijk->k].re = (1 - j) / sqrt((double)(j)* (j + 1));
			ind = d_ijk->coord1[d_ijk->k] = IndS(i, j);
			d_ijk->coord2[d_ijk->k] = IndS(i, j);
			d_ijk->coord3[d_ijk->k] = IndD(j);
			d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
			d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

			d_ijk->data[d_ijk->k].re = (1 - j) / sqrt((double)(j)* (j + 1));
			ind = d_ijk->coord1[d_ijk->k] = IndJ(i, j);
			d_ijk->coord2[d_ijk->k] = IndJ(i, j);
			d_ijk->coord3[d_ijk->k] = IndD(j);
			d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
			d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];
		}
	}

	for (int i = 0; i < N; i++)
	{
		for (int j = i + 1; j < N; j++)
		{
			for (int z = j + 1; z < N; z++)
			{
				d_ijk->data[d_ijk->k].re = 2.0 / sqrt((double)(z)* (z + 1));
				ind = d_ijk->coord1[d_ijk->k] = IndD(z);
				d_ijk->coord2[d_ijk->k] = IndS(i, j);
				d_ijk->coord3[d_ijk->k] = IndS(i, j);
				d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
				d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

				d_ijk->data[d_ijk->k].re = 2.0 / sqrt((double)(z)* (z + 1));
				ind = d_ijk->coord1[d_ijk->k] = IndS(i, j);
				d_ijk->coord2[d_ijk->k] = IndD(z);
				d_ijk->coord3[d_ijk->k] = IndS(i, j);
				d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
				d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

				d_ijk->data[d_ijk->k].re = 2.0 / sqrt((double)(z)* (z + 1));
				ind = d_ijk->coord1[d_ijk->k] = IndD(z);
				d_ijk->coord2[d_ijk->k] = IndJ(i, j);
				d_ijk->coord3[d_ijk->k] = IndJ(i, j);
				d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
				d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

				d_ijk->data[d_ijk->k].re = 2.0 / sqrt((double)(z)* (z + 1));
				ind = d_ijk->coord1[d_ijk->k] = IndJ(i, j);
				d_ijk->coord2[d_ijk->k] = IndD(z);
				d_ijk->coord3[d_ijk->k] = IndJ(i, j);
				d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
				d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

				d_ijk->data[d_ijk->k].re = 2.0 / sqrt((double)(z)* (z + 1));
				ind = d_ijk->coord1[d_ijk->k] = IndS(i, j);
				d_ijk->coord2[d_ijk->k] = IndS(i, j);
				d_ijk->coord3[d_ijk->k] = IndD(z);
				d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
				d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

				d_ijk->data[d_ijk->k].re = 2.0 / sqrt((double)(z)* (z + 1));
				ind = d_ijk->coord1[d_ijk->k] = IndJ(i, j);
				d_ijk->coord2[d_ijk->k] = IndJ(i, j);
				d_ijk->coord3[d_ijk->k] = IndD(z);
				d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
				d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];
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
						d_ijk->data[d_ijk->k].re = 1.0 / sqrt(2.0);  //Sjk*{Sij,Sik}  //and symmetric Ski, Skj
						ind = d_ijk->coord1[d_ijk->k] = IndS(j, k);
						d_ijk->coord2[d_ijk->k] = IndS(i, j);
						d_ijk->coord3[d_ijk->k] = IndS(i, k);
						d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
						d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						d_ijk->data[d_ijk->k].re = 1.0 / sqrt(2.0);  //Sik*{Sij,Sjk}  //and symmetric Skj, Ski
						ind = d_ijk->coord1[d_ijk->k] = IndS(i, k);
						d_ijk->coord2[d_ijk->k] = IndS(i, j);
						d_ijk->coord3[d_ijk->k] = IndS(j, k);
						d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
						d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						d_ijk->data[d_ijk->k].re = 1.0 / sqrt(2.0);  //Jjk*{Sij,Jik}  //and symmetric Jkj, Jki
						ind = d_ijk->coord1[d_ijk->k] = IndJ(j, k);
						d_ijk->coord2[d_ijk->k] = IndS(i, j);
						d_ijk->coord3[d_ijk->k] = IndJ(i, k);
						d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
						d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						d_ijk->data[d_ijk->k].re = 1.0 / sqrt(2.0);  //Jik*{Sij,Jjk}  //and symmetric Jki, Jkj
						ind = d_ijk->coord1[d_ijk->k] = IndJ(i, k);
						d_ijk->coord2[d_ijk->k] = IndS(i, j);
						d_ijk->coord3[d_ijk->k] = IndJ(j, k);
						d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
						d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						d_ijk->data[d_ijk->k].re = -1.0 / sqrt(2.0);  //Jjk*{Jij,Sik}  //and symmetric Ski, Jkj
						ind = d_ijk->coord1[d_ijk->k] = IndJ(j, k);
						d_ijk->coord2[d_ijk->k] = IndJ(i, j);
						d_ijk->coord3[d_ijk->k] = IndS(i, k);
						d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
						d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						d_ijk->data[d_ijk->k].re = 1.0 / sqrt(2.0);  //Jik*{Jij,Sjk}  //and symmetric Skj, Jki
						ind = d_ijk->coord1[d_ijk->k] = IndJ(i, k);
						d_ijk->coord2[d_ijk->k] = IndJ(i, j);
						d_ijk->coord3[d_ijk->k] = IndS(j, k);
						d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
						d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						d_ijk->data[d_ijk->k].re = 1.0 / sqrt(2.0);  //Sjk*{Jij,Jik}  //and symmetric Skj, Jki
						ind = d_ijk->coord1[d_ijk->k] = IndS(j, k);
						d_ijk->coord2[d_ijk->k] = IndJ(i, j);
						d_ijk->coord3[d_ijk->k] = IndJ(i, k);
						d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
						d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						d_ijk->data[d_ijk->k].re = -1.0 / sqrt(2.0);  //Sik*{Jij,Jjk}  //and symmetric Ski, Jkj
						ind = d_ijk->coord1[d_ijk->k] = IndS(i, k);
						d_ijk->coord2[d_ijk->k] = IndJ(i, j);
						d_ijk->coord3[d_ijk->k] = IndJ(j, k);
						d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
						d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];
					}
				}
				if (k < j)
				{
					if (k > i)
					{
						d_ijk->data[d_ijk->k].re = 1.0 / sqrt(2.0);
						ind = d_ijk->coord1[d_ijk->k] = IndS(k, j);
						d_ijk->coord2[d_ijk->k] = IndS(i, j);
						d_ijk->coord3[d_ijk->k] = IndS(i, k);
						d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
						d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						d_ijk->data[d_ijk->k].re = 1.0 / sqrt(2.0);
						ind = d_ijk->coord1[d_ijk->k] = IndS(i, k);
						d_ijk->coord2[d_ijk->k] = IndS(i, j);
						d_ijk->coord3[d_ijk->k] = IndS(k, j);
						d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
						d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						d_ijk->data[d_ijk->k].re = -1.0 / sqrt(2.0);
						ind = d_ijk->coord1[d_ijk->k] = IndJ(k, j);
						d_ijk->coord2[d_ijk->k] = IndS(i, j);
						d_ijk->coord3[d_ijk->k] = IndJ(i, k);
						d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
						d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						d_ijk->data[d_ijk->k].re = -1.0 / sqrt(2.0);
						ind = d_ijk->coord1[d_ijk->k] = IndJ(i, k);
						d_ijk->coord2[d_ijk->k] = IndS(i, j);
						d_ijk->coord3[d_ijk->k] = IndJ(k, j);
						d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
						d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						d_ijk->data[d_ijk->k].re = 1.0 / sqrt(2.0);
						ind = d_ijk->coord1[d_ijk->k] = IndJ(k, j);
						d_ijk->coord2[d_ijk->k] = IndJ(i, j);
						d_ijk->coord3[d_ijk->k] = IndS(i, k);
						d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
						d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						d_ijk->data[d_ijk->k].re = 1.0 / sqrt(2.0);
						ind = d_ijk->coord1[d_ijk->k] = IndJ(i, k);
						d_ijk->coord2[d_ijk->k] = IndJ(i, j);
						d_ijk->coord3[d_ijk->k] = IndS(k, j);
						d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
						d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						d_ijk->data[d_ijk->k].re = 1.0 / sqrt(2.0);
						ind = d_ijk->coord1[d_ijk->k] = IndS(k, j);
						d_ijk->coord2[d_ijk->k] = IndJ(i, j);
						d_ijk->coord3[d_ijk->k] = IndJ(i, k);
						d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
						d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						d_ijk->data[d_ijk->k].re = 1.0 / sqrt(2.0);
						ind = d_ijk->coord1[d_ijk->k] = IndS(i, k);
						d_ijk->coord2[d_ijk->k] = IndJ(i, j);
						d_ijk->coord3[d_ijk->k] = IndJ(k, j);
						d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
						d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];
					}
					if (k < i)
					{
						d_ijk->data[d_ijk->k].re = 1.0 / sqrt(2.0);
						ind = d_ijk->coord1[d_ijk->k] = IndS(k, j);
						d_ijk->coord2[d_ijk->k] = IndS(i, j);
						d_ijk->coord3[d_ijk->k] = IndS(k, i);
						d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
						d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						d_ijk->data[d_ijk->k].re = 1.0 / sqrt(2.0);
						ind = d_ijk->coord1[d_ijk->k] = IndS(k, i);
						d_ijk->coord2[d_ijk->k] = IndS(i, j);
						d_ijk->coord3[d_ijk->k] = IndS(k, j);
						d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
						d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						d_ijk->data[d_ijk->k].re = 1.0 / sqrt(2.0);
						ind = d_ijk->coord1[d_ijk->k] = IndJ(k, j);
						d_ijk->coord2[d_ijk->k] = IndS(i, j);
						d_ijk->coord3[d_ijk->k] = IndJ(k, i);
						d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
						d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						d_ijk->data[d_ijk->k].re = 1.0 / sqrt(2.0);
						ind = d_ijk->coord1[d_ijk->k] = IndJ(k, i);
						d_ijk->coord2[d_ijk->k] = IndS(i, j);
						d_ijk->coord3[d_ijk->k] = IndJ(k, j);
						d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
						d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						d_ijk->data[d_ijk->k].re = 1.0 / sqrt(2.0);
						ind = d_ijk->coord1[d_ijk->k] = IndJ(k, j);
						d_ijk->coord2[d_ijk->k] = IndJ(i, j);
						d_ijk->coord3[d_ijk->k] = IndS(k, i);
						d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
						d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						d_ijk->data[d_ijk->k].re = -1.0 / sqrt(2.0);
						ind = d_ijk->coord1[d_ijk->k] = IndJ(k, i);
						d_ijk->coord2[d_ijk->k] = IndJ(i, j);
						d_ijk->coord3[d_ijk->k] = IndS(k, j);
						d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
						d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						d_ijk->data[d_ijk->k].re = -1.0 / sqrt(2.0);
						ind = d_ijk->coord1[d_ijk->k] = IndS(k, j);
						d_ijk->coord2[d_ijk->k] = IndJ(i, j);
						d_ijk->coord3[d_ijk->k] = IndJ(k, i);
						d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
						d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

						d_ijk->data[d_ijk->k].re = 1.0 / sqrt(2.0);
						ind = d_ijk->coord1[d_ijk->k] = IndS(k, i);
						d_ijk->coord2[d_ijk->k] = IndJ(i, j);
						d_ijk->coord3[d_ijk->k] = IndJ(k, j);
						d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
						d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];
					}
				}
			}
		}
	}


	for (int j = 2; j < N; j++)
	{
		int i = 1;
		d_ijk->data[d_ijk->k].re = 2.0 / sqrt((double)(j)* (j + 1));
		ind = d_ijk->coord1[d_ijk->k] = IndD(j);
		d_ijk->coord2[d_ijk->k] = IndD(i);
		d_ijk->coord3[d_ijk->k] = IndD(i);
		d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
		d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

		d_ijk->data[d_ijk->k].re = 2.0 / sqrt((double)(j)* (j + 1));
		ind = d_ijk->coord1[d_ijk->k] = IndD(i);
		d_ijk->coord2[d_ijk->k] = IndD(i);
		d_ijk->coord3[d_ijk->k] = IndD(j);
		d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
		d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

		d_ijk->data[d_ijk->k].re = 2.0 / sqrt((double)(j)* (j + 1));
		ind = d_ijk->coord1[d_ijk->k] = IndD(i);
		d_ijk->coord2[d_ijk->k] = IndD(j);
		d_ijk->coord3[d_ijk->k] = IndD(i);
		d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
		d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];
	}

	for (int i = 2; i < N; i++)
	{
		d_ijk->data[d_ijk->k].re = 2.0 * (1 - i) / sqrt((double)(i)* (i + 1));
		ind = d_ijk->coord1[d_ijk->k] = IndD(i);
		d_ijk->coord2[d_ijk->k] = IndD(i);
		d_ijk->coord3[d_ijk->k] = IndD(i);
		d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
		d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

		for (int j = i + 1; j < N; j++)
		{
			d_ijk->data[d_ijk->k].re = 2.0 / sqrt((double)(j)* (j + 1));
			ind = d_ijk->coord1[d_ijk->k] = IndD(j);
			d_ijk->coord2[d_ijk->k] = IndD(i);
			d_ijk->coord3[d_ijk->k] = IndD(i);
			d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
			d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

			d_ijk->data[d_ijk->k].re = 2.0 / sqrt((double)(j)* (j + 1));
			ind = d_ijk->coord1[d_ijk->k] = IndD(i);
			d_ijk->coord2[d_ijk->k] = IndD(i);
			d_ijk->coord3[d_ijk->k] = IndD(j);
			d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
			d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];

			d_ijk->data[d_ijk->k].re = 2.0 / sqrt((double)(j)* (j + 1));
			ind = d_ijk->coord1[d_ijk->k] = IndD(i);
			d_ijk->coord2[d_ijk->k] = IndD(j);
			d_ijk->coord3[d_ijk->k] = IndD(i);
			d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
			d_ijk->k += sel->RowIndex[ind + 1] - sel->RowIndex[ind];
		}
	}
}

void calc_CooQs(int N_mat, Model *m, Tensor_Coordinates * f_ijk, crsMatrix *hMat, crsMatrix *&res)
{
	Tensor_Coordinates * select = new Tensor_Coordinates;
	int *hash = new int[N_mat];
	int *rowi = new int[N_mat + 1];
	dcomplex *hash_calc = new dcomplex[N_mat];
	ulli cnt;

	cnt = fijk_coord_sym(hMat, m->N + 1);
	fijk_coord_ch(select, hMat, cnt + 1, m->N + 1);

	for (ulli i = 0; i < select->k; i++)
	{
		int j = select->coord1[i];
		dcomplex val1 = select->data[i];
		dcomplex val2 = hMat->Value[hMat->RowIndex[j]];
		select->data[i].re = val1.re * val2.re - val1.im * val2.im;
		select->data[i].im = val1.re * val2.im + val1.re * val2.im;

		ulli step = N_mat;
		select->hash[i] = select->coord2[i] + step * select->coord3[i];
	}

	sort_matrix(select);


	for (int i = 0; i < N_mat; i++)
	{
		rowi[i] = -1;
	}
	for (ulli i = 0; i < select->k; i++)
	{
		if (rowi[select->coord3[i]] == -1)
		{
			rowi[select->coord3[i]] = i;
		}
	}
	rowi[N_mat] = select->k;
	for (int i = N_mat - 1; i >0; i--)
	{
		if (rowi[i] == -1)
		{
			rowi[i] = rowi[i + 1];
		}
	}
	rowi[N_mat] = select->k;

	res = new crsMatrix(N_mat, 1);

	int countNotZero = 0;
	for (int j = 0; j < N_mat; j++)
	{
		res->RowIndex[j] = countNotZero;
		memset(hash, 0, sizeof(int) * N_mat);
		int start, finish;
		start = rowi[j];
		finish = rowi[j + 1];
		for (int k = start; k < finish; k++)
		{
			if (hash[select->coord2[k]] == 0)
			{
				hash[select->coord2[k]] = 1;
				countNotZero++;
			}
		}
	}
	res->RowIndex[N_mat] = countNotZero;
	res->setNZ(countNotZero);

	countNotZero = 0;
	for (int j = 0; j < N_mat; j++)
	{
		memset(hash, 0, sizeof(int) * N_mat);
		int start, finish;
		start = rowi[j];
		finish = rowi[j + 1];
		for (int k = start; k < finish; k++)
		{
			if (hash[select->coord2[k]] == 0)
			{
				hash[select->coord2[k]] = 1;
				res->Col[countNotZero] = select->coord2[k];
				countNotZero++;
				hash_calc[select->coord2[k]] = select->data[k];
			}
			else
			{
				hash_calc[select->coord2[k]].re += select->data[k].re;
				hash_calc[select->coord2[k]].im += select->data[k].im;
			}
		}
		start = res->RowIndex[j];
		finish = res->RowIndex[j + 1];
		for (int k = start; k < finish; k++)
		{
			res->Value[k] = hash_calc[res->Col[k]];
		}
	}

	crsMatrix * R = new crsMatrix(*res);
	Transpose(*res, *R, false);
	Transpose(*R, *res, false);

	delete R;
	delete[] hash;
	delete[] hash_calc;
	delete[] rowi;

	free_matrix(select);
	delete select;
}
ulli calcZ_ijk(Tensor_Coordinates *f_ijk, Tensor_Coordinates *d_ijk, Tensor_Coordinates *&Z_ijk)
{
	ulli i1 = f_ijk->hash[0];
	ulli i2 = d_ijk->hash[0];
	ulli j1 = 0, j2 = 0;
	ulli cnt = 0;
	while ((j1 < f_ijk->k) && (j2 < d_ijk->k))
	{
		i1 = f_ijk->hash[j1];
		i2 = d_ijk->hash[j2];
		if (i1 < i2)
		{
			j1++;
		}
		else
		{
			if (i1 > i2)
			{
				j2++;
			}
			else
			{
				j1++;
				j2++;
			}
		}
		cnt++;
	}
	cnt += -(j1 - f_ijk->k) - (j2 - d_ijk->k);
	Z_ijk = create_matrix(cnt);

	j1 = 0; j2 = 0; cnt = 0;
	while ((j1 < f_ijk->k) && (j2 < d_ijk->k))
	{
		i1 = f_ijk->hash[j1];
		i2 = d_ijk->hash[j2];
		if (i1 < i2)
		{
			Z_ijk->coord1[cnt] = f_ijk->coord1[j1];
			Z_ijk->coord2[cnt] = f_ijk->coord2[j1];
			Z_ijk->coord3[cnt] = f_ijk->coord3[j1];
			Z_ijk->data[cnt] = f_ijk->data[j1];
			Z_ijk->hash[cnt] = f_ijk->hash[j1];

			j1++;
		}
		else
		{
			if (i1 > i2)
			{
				Z_ijk->coord1[cnt] = d_ijk->coord1[j2];
				Z_ijk->coord2[cnt] = d_ijk->coord2[j2];
				Z_ijk->coord3[cnt] = d_ijk->coord3[j2];
				Z_ijk->data[cnt].re = -d_ijk->data[j2].im;
				Z_ijk->data[cnt].im = d_ijk->data[j2].re;
				Z_ijk->hash[cnt] = d_ijk->hash[j2];
				j2++;
			}
			else
			{
				Z_ijk->coord1[cnt] = d_ijk->coord1[j2];
				Z_ijk->coord2[cnt] = d_ijk->coord2[j2];
				Z_ijk->coord3[cnt] = d_ijk->coord3[j2];
				Z_ijk->data[cnt] = f_ijk->data[j1];
				Z_ijk->data[cnt].re += -d_ijk->data[j2].im;
				Z_ijk->data[cnt].im += d_ijk->data[j2].re;
				Z_ijk->hash[cnt] = d_ijk->hash[j2];
				j1++;
				j2++;
			}
		}
		cnt++;
	}
	while ((j1 < f_ijk->k))
	{
		Z_ijk->coord1[cnt] = f_ijk->coord1[j1];
		Z_ijk->coord2[cnt] = f_ijk->coord2[j1];
		Z_ijk->coord3[cnt] = f_ijk->coord3[j1];
		Z_ijk->data[cnt] = f_ijk->data[j1];
		Z_ijk->hash[cnt] = f_ijk->hash[j1];
		j1++;
		cnt++;
	}
	while ((j2 < d_ijk->k))
	{
		Z_ijk->coord1[cnt] = d_ijk->coord1[j2];
		Z_ijk->coord2[cnt] = d_ijk->coord2[j2];
		Z_ijk->coord3[cnt] = d_ijk->coord3[j2];
		Z_ijk->data[cnt].re = -d_ijk->data[j2].im;
		Z_ijk->data[cnt].im = d_ijk->data[j2].re;
		Z_ijk->hash[cnt] = d_ijk->hash[j2];
		j2++;
		cnt++;
	}

	return cnt;
}

void createIndex(Tensor_Coordinates *m_ijk, ulli N_Mat, unsigned int * index)
{
	for (int i = 0; i < N_Mat + 1; i++)
	{
		index[i] = N_Mat + 1;
	}
	index[m_ijk->coord2[0]] = 0;

	for (ulli i = 1; i < m_ijk->k; i++)
	{
		if (m_ijk->coord2[i] != m_ijk->coord2[i - 1])
		{
			index[m_ijk->coord2[i]] = i;
		}
	}

	index[N_Mat] = m_ijk->k;
	for (int i = N_Mat - 1; i >= 0; i--)
	{
		if (index[i] == N_Mat + 1)
		{
			index[i] = index[i + 1];
		}
	}

}
ulli multTmpRsSTD(ulli N_mat, unsigned int * indexF, unsigned int * indexZ, Tensor_Coordinates *z_ijk, Tensor_Coordinates *f_ijk, crsMatrix * l_mat, bool swapInd, vector<map<int, dcomplex> > & mat)
{
	ulli step1 = N_mat;

	unsigned int c1, c2, j1, j2;
	ulli cnt = 0;
	dcomplex v;
	for (unsigned int i = 0; i < N_mat; i++)
	{
		c1 = (indexF[i + 1] - indexF[i]);
		c2 = (indexZ[i + 1] - indexZ[i]);
		for (unsigned int j = 0; j < c1 * c2; j++)
		{
			j1 = indexF[i] + j % c1;
			j2 = indexZ[i] + j / c1;
			dcomplex v1 = z_ijk->data[j2];
			dcomplex v2 = f_ijk->data[j1];

			int ind_i, ind_j, ind_tmp;
			ind_i = f_ijk->coord1[j1];
			ind_j = z_ijk->coord1[j2];

			if (swapInd)
			{
				ind_tmp = ind_i; ind_i = ind_j; ind_j = ind_tmp;
			}

			dcomplex v3;// = l_mat->Value[ind_i + ind_j * l_mat->N];
			v3.re = v3.im = 0.0;

			int si = l_mat->RowIndex[ind_j];
			int fi = l_mat->RowIndex[ind_j + 1];
			for (int ind_f = si; ind_f < fi; ind_f++)
			{
				if (l_mat->Col[ind_f] == ind_i)
				{
					v3 = l_mat->Value[ind_f];
					break;
				}
			}

			dcomplex tmp1;
			tmp1.re = v1.re * v2.re - v1.im * v2.im;
			tmp1.im = v1.re * v2.im + v1.im * v2.re;
			tmp1.re = tmp1.re * v3.re - tmp1.im * v3.im;
			tmp1.im = tmp1.re * v3.im + tmp1.im * v3.re;

			if ((tmp1.re != 0.0) || (tmp1.im != 0.0))
			{
				v = mat[f_ijk->coord3[j1]][z_ijk->coord3[j2]];
				v.re += tmp1.re;
				v.im += tmp1.im;

				mat[f_ijk->coord3[j1]][z_ijk->coord3[j2]] = v;
			}
		}
		cnt += c1 * c2;

	}

	return cnt;
}

void to_F_basis_for_zeros(crsMatrix * Mat, crsMatrix * vec)
{
	int N = Mat->N;
	int N_mat = N * N - 1;
	int cnt = 0;
	int *mask = new int[N_mat];
	int *col = new int[N_mat];
	dcomplex *value = new dcomplex[N_mat];
	for (int i = 0; i < N_mat; i++)
	{
		mask[i] = 0;
		value[i].re = 0.0;
		value[i].im = 0.0;
	}

	dcomplex sum;
	sum.re = 0.0;
	sum.im = 0.0;
	for (int i = 0; i < N; i++)
	{
		int start = Mat->RowIndex[i];
		int finish = Mat->RowIndex[i + 1];
		for (int j = start; j < finish; j++)
		{
			int k = Mat->Col[j];

			if (k != i)
			{
				int ii = i;
				int jj = k;
				int z = -1;
				if (ii > jj) {
					ii = k;
					jj = i;
					z = 1;
				}

				int index = ((N - 1 + N - ii) * ii) + jj - ii - 1;
				if (mask[index] != 1)
				{
					col[cnt + 0] = index + 0;
					col[cnt + 1] = index + 1;
					cnt += 2;
				}

				mask[index] = 1;
				mask[index + 1] = 1;

				value[index].re += Mat->Value[j].re / sqrt(2.0);
				value[index].im += Mat->Value[j].im / sqrt(2.0);

				value[index + 1].re += z * Mat->Value[j].im / sqrt(2.0);
				value[index + 1].im += -z * Mat->Value[j].re / sqrt(2.0);

			}
			else
			{
				if (k != 0)
				{
					int index = N * (N - 1) + k - 1;
					mask[index] = 1;
					double value_d = 1.0 / sqrt((double)((k + 0)* (k + 1)));
					value[index].re = sum.re * value_d;
					value[index].re -= Mat->Value[j].re * (k + 0) * value_d;
					value[index].im = sum.im * value_d;
					value[index].im -= Mat->Value[j].im * (k + 0) * value_d;
					col[cnt] = index;
					cnt++;
				}

				sum.re += Mat->Value[j].re;
				sum.im += Mat->Value[j].im;
			}
		}
	}

	std::sort(col, col + cnt);

	vec->NZ = cnt;
	for (int i = 0; i < vec->N + 1; i++)
	{
		vec->RowIndex[i] = 0;
	}
	for (int i = 0; i < cnt; i++)
	{
		vec->Col[i] = col[i];
		vec->Value[i] = value[col[i]];
		vec->RowIndex[col[i] + 1] ++;
	}
	for (int i = 0; i < vec->N; i++)
	{
		vec->RowIndex[i + 1] = vec->RowIndex[i] + vec->RowIndex[i + 1];
	}

	delete[] mask;
	delete[] col;
	delete[] value;
}
void to_F_basis(crsMatrix * Mat, crsMatrix * vec)
{
	int N = Mat->N;
	int N_mat = N * N - 1;
	int cnt = 0;
	int *mask = new int[N_mat];
	int *col = new int[N_mat];
	dcomplex *value = new dcomplex[N_mat];
	for (int i = 0; i < N_mat; i++)
	{
		mask[i] = 0;
		value[i].re = 0.0;
		value[i].im = 0.0;
	}

	dcomplex sum;
	sum.re = 0.0;
	sum.im = 0.0;
	for (int i = 0; i < N; i++)
	{
		int start = Mat->RowIndex[i];
		int finish = Mat->RowIndex[i + 1];
		for (int j = start; j < finish; j++)
		{
			int k = Mat->Col[j];

			if (k != i)
			{
				int ii = i;
				int jj = k;
				int z = -1;
				if (ii > jj) {
					ii = k;
					jj = i;
					z = 1;
				}

				int index = ((N - 1 + N - ii) * ii) + (jj - ii - 1) * 2;
				if (mask[index] != 1)
				{
					col[cnt + 0] = index + 0;
					col[cnt + 1] = index + 1;
					cnt += 2;
				}

				mask[index] = 1;
				mask[index + 1] = 1;

				value[index].re += Mat->Value[j].re / sqrt(2.0);
				value[index].im += Mat->Value[j].im / sqrt(2.0);

				value[index + 1].re += z * Mat->Value[j].im / sqrt(2.0);
				value[index + 1].im += -z * Mat->Value[j].re / sqrt(2.0);

			}
			else
			{
				if (k != 0)
				{
					int index = N * (N - 1) + k - 1;
					mask[index] = 1;
					double value_d = 1.0 / sqrt((double)((k + 0)* (k + 1)));
					value[index].re = sum.re * value_d;
					value[index].re -= Mat->Value[j].re * (k + 0) * value_d;
					value[index].im = sum.im * value_d;
					value[index].im -= Mat->Value[j].im * (k + 0) * value_d;
					col[cnt] = index;
					cnt++;
				}

				sum.re += Mat->Value[j].re;
				sum.im += Mat->Value[j].im;
			}
		}
	}

	for (int i = 0; i < cnt; i++)
	{
		if ((value[col[i]].re == 0.0) && (value[col[i]].im == 0.0))
		{
			cnt--;
			col[i] = col[cnt];
			i--;
		}
	}
	std::sort(col, col + cnt);

	vec->NZ = cnt;
	for (int i = 0; i < vec->N + 1; i++)
	{
		vec->RowIndex[i] = 0;
	}
	for (int i = 0; i < cnt; i++)
	{
		vec->Col[i] = col[i];
		vec->Value[i] = value[col[i]];
		vec->RowIndex[col[i] + 1] ++;
	}
	for (int i = 0; i < vec->N; i++)
	{
		vec->RowIndex[i + 1] = vec->RowIndex[i] + vec->RowIndex[i + 1];
	}

	delete[] mask;
	delete[] col;
	delete[] value;
}

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

void complex_to_real(dcomplex *mat, int N)
{
	double *value = (double *)(mat);
	for (int i = 0; i <N; i++)
	{
		value[i] = mat[i].re;
	}
}
void real_to_complex(dcomplex *mat, int N)
{
	double *value = (double *)(mat);
	for (int i = N - 1; i >= 0; i--)
	{
		mat[i].re = value[i];
		mat[i].im = 0;
	}
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

	model->h_J = new crsMatrix(model->N_mat, model->N_mat);
	model->h_h_z = new crsMatrix(model->N_mat, model->N_mat);
	model->h_h_x = new crsMatrix(model->N_mat, model->N_mat);
	model->h_s_x = new crsMatrix(model->N_mat, model->N_mat);

	model->H_J = new crsMatrix(model->N_mat, model->N_mat);
	model->H_h_z = new crsMatrix(model->N_mat, model->N_mat);
	model->H_h_x = new crsMatrix(model->N_mat, model->N_mat);
	model->H_s_x = new crsMatrix(model->N_mat, model->N_mat);

	model->H0 = NULL;
	model->H1 = NULL;

	model->Q_J_s = NULL;
	model->Q_h_z_s = NULL;
	model->Q_h_x_s = NULL;
	model->Q_s_x_s = NULL;

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

	model->f_ijk = NULL;

	model->l_mat = NULL;

	return model;
}
void freeModel(Model * model)
{
	if (model->h_J != NULL)
	{
		delete model->h_J;
	}
	if (model->h_h_z != NULL)
	{
		delete model->h_h_z;
	}
	if (model->h_h_x != NULL)
	{
		delete model->h_h_x;
	}
	if (model->h_s_x != NULL)
	{
		delete model->h_s_x;
	}

	if (model->H_J != NULL)
	{
		delete model->H_J;
	}
	if (model->H_h_z != NULL)
	{
		delete model->H_h_z;
	}
	if (model->H_h_x != NULL)
	{
		delete model->H_h_x;
	}
	if (model->H_s_x != NULL)
	{
		delete model->H_s_x;
	}

	if (model->H0 != NULL)
	{
		delete model->H0;
	}
	if (model->H1 != NULL)
	{
		delete model->H1;
	}

	if (model->Q_J_s != NULL)
	{
		delete model->Q_J_s;
	}
	if (model->Q_h_z_s != NULL)
	{
		delete model->Q_h_z_s;
	}
	if (model->Q_h_x_s != NULL)
	{
		delete model->Q_h_x_s;
	}
	if (model->Q_s_x_s != NULL)
	{
		delete model->Q_s_x_s;
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

	if (model->l_mat != NULL)
	{
		delete model->l_mat;
	}

	if (model->f_ijk != NULL)
	{
		free_matrix(model->f_ijk);
		delete model->f_ijk;
	}

	delete model;
}

crsMatrix * create_sigma_i_x_matrix(Model * m, ConfigParam &cp, int bit_id)
{
	int N = m->N;

	int * cols_ids = new int[N + 1];
	for (int col_id = 0; col_id < N + 1; col_id++)
	{
		cols_ids[col_id] = 0;
	}
	fill_cols_h_x(cols_ids, cp.N, bit_id);

	crsMatrix * sigma_i_x = new crsMatrix(N + 1, N + 1);

	for (int i = 0; i < N + 1; i++)
	{
		sigma_i_x->RowIndex[i] = i;
		sigma_i_x->Col[i] = cols_ids[i];
		sigma_i_x->Value[i].re = 1.0;
	}
	sigma_i_x->RowIndex[N + 1] = N + 1;

	delete[] cols_ids;

	return sigma_i_x;
}
crsMatrix * create_sigma_i_y_matrix(Model * m, ConfigParam &cp, int bit_id)
{
	int N = m->N;

	int * cols_ids = new int[N + 1];
	for (int col_id = 0; col_id < N + 1; col_id++)
	{
		cols_ids[col_id] = 0;
	}
	fill_cols_h_x(cols_ids, cp.N, bit_id);

	crsMatrix * sigma_i_y = new crsMatrix(N + 1, N + 1);

	for (int i = 0; i < N + 1; i++)
	{
		sigma_i_y->RowIndex[i] = i;
		sigma_i_y->Col[i] = cols_ids[i];

		if (sigma_i_y->Col[i] > i)
		{
			sigma_i_y->Value[i].im = -1.0;
		}
		else
		{
			sigma_i_y->Value[i].im = 1.0;
		}

	}
	sigma_i_y->RowIndex[N + 1] = N + 1;

	delete[] cols_ids;

	return sigma_i_y;
}
crsMatrix * create_sigma_i_z_matrix(Model * m, ConfigParam &cp, int bit_id)
{
	int N = m->N;

	crsMatrix * sigma_i_z = new crsMatrix(N + 1, N + 1);
	int * RowIndex = sigma_i_z->RowIndex;
	int * Col = sigma_i_z->Col;
	dcomplex * Value = sigma_i_z->Value;

	double * diag_array = new double[N + 1];

	for (int state_id = 0; state_id < N + 1; state_id++)
	{
		int bit = bit_at(state_id, bit_id);

		if (bit == 0)
		{
			diag_array[state_id] = 1.0;
		}
		else
		{
			diag_array[state_id] = -1.0;
		}
	}

	int i, j;
	RowIndex[0] = 0;
	for (i = 0; i < N + 1; i++)
	{
		RowIndex[i] = i;
		Col[i] = i;
	}
	RowIndex[i] = i;

	for (i = 0; i < N + 1; i++)
	{
		Value[i].re = diag_array[i];
	}

	delete[] diag_array;

	return sigma_i_z;
}

crsMatrix * create_H_J_matrix(Model * m, RunParam &rp, ConfigParam &cp, MainData &md)
{
	int N = m->N;

	crsMatrix * H_J = new crsMatrix(N + 1, N + 1);
	int * RowIndex = H_J->RowIndex;
	int * Col = H_J->Col;
	dcomplex * Value = H_J->Value;

	double * diag_array = new double[N + 1];

	for (int state_id = 0; state_id < N + 1; state_id++)
	{
		diag_array[state_id] = 0.0;

		for (int bit_id = 0; bit_id < cp.N - 1; bit_id++)
		{
			int first_bit = bit_at(state_id, bit_id);
			int second_bit = bit_at(state_id, bit_id + 1);

			if (first_bit == second_bit)
			{
				diag_array[state_id] += -(1.0 * md.disorder_J[bit_id]);
			}
			else
			{
				diag_array[state_id] += -(-1.0 * md.disorder_J[bit_id]);
			}
		}
	}

	int i, j;
	RowIndex[0] = 0;
	for (i = 0; i < N + 1; i++)
	{
		RowIndex[i] = i;
		Col[i] = i;
	}
	RowIndex[i] = i;

	for (i = 0; i < N + 1; i++)
	{
		Value[i].re = diag_array[i];
	}

	delete[] diag_array;

	if (rp.debug == 1)
	{
		string fn = rp.path + "H_J" + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, H_J, 16, false);
	}

	return H_J;
}
void init_h_J_vector(Model * m, RunParam &rp, ConfigParam &cp, MainData &md)
{
	crsMatrix * H_J = create_H_J_matrix(m, rp, cp, md), *res;
	m->H_J = H_J;
	int N = m->N;
	crsMatrix * h_J = m->h_J;

	to_F_basis(H_J, h_J);
}

crsMatrix * create_H_h_z_matrix(Model * m, RunParam &rp, ConfigParam &cp, MainData &md)
{
	int N = m->N;

	crsMatrix * H_h_z = new crsMatrix(N + 1, N + 1);
	int * RowIndex = H_h_z->RowIndex;
	int * Col = H_h_z->Col;
	dcomplex * Value = H_h_z->Value;

	double * diag_array = new double[N + 1];

	for (int state_id = 0; state_id < N + 1; state_id++)
	{
		diag_array[state_id] = 0.0;

		for (int bit_id = 0; bit_id < cp.N; bit_id++)
		{
			int bit = bit_at(state_id, bit_id);
			
			if (bit == 0)
			{
				diag_array[state_id] += -(1.0 * md.disorder_h_z[bit_id]);
			}
			else
			{
				diag_array[state_id] += -(-1.0 * md.disorder_h_z[bit_id]);
			}
		}
	}

	int i, j;
	RowIndex[0] = 0;
	for (i = 0; i < N + 1; i++)
	{
		RowIndex[i] = i;
		Col[i] = i;
	}
	RowIndex[i] = i;

	for (i = 0; i < N + 1; i++)
	{
		Value[i].re = diag_array[i];
	}

	delete[] diag_array;

	if (rp.debug == 1)
	{
		string fn = rp.path + "H_h_z" + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, H_h_z, 16, false);
	}

	return H_h_z;
}
void init_h_h_z_vector(Model * m, RunParam &rp, ConfigParam &cp, MainData &md)
{
	crsMatrix * H_h_z = create_H_h_z_matrix(m, rp, cp, md), *res;
	m->H_h_z = H_h_z;
	int N = m->N;

	crsMatrix * h_h_z = m->h_h_z;

	to_F_basis(H_h_z, h_h_z);
}

crsMatrix * create_H_h_x_matrix(Model * m, RunParam &rp, ConfigParam &cp, MainData &md)
{
	int N = m->N;

	crsMatrix * H_h_x = new crsMatrix(N + 1, N + 1);
	for (int i = 0; i < N + 1; i++)
	{
		H_h_x->RowIndex[i] = i;
		H_h_x->Col[i] = i;
	}
	H_h_x->RowIndex[N + 1] = N + 1;


	for (int bit_id = 0; bit_id < cp.N; bit_id++)
	{
		int * cols_ids = new int[N + 1];
		for (int col_id = 0; col_id < N + 1; col_id++)
		{
			cols_ids[col_id] = 0;
		}
		fill_cols_h_x(cols_ids, cp.N, bit_id);

		crsMatrix * tmp_mtx_1 = new crsMatrix(N + 1, N + 1);

		for (int i = 0; i < N + 1; i++)
		{
			tmp_mtx_1->RowIndex[i] = i;
			tmp_mtx_1->Col[i] = cols_ids[i];
			tmp_mtx_1->Value[i].re = -md.disorder_h_x[bit_id];
		}
		tmp_mtx_1->RowIndex[N + 1] = N + 1;


		crsMatrix * tmp_mtx_2 = new crsMatrix(*H_h_x);

		dcomplex beta;
		beta.re = 1.0;
		beta.im = 0.0;

		SparseMKLAdd(*tmp_mtx_1, beta, *tmp_mtx_2, *H_h_x);

		delete tmp_mtx_1;
		delete tmp_mtx_2;

		delete[] cols_ids;
	}

	if (rp.debug == 1)
	{
		string fn = rp.path + "H_h_x" + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, H_h_x, 16, false);
	}

	return H_h_x;
}
void init_h_h_x_vector(Model * m, RunParam &rp, ConfigParam &cp, MainData &md)
{
	crsMatrix * H_h_x = create_H_h_x_matrix(m, rp, cp, md), *res;
	m->H_h_x = H_h_x;
	int N = m->N;
	crsMatrix * h_h_x = m->h_h_x;

	if (fabs(cp.h_x) > 1.0e-14)
	{
		to_F_basis(H_h_x, h_h_x);
	}
	else
	{
		to_F_basis_for_zeros(H_h_x, h_h_x);
	}
}

crsMatrix * create_H_s_x_matrix(Model * m, RunParam &rp, ConfigParam &cp, MainData &md)
{
	int N = m->N;

	crsMatrix * H_s_x = new crsMatrix(N + 1, N + 1);
	for (int i = 0; i < N + 1; i++)
	{
		H_s_x->RowIndex[i] = i;
		H_s_x->Col[i] = i;
	}
	H_s_x->RowIndex[N + 1] = N + 1;


	for (int bit_id = 0; bit_id < cp.N; bit_id++)
	{
		int * cols_ids = new int[N + 1];
		for (int col_id = 0; col_id < N + 1; col_id++)
		{
			cols_ids[col_id] = 0;
		}
		fill_cols_h_x(cols_ids, cp.N, bit_id);

		crsMatrix * tmp_mtx_1 = new crsMatrix(N + 1, N + 1);

		for (int i = 0; i < N + 1; i++)
		{
			tmp_mtx_1->RowIndex[i] = i;
			tmp_mtx_1->Col[i] = cols_ids[i];
			tmp_mtx_1->Value[i].re = 1.0;
		}
		tmp_mtx_1->RowIndex[N + 1] = N + 1;


		crsMatrix * tmp_mtx_2 = new crsMatrix(*H_s_x);

		dcomplex beta;
		beta.re = 1.0;
		beta.im = 0.0;

		SparseMKLAdd(*tmp_mtx_1, beta, *tmp_mtx_2, *H_s_x);

		delete tmp_mtx_1;
		delete tmp_mtx_2;

		delete[] cols_ids;
	}

	if (rp.debug == 1)
	{
		string fn = rp.path + "H_s_x" + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, H_s_x, 16, false);
	}

	return H_s_x;
}
void init_h_s_x_vector(Model * m, RunParam &rp, ConfigParam &cp, MainData &md)
{
	crsMatrix * H_s_x = create_H_s_x_matrix(m, rp, cp, md), *res;
	m->H_s_x = H_s_x;
	int N = m->N;
	crsMatrix * h_s_x = m->h_s_x;

	to_F_basis(H_s_x, h_s_x);
}

void init_H0(Model * m)
{
	int N_mat = m->N_mat;
	crsMatrix * H0 = new crsMatrix();

	crsMatrix * subSum1 = new crsMatrix();

	dcomplex sum;
	sum.re = 1.0;
	sum.im = 0.0;

	SparseMKLAdd(*(m->H_J), sum, *(m->H_h_z), *subSum1);
	SparseMKLAdd(*subSum1, sum, *(m->H_h_x), *H0);

	m->H0 = H0;
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

void gamma_mult(crsMatrix * res, double gamma)
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

	crsMatrix * mat;

	mat = new crsMatrix((N + 1), (N + 1) / 2);
	int * RowIndex = mat->RowIndex;
	int * Col = mat->Col;
	dcomplex * Value = mat->Value;

	init_diss1_A1_rows(RowIndex, cp.N, diss_id, cp);
	init_diss1_A1_cols(Col, cp.N, diss_id, cp);
	init_diss1_A1_vals(Value, cp.N, diss_id, cp);

	if (rp.debug == 1)
	{
		string fn = rp.path + "A1_" + to_string(diss_id) + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, mat, 16, false);
	}

	return mat;
}
crsMatrix * create_A2_diss1_matrix(Model * m, int diss_id, RunParam &rp, ConfigParam &cp)
{
	int N = m->N;

	crsMatrix * mat;

	mat = new crsMatrix((N + 1), (N + 1) / 2);
	int * RowIndex = mat->RowIndex;
	int * Col = mat->Col;
	dcomplex * Value = mat->Value;

	init_diss1_A2_rows(RowIndex, cp.N, diss_id, cp);
	init_diss1_A2_cols(Col, cp.N, diss_id, cp);
	init_diss1_A2_vals(Value, cp.N, diss_id, cp);

	if (rp.debug == 1)
	{
		string fn = rp.path + "A2_" + to_string(diss_id) + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, mat, 16, false);
	}

	return mat;
}
crsMatrix * create_A1_diss2_matrix(Model * m, int diss_id, RunParam &rp, ConfigParam &cp)
{
	int N = m->N;

	crsMatrix * mat = new crsMatrix();;

	crsMatrix * sigma_x_1 = create_sigma_i_x_matrix(m, cp, diss_id);
	crsMatrix * sigma_x_2 = create_sigma_i_x_matrix(m, cp, diss_id + 1);

	crsMatrix * sigma_z_1 = create_sigma_i_z_matrix(m, cp, diss_id);
	crsMatrix * sigma_z_2 = create_sigma_i_z_matrix(m, cp, diss_id + 1);

	dcomplex sum;
	sum.re = 1.0;
	sum.im = 0.0;

	crsMatrix * sub_sum_1 = new crsMatrix();
	SparseMKLAdd(*(sigma_x_1), sum, *(sigma_x_2), *sub_sum_1);

	dcomplex sub;
	sub.re = -1.0;
	sub.im = 0.0;

	crsMatrix * sub_sum_2 = new crsMatrix();
	SparseMKLAdd(*(sigma_z_1), sub, *(sigma_z_2), *sub_sum_2);

	SparseMKLMult(*(sub_sum_1), *(sub_sum_2), *mat);

	delete sigma_x_1;
	delete sigma_x_2;
	delete sigma_z_1;
	delete sigma_z_2;
	delete sub_sum_1;
	delete sub_sum_2;

	if (rp.debug == 1)
	{
		string fn = rp.path + "A1_" + to_string(diss_id) + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, mat, 16, false);
	}

	return mat;
}
crsMatrix * create_A2_diss2_matrix(Model * m, int diss_id, RunParam &rp, ConfigParam &cp)
{
	int N = m->N;

	crsMatrix * mat;

	mat = new crsMatrix((N + 1), (N + 1));
	int * RowIndex = mat->RowIndex;
	int * Col = mat->Col;
	dcomplex * Value = mat->Value;

	for (int state_id_1 = 0; state_id_1 < (N + 1); state_id_1++)
	{
		Col[state_id_1] = state_id_1;
		RowIndex[state_id_1] = state_id_1;
		Value[state_id_1].re = 0.0;
	}
	RowIndex[(N + 1)] = (N + 1);

	if (rp.debug == 1)
	{
		string fn = rp.path + "A2_" + to_string(diss_id) + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, mat, 16, false);
	}

	return mat;
}
crsMatrix * create_A1_diss0_matrix(Model * m, int diss_id, RunParam &rp, ConfigParam &cp)
{
	int N = m->N;

	crsMatrix * mat;

	mat = new crsMatrix((N + 1), (N + 1));
	int * RowIndex = mat->RowIndex;
	int * Col = mat->Col;
	dcomplex * Value = mat->Value;

	for (int state_id_1 = 0; state_id_1 < (N + 1); state_id_1++)
	{
		Col[state_id_1] = state_id_1;
		RowIndex[state_id_1] = state_id_1;

		if (bit_at(state_id_1, diss_id) == 1)
		{
			Value[state_id_1].re = 1.0;
		}
		else
		{
			Value[state_id_1].re = 0.0;
		}
	}
	RowIndex[(N + 1)] = (N + 1);

	if (rp.debug == 1)
	{
		string fn = rp.path + "A2_" + to_string(diss_id) + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, mat, 16, false);
	}

	return mat;
}
crsMatrix * create_A2_diss0_matrix(Model * m, int diss_id, RunParam &rp, ConfigParam &cp)
{
	int N = m->N;

	crsMatrix * mat;

	mat = new crsMatrix((N + 1), (N + 1));
	int * RowIndex = mat->RowIndex;
	int * Col = mat->Col;
	dcomplex * Value = mat->Value;

	for (int state_id_1 = 0; state_id_1 < (N + 1); state_id_1++)
	{
		Col[state_id_1] = state_id_1;
		RowIndex[state_id_1] = state_id_1;
		Value[state_id_1].re = 0.0;
	}
	RowIndex[(N + 1)] = (N + 1);

	if (rp.debug == 1)
	{
		string fn = rp.path + "A2_" + to_string(diss_id) + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, mat, 16, false);
	}

	return mat;
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

	int N_mat = m->N_mat;

	crsMatrix * result_a_matrix = NULL;

	double time = omp_get_wtime();
	double init_time = time;

	crsMatrix * A1 = create_A1_diss1_matrix(m, 0, rp, cp);
	crsMatrix * A2 = create_A2_diss1_matrix(m, 0, rp, cp);

	crsMatrix * a1_mat = new crsMatrix(N_mat, N_mat);
	crsMatrix * a2_mat = new crsMatrix(N_mat, N_mat);

	to_F_basis(A1, a1_mat);
	to_F_basis(A2, a2_mat);

	vector<pair<int, dcomplex> > * a_std = create_a_std_matrix(a1_mat, a2_mat, N_mat);

	delete a1_mat;
	delete a2_mat;

	delete A1;
	delete A2;

	crsMatrix * a1_i_a2_mat = NULL;
	a1_i_a2_mat = stdToCrs(a_std, N_mat);

	delete[] a_std;

	gamma_mult(a1_i_a2_mat, cp.g);
	result_a_matrix = new crsMatrix(*a1_i_a2_mat);

	delete a1_i_a2_mat;

	time = omp_get_wtime() - init_time;
	cout << "time of a_" << "0 : " << time << endl << endl;

	int diss_num = cp.N - 1;

	if (N > 1)
	{
		for (int dissId = 1; dissId < diss_num; dissId++)
		{
			crsMatrix * A1 = create_A1_diss1_matrix(m, dissId, rp, cp);
			crsMatrix * A2 = create_A2_diss1_matrix(m, dissId, rp, cp);

			crsMatrix * a1_mat = new crsMatrix(N_mat, N_mat);
			crsMatrix * a2_mat = new crsMatrix(N_mat, N_mat);

			to_F_basis(A1, a1_mat);
			to_F_basis(A2, a2_mat);

			vector<pair<int, dcomplex> > * a_std = create_a_std_matrix(a1_mat, a2_mat, N_mat);

			delete a1_mat;
			delete a2_mat;

			delete A1;
			delete A2;

			crsMatrix * a1_i_a2_mat = NULL;
			a1_i_a2_mat = stdToCrs(a_std, N_mat);

			delete[] a_std;

			gamma_mult(a1_i_a2_mat, cp.g);

			crsMatrix * a_mat_tmp = new crsMatrix(*result_a_matrix);
			delete result_a_matrix;
			result_a_matrix = new crsMatrix();

			dcomplex beta;
			beta.re = 1.0;
			beta.im = 0.0;

			SparseMKLAdd(*a_mat_tmp, beta, *a1_i_a2_mat, *result_a_matrix);

			delete a_mat_tmp;
			delete a1_i_a2_mat;

			time = omp_get_wtime() - init_time;
			cout << "time of a_" << dissId << " : " << time << endl << endl;
		}
	}

	if (rp.debug == 1)
	{
		string fn = rp.path + "a" + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, result_a_matrix, 16, false);
	}

	m->l_mat = new crsMatrix(*result_a_matrix);
	delete result_a_matrix;
}
void init_diss_2(Model * m, RunParam &rp, ConfigParam &cp, MainData &md)
{
	int N = m->N;

	int N_mat = m->N_mat;

	crsMatrix * result_a_matrix = NULL;

	double time = omp_get_wtime();
	double init_time = time;

	crsMatrix * A1 = create_A1_diss2_matrix(m, 0, rp, cp);
	crsMatrix * A2 = create_A2_diss2_matrix(m, 0, rp, cp);

	crsMatrix * a1_mat = new crsMatrix(N_mat, N_mat);
	crsMatrix * a2_mat = new crsMatrix(N_mat, N_mat);

	to_F_basis(A1, a1_mat);
	to_F_basis(A2, a2_mat);

	vector<pair<int, dcomplex> > * a_std = create_a_std_matrix(a1_mat, a2_mat, N_mat);

	delete a1_mat;
	delete a2_mat;

	delete A1;
	delete A2;

	crsMatrix * a1_i_a2_mat = NULL;
	a1_i_a2_mat = stdToCrs(a_std, N_mat);

	delete[] a_std;

	gamma_mult(a1_i_a2_mat, cp.g);

	result_a_matrix = new crsMatrix(*a1_i_a2_mat);

	delete a1_i_a2_mat;

	time = omp_get_wtime() - init_time;
	cout << "time of a_" << "0 : " << time << endl << endl;

	int diss_num = cp.N - 1;

	if (N > 1)
	{
		for (int dissId = 1; dissId < diss_num; dissId++)
		{
			crsMatrix * A1 = create_A1_diss2_matrix(m, dissId, rp, cp);
			crsMatrix * A2 = create_A2_diss2_matrix(m, dissId, rp, cp);

			crsMatrix * a1_mat = new crsMatrix(N_mat, N_mat);
			crsMatrix * a2_mat = new crsMatrix(N_mat, N_mat);

			to_F_basis(A1, a1_mat);
			to_F_basis(A2, a2_mat);

			vector<pair<int, dcomplex> > * a_std = create_a_std_matrix(a1_mat, a2_mat, N_mat);

			delete a1_mat;
			delete a2_mat;

			delete A1;
			delete A2;

			crsMatrix * a1_i_a2_mat = NULL;
			a1_i_a2_mat = stdToCrs(a_std, N_mat);

			delete[] a_std;

			gamma_mult(a1_i_a2_mat, cp.g);

			crsMatrix * a_mat_tmp = new crsMatrix(*result_a_matrix);
			delete result_a_matrix;
			result_a_matrix = new crsMatrix();

			dcomplex beta;
			beta.re = 1.0;
			beta.im = 0.0;

			SparseMKLAdd(*a_mat_tmp, beta, *a1_i_a2_mat, *result_a_matrix);

			delete a_mat_tmp;
			delete a1_i_a2_mat;

			time = omp_get_wtime() - init_time;
			cout << "time of a_" << dissId << " : " << time << endl << endl;
		}
	}

	if (rp.debug == 1)
	{
		string fn = rp.path + "a" + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, result_a_matrix, 16, false);
	}

	m->l_mat = new crsMatrix(*result_a_matrix);
	delete result_a_matrix;
}
void init_diss_3(Model * m, RunParam &rp, ConfigParam &cp, MainData &md)
{
	int N = m->N;

	int N_mat = m->N_mat;

	crsMatrix * result_a_matrix = NULL;

	double time = omp_get_wtime();
	double init_time = time;

	crsMatrix * A1 = create_A1_diss1_matrix(m, 0, rp, cp);
	crsMatrix * A2 = create_A2_diss1_matrix(m, 0, rp, cp);

	crsMatrix * a1_mat = new crsMatrix(N_mat, N_mat);
	crsMatrix * a2_mat = new crsMatrix(N_mat, N_mat);

	to_F_basis(A1, a1_mat);
	to_F_basis(A2, a2_mat);

	vector<pair<int, dcomplex> > * a_std = create_a_std_matrix(a1_mat, a2_mat, N_mat);

	delete a1_mat;
	delete a2_mat;

	delete A1;
	delete A2;

	crsMatrix * a1_i_a2_mat = NULL;
	a1_i_a2_mat = stdToCrs(a_std, N_mat);

	delete[] a_std;

	gamma_mult(a1_i_a2_mat, cp.g);
	result_a_matrix = new crsMatrix(*a1_i_a2_mat);

	delete a1_i_a2_mat;

	time = omp_get_wtime() - init_time;
	cout << "time of a_" << "0 : " << time << endl << endl;

	int diss_num = cp.N - 1;

	if (N > 1)
	{
		for (int dissId = 1; dissId < diss_num; dissId++)
		{
			crsMatrix * A1 = create_A1_diss1_matrix(m, dissId, rp, cp);
			crsMatrix * A2 = create_A2_diss1_matrix(m, dissId, rp, cp);

			crsMatrix * a1_mat = new crsMatrix(N_mat, N_mat);
			crsMatrix * a2_mat = new crsMatrix(N_mat, N_mat);

			to_F_basis(A1, a1_mat);
			to_F_basis(A2, a2_mat);
			vector<pair<int, dcomplex> > * a_std = create_a_std_matrix(a1_mat, a2_mat, N_mat);

			delete a1_mat;
			delete a2_mat;

			delete A1;
			delete A2;

			crsMatrix * a1_i_a2_mat = NULL;
			a1_i_a2_mat = stdToCrs(a_std, N_mat);

			delete[] a_std;

			gamma_mult(a1_i_a2_mat, cp.g);

			crsMatrix * a_mat_tmp = new crsMatrix(*result_a_matrix);
			delete result_a_matrix;
			result_a_matrix = new crsMatrix();

			dcomplex beta;
			beta.re = 1.0;
			beta.im = 0.0;

			SparseMKLAdd(*a_mat_tmp, beta, *a1_i_a2_mat, *result_a_matrix);

			delete a_mat_tmp;
			delete a1_i_a2_mat;

			time = omp_get_wtime() - init_time;
			cout << "time of a_" << dissId << " : " << time << endl << endl;

		}

		int diss_num = cp.N;

		for (int dissId = 0; dissId < diss_num; dissId++)
		{
			crsMatrix * A1 = create_A1_diss0_matrix(m, dissId, rp, cp);
			crsMatrix * A2 = create_A2_diss0_matrix(m, dissId, rp, cp);

			crsMatrix * a1_mat = new crsMatrix(N_mat, N_mat);
			crsMatrix * a2_mat = new crsMatrix(N_mat, N_mat);

			to_F_basis(A1, a1_mat);
			to_F_basis(A2, a2_mat);

			vector<pair<int, dcomplex> > * a_std = create_a_std_matrix(a1_mat, a2_mat, N_mat);

			delete a1_mat;
			delete a2_mat;

			delete A1;
			delete A2;

			crsMatrix * a1_i_a2_mat = NULL;
			a1_i_a2_mat = stdToCrs(a_std, N_mat);

			delete[] a_std;

			gamma_mult(a1_i_a2_mat, cp.g_add);

			crsMatrix * a_mat_tmp = new crsMatrix(*result_a_matrix);
			delete result_a_matrix;
			result_a_matrix = new crsMatrix();

			dcomplex beta;
			beta.re = 1.0;
			beta.im = 0.0;

			SparseMKLAdd(*a_mat_tmp, beta, *a1_i_a2_mat, *result_a_matrix);

			delete a_mat_tmp;
			delete a1_i_a2_mat;

			time = omp_get_wtime() - init_time;
			cout << "time of a_" << dissId << " : " << time << endl << endl;

		}
	}

	if (rp.debug == 1)
	{
		string fn = rp.path + "a" + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, result_a_matrix, 16, false);
	}

	m->l_mat = new crsMatrix(*result_a_matrix);
	delete result_a_matrix;
}
void init_diss_4(Model * m, RunParam &rp, ConfigParam &cp, MainData &md)
{
	int N = m->N;

	int N_mat = m->N_mat;

	crsMatrix * result_a_matrix = NULL;

	double time = omp_get_wtime();
	double init_time = time;

	crsMatrix * A1 = create_A1_diss2_matrix(m, 0, rp, cp);
	crsMatrix * A2 = create_A2_diss2_matrix(m, 0, rp, cp);

	crsMatrix * a1_mat = new crsMatrix(N_mat, N_mat);
	crsMatrix * a2_mat = new crsMatrix(N_mat, N_mat);

	to_F_basis(A1, a1_mat);
	to_F_basis(A2, a2_mat);

	vector<pair<int, dcomplex> > * a_std = create_a_std_matrix(a1_mat, a2_mat, N_mat);

	delete a1_mat;
	delete a2_mat;

	delete A1;
	delete A2;

	crsMatrix * a1_i_a2_mat = NULL;
	a1_i_a2_mat = stdToCrs(a_std, N_mat);

	delete[] a_std;

	gamma_mult(a1_i_a2_mat, cp.g);
	result_a_matrix = new crsMatrix(*a1_i_a2_mat);

	delete a1_i_a2_mat;

	time = omp_get_wtime() - init_time;
	cout << "time of a_" << "0 : " << time << endl << endl;

	int diss_num = cp.N - 1;

	if (N > 1)
	{
		for (int dissId = 1; dissId < diss_num; dissId++)
		{
			crsMatrix * A1 = create_A1_diss2_matrix(m, dissId, rp, cp);
			crsMatrix * A2 = create_A2_diss2_matrix(m, dissId, rp, cp);

			crsMatrix * a1_mat = new crsMatrix(N_mat, N_mat);
			crsMatrix * a2_mat = new crsMatrix(N_mat, N_mat);

			to_F_basis(A1, a1_mat);
			to_F_basis(A2, a2_mat);

			vector<pair<int, dcomplex> > * a_std = create_a_std_matrix(a1_mat, a2_mat, N_mat);

			delete a1_mat;
			delete a2_mat;

			delete A1;
			delete A2;

			crsMatrix * a1_i_a2_mat = NULL;
			a1_i_a2_mat = stdToCrs(a_std, N_mat);

			delete[] a_std;

			gamma_mult(a1_i_a2_mat, cp.g);

			crsMatrix * a_mat_tmp = new crsMatrix(*result_a_matrix);
			delete result_a_matrix;
			result_a_matrix = new crsMatrix();

			dcomplex beta;
			beta.re = 1.0;
			beta.im = 0.0;

			SparseMKLAdd(*a_mat_tmp, beta, *a1_i_a2_mat, *result_a_matrix);

			delete a_mat_tmp;
			delete a1_i_a2_mat;

			time = omp_get_wtime() - init_time;
			cout << "time of a_" << dissId << " : " << time << endl << endl;

		}

		int diss_num = cp.N;

		for (int dissId = 0; dissId < diss_num; dissId++)
		{
			crsMatrix * A1 = create_A1_diss0_matrix(m, dissId, rp, cp);
			crsMatrix * A2 = create_A2_diss0_matrix(m, dissId, rp, cp);

			crsMatrix * a1_mat = new crsMatrix(N_mat, N_mat);
			crsMatrix * a2_mat = new crsMatrix(N_mat, N_mat);

			to_F_basis(A1, a1_mat);
			to_F_basis(A2, a2_mat);

			vector<pair<int, dcomplex> > * a_std = create_a_std_matrix(a1_mat, a2_mat, N_mat);

			delete a1_mat;
			delete a2_mat;

			delete A1;
			delete A2;

			crsMatrix * a1_i_a2_mat = NULL;
			a1_i_a2_mat = stdToCrs(a_std, N_mat);

			delete[] a_std;

			gamma_mult(a1_i_a2_mat, cp.g_add);

			crsMatrix * a_mat_tmp = new crsMatrix(*result_a_matrix);
			delete result_a_matrix;
			result_a_matrix = new crsMatrix();

			dcomplex beta;
			beta.re = 1.0;
			beta.im = 0.0;

			SparseMKLAdd(*a_mat_tmp, beta, *a1_i_a2_mat, *result_a_matrix);

			delete a_mat_tmp;
			delete a1_i_a2_mat;

			time = omp_get_wtime() - init_time;
			cout << "time of a_" << dissId << " : " << time << endl << endl;

		}
	}

	if (rp.debug == 1)
	{
		string fn = rp.path + "a" + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, result_a_matrix, 16, false);
	}

	m->l_mat = new crsMatrix(*result_a_matrix);
	delete result_a_matrix;
}

void calc_Q_J_s(Model * m)
{
	int N_mat = m->N_mat;
	crsMatrix * h_J = m->h_J;

	crsMatrix * Q_J_s;
	calc_CooQs(N_mat, m, m->f_ijk, h_J, Q_J_s);

	m->Q_J_s = Q_J_s;
}
void calc_Q_h_z_s(Model * m)
{
	int N_mat = m->N_mat;
	crsMatrix * h_h_z = m->h_h_z;

	crsMatrix * Q_h_z_s;
	calc_CooQs(N_mat, m, m->f_ijk, h_h_z, Q_h_z_s);

	m->Q_h_z_s = Q_h_z_s;
}
void calc_Q_h_x_s(Model * m)
{
	int N_mat = m->N_mat;
	crsMatrix * h_h_x = m->h_h_x;

	crsMatrix * Q_h_x_s;
	calc_CooQs(N_mat, m, m->f_ijk, h_h_x, Q_h_x_s);

	m->Q_h_x_s = Q_h_x_s;
}
void calc_Q_s_x_s(Model * m)
{
	int N_mat = m->N_mat;
	crsMatrix * h_s_x = m->h_s_x;

	crsMatrix * Q_s_x_s;
	calc_CooQs(N_mat, m, m->f_ijk, h_s_x, Q_s_x_s);

	m->Q_s_x_s = Q_s_x_s;
}

void calcKs(Model *m)
{
	int N = m->N;
	int N_mat = m->N_mat;
	dcomplex  *Ks = m->Ks;

	crsMatrix * l_mat = m->l_mat;

	m->f_ijk = new Tensor_Coordinates;
	fijk_coord(m->f_ijk, m->N + 1);
	Tensor_Coordinates * f_ijk = m->f_ijk;

	for (int i = 0; i < N_mat; i++)
	{
		Ks[i].re = 0.0;
		Ks[i].im = 0.0;
	}

	ulli uN = l_mat->N;
	ulli uN2 = uN * uN;
	for (ulli i = 0; i < f_ijk->k; i++)
	{
		int ii = f_ijk->coord1[i];
		int ji = f_ijk->coord2[i];
		int ki = f_ijk->coord3[i];
		f_ijk->hash[i] = ji * uN + ii;
	}
	sort_matrix(f_ijk);

	dcomplex * hash = new dcomplex[l_mat->N];

	for (ulli i = 0; i < f_ijk->k; i++)
	{
		int sji = f_ijk->coord2[i];

		memset(hash, 0, sizeof(dcomplex) * l_mat->N);
		int start = l_mat->RowIndex[sji];
		int finish = l_mat->RowIndex[sji + 1];
		for (int ind = start; ind < finish; ind++)
		{
			hash[l_mat->Col[ind]] = l_mat->Value[ind];
		}

		while ((i < f_ijk->k) && (f_ijk->coord2[i] == sji))
		{
			int ii = f_ijk->coord1[i];
			int ki = f_ijk->coord3[i];
			int ji = f_ijk->coord2[i];

			dcomplex v1 = hash[ii];
			dcomplex v2 = f_ijk->data[i];
			dcomplex tmp;
			tmp.re = v1.re * v2.re - v1.im * v2.im;
			tmp.im = v1.re * v2.im + v1.im * v2.re;

			Ks[ki].re += tmp.re;
			Ks[ki].im += tmp.im;
			i++;
		}
		i--;
	}

	dcomplex val;
	for (int i = 0; i < N_mat; i++)
	{
		Ks[i].re *= (1.0) / (N + 1);
		Ks[i].im *= (-1.0) / (N + 1);
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
	crsMatrix * Rs;

	Tensor_Coordinates * subF_ijk;
	Tensor_Coordinates * subD_ijk;
	Tensor_Coordinates * subZ_ijk;

	ulli cnt;

	unsigned int *indexZ = new unsigned int[N_mat + 1];
	unsigned int *indexF = new unsigned int[N_mat + 1];

	subD_ijk = new Tensor_Coordinates;
	dijk_coord(subD_ijk, m->N + 1);


	subF_ijk = m->f_ijk;

	ulli uN = m->l_mat->N;
	ulli uN2 = uN * uN;
	for (ulli i = 0; i < subF_ijk->k; i++)
	{
		int ii = subF_ijk->coord1[i];
		int ji = subF_ijk->coord2[i];
		int ki = subF_ijk->coord3[i];
		subF_ijk->hash[i] = ji * uN2 + ii * uN + ki;
	}
	for (ulli i = 0; i < subD_ijk->k; i++)
	{
		int ii = subD_ijk->coord1[i];
		int ji = subD_ijk->coord2[i];
		int ki = subD_ijk->coord3[i];
		subD_ijk->hash[i] = ji * uN2 + ii * uN + ki;
	}
	sort_matrix(subF_ijk);
	sort_matrix(subD_ijk);

	calcZ_ijk(subF_ijk, subD_ijk, subZ_ijk);

	vector<map<int, dcomplex> > mat(N_mat + 1);

	createIndex(subZ_ijk, N_mat, indexZ);
	createIndex(subF_ijk, N_mat, indexF);

	multTmpRsSTD(N_mat, indexF, indexZ, subZ_ijk, subF_ijk, m->l_mat, false, mat);

	for (ulli i = 0; i < subZ_ijk->k; i++)
	{
		subZ_ijk->data[i].im = -subZ_ijk->data[i].im;
	}

	multTmpRsSTD(N_mat, indexF, indexZ, subZ_ijk, subF_ijk, m->l_mat, true, mat);

	crsMatrix * Rs_std = stdToCrs(mat, N_mat);

	free_matrix(subD_ijk);
	free_matrix(subZ_ijk);
	delete subD_ijk;
	delete subZ_ijk;

	Rs = Rs_std;

	delete[] indexZ;
	delete[] indexF;


	double mm = -1.0 / 4.0;
	for (int i = 0; i < Rs->NZ; i++)
	{
		Rs->Value[i].re *= mm;
		Rs->Value[i].im *= mm;
	}
	m->Rs = Rs;
}
void calc_G_0_s(Model *m)
{
	int N_mat = m->N_mat;
	crsMatrix * G_0_s = new crsMatrix();

	crsMatrix * subSum1 = new crsMatrix();
	crsMatrix * subSum2 = new crsMatrix();

	dcomplex sum;
	sum.re = 1.0;
	sum.im = 0.0;

	SparseMKLAdd(*(m->Q_J_s), sum, *(m->Q_h_z_s), *subSum1);
	SparseMKLAdd(*(m->Q_h_x_s), sum, *(m->Rs), *subSum2);

	SparseMKLAdd(*(subSum1), sum, *(subSum2), *G_0_s);

	m->G_0_s = G_0_s;

	delete subSum1;
	delete subSum2;
}
void calc_G_1_s(Model *m)
{
	int N_mat = m->N_mat;
	crsMatrix * G_1_s = new crsMatrix();

	dcomplex sum;
	sum.re = 1.0;
	sum.im = 0.0;

	SparseMKLAdd(*(m->Q_s_x_s), sum, *(m->Rs), *G_1_s);

	m->G_1_s = G_1_s;
}

void multMatVec(crsMatrix *mat, double * x, double * res)
{
	char trans = 'n';
	double *value = (double *)(mat->Value);
	mkl_dcsrgemv(&trans, &(mat->N), value, mat->RowIndex, mat->Col, x, res);
}
void multMatVec_complex(crsMatrix *mat, dcomplex * x, dcomplex * res)
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
			res[i].re += v1.re * v2.re;
			res[i].im = 0.0;
		}
	}
}

void calcVectValue_t0(double h, Model * m, double *x, double * res, double * tmp)
{
	int i;
	int N_mat = m->N_mat;
	crsMatrix * G_0_s = m->G_0_s;
	double  * Ks = (double *)m->Ks;

	multMatVec(G_0_s, x, tmp);

	for (i = 0; i < N_mat; i++)
	{
		res[i] = (tmp[i] - Ks[i]) * h;
	}
}
void calcVectValue_t1(double h, Model * m, double *x, double * res, double * tmp)
{
	int i;
	int N_mat = m->N_mat;
	crsMatrix * G_1_s = m->G_1_s;
	double  * Ks = (double *)m->Ks;

	multMatVec(G_1_s, x, tmp);

	for (i = 0; i < N_mat; i++)
	{
		res[i] = (tmp[i] - Ks[i]) * h;
	}
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
	else if (cp.int_ist == 2)
	{
		double * init_vec = new double[md.size];

		for (int s_id = 0; s_id < md.size; s_id++)
		{
			init_vec[s_id] = 1.0;
			for (int bit_id = 0; bit_id < cp.N; bit_id++)
			{
				if (bit_at(s_id, bit_id) == 0)
				{
					init_vec[s_id] *= cos(PI * 0.125);
				}
				else
				{
					init_vec[s_id] *= sin(PI * 0.125);
				}
			}
		}

		for (int s_id_1 = 0; s_id_1 < md.size; s_id_1++)
		{
			for (int s_id_2 = 0; s_id_2 < md.size; s_id_2++)
			{
				mtx[s_id_1 * md.size + s_id_2].re = init_vec[s_id_1] * init_vec[s_id_2];
				mtx[s_id_1 * md.size + s_id_2].im = 0.0;
			}
		}

		delete[] init_vec;
	}
	else if (cp.int_ist == 3)
	{
		int random = 0;
		VSLStreamStatePtr stream;
		vslNewStream(&stream, VSL_BRNG_MCG31, 77778888);
		vslLeapfrogStream(stream, cp.seed, cp.max_num_seeds);
		viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, 1, &random, 0, pow(2, cp.N));

		cout << "Random state: " << random << endl;

		string fn = rp.path + "init_state" + file_name_suffix(cp, 4);
		save_int_data(fn, &random, 1, false);

		mtx[random * md.size + random].re = 1.0;
		mtx[random * md.size + random].im = 0.0;
	}
	else
	{
		stringstream msg;
		msg << "wrong int_ist value: " << cp.int_ist << endl;
		Error(msg.str());
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

void before(Model *m)
{
	complex_to_real(m->G_0_s->Value, m->G_0_s->NZ);
	complex_to_real(m->G_1_s->Value, m->G_1_s->NZ);

	complex_to_real(m->Q_J_s->Value, m->Q_J_s->NZ);
	complex_to_real(m->Q_h_z_s->Value, m->Q_h_z_s->NZ);
	complex_to_real(m->Q_h_x_s->Value, m->Q_h_x_s->NZ);
	complex_to_real(m->Q_s_x_s->Value, m->Q_s_x_s->NZ);

	complex_to_real(m->Ks, m->N_mat);
	complex_to_real(m->RhoF, m->N_mat);

	toOneBase(*(m->G_0_s));
	toOneBase(*(m->G_1_s));

	toOneBase(*(m->Q_J_s));
	toOneBase(*(m->Q_h_z_s));
	toOneBase(*(m->Q_h_x_s));
	toOneBase(*(m->Q_s_x_s));
}

void after(Model *m)
{
	toZeroBase(*(m->G_0_s));
	toZeroBase(*(m->G_1_s));

	toZeroBase(*(m->Q_J_s));
	toZeroBase(*(m->Q_h_z_s));
	toZeroBase(*(m->Q_h_x_s));
	toZeroBase(*(m->Q_s_x_s));

	real_to_complex(m->G_0_s->Value, m->G_0_s->NZ);
	real_to_complex(m->G_1_s->Value, m->G_1_s->NZ);

	real_to_complex(m->Q_J_s->Value, m->Q_J_s->NZ);
	real_to_complex(m->Q_h_z_s->Value, m->Q_h_z_s->NZ);
	real_to_complex(m->Q_h_x_s->Value, m->Q_h_x_s->NZ);
	real_to_complex(m->Q_s_x_s->Value, m->Q_s_x_s->NZ);

	real_to_complex(m->Ks, m->N_mat);
	real_to_complex(m->RhoF, m->N_mat);
}

void calcODE_std(Model *m, RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd)
{
	int i;
	int N_mat = m->N_mat;
	double * RhoF = (double *)(m->RhoF);
	double * prevRhoF = (double *)(m->prevRhoF);

	double eps = 1.0e-10;

	for (i = 0; i < N_mat; i++)
	{
		prevRhoF[i] = RhoF[i];
	}

	double * k1 = pd.k1;
	double * k2 = pd.k2;
	double * k3 = pd.k3;
	double * k4 = pd.k4;
	double * val = pd.val;
	double * tmp = pd.tmp;

	double step_t_0 = md.step_t_0;
	double step_t_1 = md.step_t_1;

	double curr_time = 0.0;
	int dump_id = 0;

	calcRho(m);

	init_sign_sigma_z_start(m, cp, md, pd);

	characteristics_std(m, rp, cp, md, pd, dump_id);
	dump_id++;
	
	before(m);

	for (int period = 1; period <= cp.num_periods; period++)
	{
		for (int t_0_step_id = 0; t_0_step_id < cp.num_steps_t_0; t_0_step_id++)
		{
			calcVectValue_t0(step_t_0, m, RhoF, k1, tmp);
			for (i = 0; i < N_mat; i++)
			{
				val[i] = RhoF[i] + k1[i] / 2.0;
			}
			calcVectValue_t0(step_t_0, m, val, k2, tmp);
			for (i = 0; i < N_mat; i++)
			{
				val[i] = RhoF[i] + k2[i] / 2.0;
			}
			calcVectValue_t0(step_t_0, m, val, k3, tmp);
			for (i = 0; i < N_mat; i++)
			{
				val[i] = RhoF[i] + k3[i];
			}
			calcVectValue_t0(step_t_0, m, val, k4, tmp);

			for (i = 0; i < N_mat; i++)
			{
				RhoF[i] += (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
			}
		}

		for (int t_1_step_id = 0; t_1_step_id < cp.num_steps_t_1; t_1_step_id++)
		{
			calcVectValue_t1(step_t_1, m, RhoF, k1, tmp);
			for (i = 0; i < N_mat; i++)
			{
				val[i] = RhoF[i] + k1[i] / 2.0;
			}
			calcVectValue_t1(step_t_1, m, val, k2, tmp);
			for (i = 0; i < N_mat; i++)
			{
				val[i] = RhoF[i] + k2[i] / 2.0;
			}
			calcVectValue_t1(step_t_1, m, val, k3, tmp);
			for (i = 0; i < N_mat; i++)
			{
				val[i] = RhoF[i] + k3[i];
			}
			calcVectValue_t1(step_t_1, m, val, k4, tmp);

			for (i = 0; i < N_mat; i++)
			{
				RhoF[i] += (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
			}
		}

		curr_time = period * md.T;

		if (period == pd.dump_periods[dump_id])
		{

			if (rp.ipp == 1)
			{
				cout << endl << "Dump period: " << period << endl;
			}

			after(m);
			calcRho(m);
			characteristics_std(m, rp, cp, md, pd, dump_id);
			before(m);

			dump_id++;
		}
	}

	after(m);
}
void calcODE_rate(Model *m, RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd)
{
	int i;
	int N_mat = m->N_mat;
	double * RhoF = (double *)(m->RhoF);
	double * prevRhoF = (double *)(m->prevRhoF);

	double eps = 1.0e-10;

	for (i = 0; i < N_mat; i++)
	{
		prevRhoF[i] = RhoF[i];
	}

	double * k1 = pd.k1;
	double * k2 = pd.k2;
	double * k3 = pd.k3;
	double * k4 = pd.k4;
	double * val = pd.val;
	double * tmp = pd.tmp;

	double step_t_0 = md.step_t_0;
	double step_t_1 = md.step_t_1;

	double curr_time = 0.0;
	int dump_id = 0;

	calcRho(m);

	init_sign_sigma_z_start(m, cp, md, pd);

	characteristics_std(m, rp, cp, md, pd, dump_id);
	dump_id++;

	before(m);

	int period = 1;

	while (period <= cp.num_periods)
	{
		for (int t_0_step_id = 0; t_0_step_id < cp.num_steps_t_0; t_0_step_id++)
		{
			calcVectValue_t0(step_t_0, m, RhoF, k1, tmp);
			for (i = 0; i < N_mat; i++)
			{
				val[i] = RhoF[i] + k1[i] / 2.0;
			}
			calcVectValue_t0(step_t_0, m, val, k2, tmp);
			for (i = 0; i < N_mat; i++)
			{
				val[i] = RhoF[i] + k2[i] / 2.0;
			}
			calcVectValue_t0(step_t_0, m, val, k3, tmp);
			for (i = 0; i < N_mat; i++)
			{
				val[i] = RhoF[i] + k3[i];
			}
			calcVectValue_t0(step_t_0, m, val, k4, tmp);

			for (i = 0; i < N_mat; i++)
			{
				RhoF[i] += (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
			}

			after(m);
			calcRho(m);
			if (period - 1 < cp.rate_num_periods)
			{
				characteristics_rate(m, rp, cp, md, pd, 0, period);
			}
			else
			{
				characteristics_rate(m, rp, cp, md, pd, 1, period);
			}
			before(m);
		}

		for (int t_1_step_id = 0; t_1_step_id < cp.num_steps_t_1; t_1_step_id++)
		{
			calcVectValue_t1(step_t_1, m, RhoF, k1, tmp);
			for (i = 0; i < N_mat; i++)
			{
				val[i] = RhoF[i] + k1[i] / 2.0;
			}
			calcVectValue_t1(step_t_1, m, val, k2, tmp);
			for (i = 0; i < N_mat; i++)
			{
				val[i] = RhoF[i] + k2[i] / 2.0;
			}
			calcVectValue_t1(step_t_1, m, val, k3, tmp);
			for (i = 0; i < N_mat; i++)
			{
				val[i] = RhoF[i] + k3[i];
			}
			calcVectValue_t1(step_t_1, m, val, k4, tmp);

			for (i = 0; i < N_mat; i++)
			{
				RhoF[i] += (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
			}

			after(m);
			calcRho(m);
			if (period - 1 < cp.rate_num_periods)
			{
				characteristics_rate(m, rp, cp, md, pd, 0, period);
			}
			else
			{
				characteristics_rate(m, rp, cp, md, pd, 1, period);
			}
			before(m);
		}

		if((period - 1) % cp.rate_num_periods == 0)
		{ 
			characteristics_rate(m, rp, cp, md, pd, 2, period);
		}

		curr_time = period * md.T;

		if (period == pd.dump_periods[dump_id])
		{
			after(m);
			calcRho(m);
			characteristics_std(m, rp, cp, md, pd, dump_id);
			before(m);

			dump_id++;
		}

		if (rp.ipp == 1)
		{
			cout << endl << "Period: " << period << endl;
		}

		period++;
	}

	after(m);
}
void calcODE_deep(Model *m, RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd)
{
	string fn;

	int dupm_step_t_0 = cp.num_steps_t_0 / cp.int_dn;
	int dupm_step_t_1 = cp.num_steps_t_1 / cp.int_dn;

	int i;
	int N_mat = m->N_mat;
	double * RhoF = (double *)(m->RhoF);
	double * prevRhoF = (double *)(m->prevRhoF);

	double eps = 1.0e-10;

	for (i = 0; i < N_mat; i++)
	{
		prevRhoF[i] = RhoF[i];
	}

	double * k1 = pd.k1;
	double * k2 = pd.k2;
	double * k3 = pd.k3;
	double * k4 = pd.k4;
	double * val = pd.val;
	double * tmp = pd.tmp;

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

	before(m);

	for (int period = 1; period <= cp.num_periods; period++)
	{
		for (int t_0_step_id = 0; t_0_step_id < cp.num_steps_t_0; t_0_step_id++)
		{
			calcVectValue_t0(step_t_0, m, RhoF, k1, tmp);
			for (i = 0; i < N_mat; i++)
			{
				val[i] = RhoF[i] + k1[i] / 2.0;
			}
			calcVectValue_t0(step_t_0, m, val, k2, tmp);
			for (i = 0; i < N_mat; i++)
			{
				val[i] = RhoF[i] + k2[i] / 2.0;
			}
			calcVectValue_t0(step_t_0, m, val, k3, tmp);
			for (i = 0; i < N_mat; i++)
			{
				val[i] = RhoF[i] + k3[i];
			}
			calcVectValue_t0(step_t_0, m, val, k4, tmp);

			for (i = 0; i < N_mat; i++)
			{
				RhoF[i] += (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
			}

			if (rp.debug == 1)
			{
				after(m);
				calcRho(m);
				before(m);

				fn = rp.path + "rho" + "_" + to_string(dump_id) + file_name_suffix(cp, 4);
				save_sparse_complex_mtx(fn, m->Rho, 16, false);
			}

			if (t_0_step_id % dupm_step_t_0 == 0)
			{
				after(m);
				calcRho(m);
				
				characteristics_deep(m, rp, cp, md, pd, dump_id);
				before(m);
				dump_id++;
			}		
		}

		for (int t_1_step_id = 0; t_1_step_id < cp.num_steps_t_1; t_1_step_id++)
		{
			calcVectValue_t1(step_t_1, m, RhoF, k1, tmp);
			for (i = 0; i < N_mat; i++)
			{
				val[i] = RhoF[i] + k1[i] / 2.0;
			}
			calcVectValue_t1(step_t_1, m, val, k2, tmp);
			for (i = 0; i < N_mat; i++)
			{
				val[i] = RhoF[i] + k2[i] / 2.0;
			}
			calcVectValue_t1(step_t_1, m, val, k3, tmp);
			for (i = 0; i < N_mat; i++)
			{
				val[i] = RhoF[i] + k3[i];
			}
			calcVectValue_t1(step_t_1, m, val, k4, tmp);

			for (i = 0; i < N_mat; i++)
			{
				RhoF[i] += (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
			}

			if (rp.debug == 1)
			{
				after(m);
				calcRho(m);
				before(m);

				fn = rp.path + "rho" + "_" + to_string(dump_id) + file_name_suffix(cp, 4);
				save_sparse_complex_mtx(fn, m->Rho, 16, false);
			}

			if (t_1_step_id % dupm_step_t_1 == 0)
			{
				after(m);
				calcRho(m);

				characteristics_deep(m, rp, cp, md, pd, dump_id);
				before(m);
				dump_id++;
			}
		}

		curr_time = period * md.T;

		if (rp.ipp == 1)
		{
			cout << endl << "Dump period: " << period << endl;
		}
	}

	after(m);
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
	int is_dump_adr = cp.int_is_dump_adr;

	MKL_Complex16 * rho_in_d = new MKL_Complex16[md.size * md.size];

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

	if (is_dump_adr == 1)
	{
		double * adr_avg = new double[md.size];

		for (int st_id = 0; st_id < md.size; st_id++)
		{
			adr_avg[st_id] = rho_in_d[st_id * md.size + st_id].real;
		}

		string fn = rp.path + "adr_avg" + file_name_suffix(cp, 4);

		int append = 0;
		if (dump_id > 0)
		{
			append = 1;
		}

		save_double_data(fn, adr_avg, md.size, 16, append);

		delete[] adr_avg;
	}

	for (int spin_id = 0; spin_id < cp.N; spin_id++)
	{
		crsMatrix * sigma_i_z = create_sigma_i_z_matrix(m, cp, spin_id);

		crsMatrix * mult_z = new crsMatrix();

		SparseMKLMult(*(sigma_i_z), *(m->Rho), *(mult_z));

		dcomplex curr_trace_z = trace(*mult_z);

		int current_period = pd.dump_periods[dump_id];
		int curr_sign = 0;
		if (current_period % 2 == 0)
		{
			curr_sign = 1;
		}
		else
		{
			curr_sign = -1;
		}

		double val = double(curr_sign) * curr_trace_z.re * double(pd.sign_sigma_z_start[spin_id]);

		delete sigma_i_z;
		delete mult_z;

		crsMatrix * sigma_i_x = create_sigma_i_x_matrix(m, cp, spin_id);
		crsMatrix * mult_x = new crsMatrix();
		SparseMKLMult(*(sigma_i_x), *(m->Rho), *(mult_x));
		dcomplex curr_trace_x = trace(*mult_x);
		delete sigma_i_x;
		delete mult_x;

		crsMatrix * sigma_i_y = create_sigma_i_y_matrix(m, cp, spin_id);

		crsMatrix * mult_y = new crsMatrix();
		SparseMKLMult(*(sigma_i_y), *(m->Rho), *(mult_y));
		dcomplex curr_trace_y = trace(*mult_y);
		delete sigma_i_y;
		delete mult_y;

		pd.magnetization[spin_id][dump_id] = val;
		pd.sigma_x[spin_id][dump_id] = curr_trace_x.re;
		pd.sigma_y[spin_id][dump_id] = curr_trace_y.re;
		pd.sigma_z[spin_id][dump_id] = curr_trace_z.re;
	}

	delete[] rho_in_d;	
}
void characteristics_rate(Model *m, RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd, int type, int period)
{
	for (int spin_id = 0; spin_id < cp.N; spin_id++)
	{
		if (type == 2)
		{
			if (pd.rate_flag[spin_id] == 0)
			{
				if (pd.rate_ampl_curr[spin_id] < cp.rate_ampl_part * pd.rate_ampl_start[spin_id])
				{
					pd.rate_result[spin_id] = period;
					pd.rate_flag[spin_id] = 1;
				}
			}

			pd.rate_max_curr[spin_id] = -1.0e16;
			pd.rate_min_curr[spin_id] = 1.0e16;
			pd.rate_ampl_curr[spin_id] = 0.0;
		}
		else
		{
			crsMatrix * sigma_i_z = create_sigma_i_z_matrix(m, cp, spin_id);
			crsMatrix * mult_z = new crsMatrix();
			SparseMKLMult(*(sigma_i_z), *(m->Rho), *(mult_z));
			dcomplex curr_trace_z = trace(*mult_z);
			delete sigma_i_z;
			delete mult_z;

			double sigma_z_curr = curr_trace_z.re;

			if (type == 0)
			{
				if (sigma_z_curr > pd.rate_max_curr[spin_id])
				{
					pd.rate_max_curr[spin_id] = sigma_z_curr;
					pd.rate_ampl_start[spin_id] = pd.rate_max_curr[spin_id] - pd.rate_min_curr[spin_id];
					pd.rate_ampl_curr[spin_id] = pd.rate_ampl_start[spin_id];
				}
				if (sigma_z_curr < pd.rate_min_curr[spin_id])
				{
					pd.rate_min_curr[spin_id] = sigma_z_curr;
					pd.rate_ampl_start[spin_id] = pd.rate_max_curr[spin_id] - pd.rate_min_curr[spin_id];
					pd.rate_ampl_curr[spin_id] = pd.rate_ampl_start[spin_id];
				}
			}
			else if (type == 1)
			{
				if (sigma_z_curr > pd.rate_max_curr[spin_id])
				{
					pd.rate_max_curr[spin_id] = sigma_z_curr;
					pd.rate_ampl_curr[spin_id] = pd.rate_max_curr[spin_id] - pd.rate_min_curr[spin_id];
				}
				if (sigma_z_curr < pd.rate_min_curr[spin_id])
				{
					pd.rate_min_curr[spin_id] = sigma_z_curr;
					pd.rate_ampl_curr[spin_id] = pd.rate_max_curr[spin_id] - pd.rate_min_curr[spin_id];
				}

			}
		}
	}
}
void characteristics_deep(Model *m, RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd, int dump_id)
{
	int is_dump_adr = cp.int_is_dump_adr;

	MKL_Complex16 * rho_in_d = new MKL_Complex16[md.size * md.size];

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

	if (is_dump_adr == 1)
	{
		double * adr_avg = new double[md.size];

		for (int st_id = 0; st_id < md.size; st_id++)
		{
			adr_avg[st_id] = rho_in_d[st_id * md.size + st_id].real;
		}

		string fn = rp.path + "adr_avg" + file_name_suffix(cp, 4);

		int append = 0;
		if (dump_id > 0)
		{
			append = 1;
		}

		save_double_data(fn, adr_avg, md.size, 16, append);

		delete[] adr_avg;
	}

	delete[] rho_in_d;

	string fn;

	for (int spin_id = 0; spin_id < cp.N; spin_id++)
	{
		crsMatrix * sigma_i_z = create_sigma_i_z_matrix(m, cp, spin_id);

		if (rp.debug)
		{
			fn = rp.path + "sigma_z_mtx_" + to_string(spin_id) + "_" + to_string(dump_id) + file_name_suffix(cp, 4);
			save_sparse_complex_mtx(fn, sigma_i_z, 16, false);
		}

		crsMatrix * mult_z = new crsMatrix();
		SparseMKLMult(*(sigma_i_z), *(m->Rho), *(mult_z));
		dcomplex curr_trace_z = trace(*mult_z);
		delete sigma_i_z;
		delete mult_z;

		crsMatrix * sigma_i_x = create_sigma_i_x_matrix(m, cp, spin_id);

		if (rp.debug)
		{
			fn = rp.path + "sigma_x_mtx_" + to_string(spin_id) + "_" + to_string(dump_id) + file_name_suffix(cp, 4);
			save_sparse_complex_mtx(fn, sigma_i_x, 16, false);
		}

		crsMatrix * mult_x = new crsMatrix();
		SparseMKLMult(*(sigma_i_x), *(m->Rho), *(mult_x));
		dcomplex curr_trace_x = trace(*mult_x);
		delete sigma_i_x;
		delete mult_x;

		crsMatrix * sigma_i_y = create_sigma_i_y_matrix(m, cp, spin_id);

		if (rp.debug)
		{
			fn = rp.path + "sigma_y_mtx_" + to_string(spin_id) + "_" + to_string(dump_id) + file_name_suffix(cp, 4);
			save_sparse_complex_mtx(fn, sigma_i_y, 16, false);
		}

		crsMatrix * mult_y = new crsMatrix();
		SparseMKLMult(*(sigma_i_y), *(m->Rho), *(mult_y));
		dcomplex curr_trace_y = trace(*mult_y);
		delete sigma_i_y;
		delete mult_y;

		pd.sigma_x[spin_id][dump_id] = curr_trace_x.re;
		pd.sigma_y[spin_id][dump_id] = curr_trace_y.re;
		pd.sigma_z[spin_id][dump_id] = curr_trace_z.re;
	}
}

void init_sign_sigma_z_start(Model* m, ConfigParam &cp, MainData &md, PropData &pd)
{
	for (int spin_id = 0; spin_id < cp.N; spin_id++)
	{
		crsMatrix * sigma_i_z = create_sigma_i_z_matrix(m, cp, spin_id);

		crsMatrix * mult = new crsMatrix();

		SparseMKLMult(*(sigma_i_z), *(m->Rho), *(mult));

		dcomplex curr_trace = trace(*mult);

		if (curr_trace.re > 0)
		{
			pd.sign_sigma_z_start[spin_id] = 1;
		}
		else
		{
			pd.sign_sigma_z_start[spin_id] = -1;
		}

		delete sigma_i_z;
		delete mult;
	}	
}

void init_diss1_A1_rows(int * rows, int N, int i, ConfigParam &cp)
{
	int size = pow(2, N);
	int shift = pow(2, i);
	int chunk = 2 * shift;
	int area = pow(2, i + 2);

	rows[0] = 0;

	for (int st_id = 1; st_id < size + 1; st_id++)
	{
		if ((st_id % area >= shift + 1) && (st_id % area <= shift + chunk))
		{
			rows[st_id] = rows[st_id - 1] + 1;
		}
		else
		{
			rows[st_id] = rows[st_id - 1];
		}
	}
}
void init_diss1_A1_cols(int * cols, int N, int i, ConfigParam &cp)
{
	int size = pow(2, N);
	int shift = pow(2, i);
	int chunk = 2 * shift;
	int area = pow(2, i + 2);

	int num_filled = 0;

	for (int st_id = 0; st_id < size; st_id++)
	{
		if ((st_id % area >= shift) && (st_id % area < shift + chunk))
		{
			cols[num_filled] = st_id;
			num_filled++;
		}
	}
}
void init_diss1_A1_vals(dcomplex * vals, int N, int i, ConfigParam &cp)
{
	int size = pow(2, N);

	int num_filled = 0;
	for (int st_id = 0; st_id < size; st_id++)
	{
		if (bit_at(st_id, i) > bit_at(st_id, i + 1))
		{
			vals[num_filled].re = 1.0;
			vals[num_filled].im = 0.0;

			num_filled++;
		}

		if (bit_at(st_id, i) < bit_at(st_id, i + 1))
		{
			vals[num_filled].re = -1.0;
			vals[num_filled].im = 0.0;

			num_filled++;
		}
	}
}
void init_diss1_A2_rows(int * rows, int N, int i, ConfigParam &cp)
{
	int size = pow(2, N);
	int shift = pow(2, i);
	int chunk = 2 * shift;
	int area = pow(2, i + 2);

	rows[0] = 0;

	for (int st_id = 1; st_id < size + 1; st_id++)
	{
		if ((st_id % area >= shift + 1) && (st_id % area <= shift + chunk))
		{
			rows[st_id] = rows[st_id - 1] + 1;
		}
		else
		{
			rows[st_id] = rows[st_id - 1];
		}
	}
}
void init_diss1_A2_cols(int * cols, int N, int i, ConfigParam &cp)
{
	int size = pow(2, N);
	int shift = pow(2, i);
	int chunk = 2 * shift;
	int area = pow(2, i + 2);

	int num_filled = 0;

	for (int st_id = 0; st_id < size; st_id++)
	{
		int mod = st_id % area;

		if ((mod >= shift) && (mod < shift + chunk))
		{
			if (mod < 2 * shift)
			{
				cols[num_filled] = st_id + shift;
				num_filled++;
			}
			else
			{
				cols[num_filled] = st_id - shift;
				num_filled++;
			}
		}
	}
}
void init_diss1_A2_vals(dcomplex * vals, int N, int i, ConfigParam &cp)
{
	int size = pow(2, N);

	int num_filled = 0;
	for (int st_id = 0; st_id < size; st_id++)
	{
		if (bit_at(st_id, i) > bit_at(st_id, i + 1))
		{
			vals[num_filled].re = sin(-cp.dp);
			vals[num_filled].im = -cos(-cp.dp);

			num_filled++;
		}

		if (bit_at(st_id, i) < bit_at(st_id, i + 1))
		{
			vals[num_filled].re = -sin(cp.dp);
			vals[num_filled].im = cos(cp.dp);

			num_filled++;
		}
	}
}

void f_basis_init(Model* model, RunParam &rp, ConfigParam &cp, MainData &md)
{
	string fn;

	double time = omp_get_wtime();
	double init_time = time;

	time = omp_get_wtime() - init_time;
	cout << "Time of createModel: " << time << endl << endl;

	init_h_J_vector(model, rp, cp, md);
	time = omp_get_wtime() - init_time;
	cout << "Time of init_h_J_vector: " << time << endl << endl;
	if (rp.debug == 1)
	{
		fn = rp.path + "h_j_vec" + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, model->h_J, 16, false);
	}

	init_h_h_z_vector(model, rp, cp, md);
	time = omp_get_wtime() - init_time;
	cout << "Time of init_h_h_z_vector: " << time << endl << endl;
	if (rp.debug == 1)
	{
		fn = rp.path + "h_h_z_vec" + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, model->h_h_z, 16, false);
	}

	init_h_h_x_vector(model, rp, cp, md);
	time = omp_get_wtime() - init_time;
	cout << "Time of init_h_h_x_vector: " << time << endl << endl;
	if (rp.debug == 1)
	{
		fn = rp.path + "h_h_x_vec" + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, model->h_h_x, 16, false);
	}

	init_h_s_x_vector(model, rp, cp, md);
	time = omp_get_wtime() - init_time;
	cout << "Time of init_h_s_x_vector: " << time << endl << endl;
	if (rp.debug == 1)
	{
		fn = rp.path + "h_s_x_vec" + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, model->h_s_x, 16, false);
	}

	init_H0(model);

	if (rp.issmtx == 1)
	{
		fn = rp.path + "H0" + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, model->H0, 16, false);

		fn = rp.path + "H1" + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, model->H_s_x, 16, false);
	}

	if (cp.dt == 1)
	{
		init_diss_1(model, rp, cp, md);
	}
	else if (cp.dt == 2)
	{
		init_diss_2(model, rp, cp, md);
	}
	else if (cp.dt == 3)
	{
		init_diss_3(model, rp, cp, md);
	}
	else if (cp.dt == 4)
	{
		init_diss_4(model, rp, cp, md);
	}
	else
	{
		init_diss_1(model, rp, cp, md);
	}

	time = omp_get_wtime() - init_time;
	cout << "time of init_diss_" << to_string(cp.dt) << ": " << time << endl << endl;

	calc_Q_J_s(model);
	time = omp_get_wtime() - init_time;
	cout << "time of calc_Q_J_s: " << time << endl << endl;
	if (rp.debug == 1)
	{
		fn = rp.path + "Q_J_s" + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, model->Q_J_s, 16, false);
	}

	calc_Q_h_z_s(model);
	time = omp_get_wtime() - init_time;
	cout << "time of calc_Q_h_z_s: " << time << endl << endl;
	if (rp.debug == 1)
	{
		fn = rp.path + "Q_h_z_s" + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, model->Q_h_z_s, 16, false);
	}

	calc_Q_h_x_s(model);
	time = omp_get_wtime() - init_time;
	cout << "time of calc_Q_h_x_s: " << time << endl << endl;
	if (rp.debug == 1)
	{
		fn = rp.path + "Q_h_x_s" + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, model->Q_h_x_s, 16, false);
	}

	calc_Q_s_x_s(model);
	time = omp_get_wtime() - init_time;
	cout << "time of calc_Q_s_x_s: " << time << endl << endl;
	if (rp.debug == 1)
	{
		fn = rp.path + "Q_s_x_s" + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, model->Q_s_x_s, 16, false);
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

	calc_G_0_s(model);
	time = omp_get_wtime() - init_time;
	cout << "time of calc_G_0_s: " << time << endl << endl;

	if (rp.debug == 1)
	{
		fn = rp.path + "G_0_s" + file_name_suffix(cp, 4);
		save_sparse_complex_mtx(fn, model->G_0_s, 16, false);
	}

	calc_G_1_s(model);
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

