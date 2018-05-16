#include "CalcRs.h"
#include "sortTensor.h"
#include <omp.h>
#include <math.h>
#include "coef_coord.h"
#include "stdToCrs.h"

#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#endif

unsigned int countSub_ijk(Tensor_Coordinates *m_ijk, crsMatrix * Ls, unsigned int from, unsigned int to)
{
	unsigned int cnt = 0;
	for (unsigned int i = from; i < to; i++)
	{
		unsigned int j = m_ijk->coord1[i];
		if ((Ls->RowIndex[j + 1] - Ls->RowIndex[j]) > 0)
		{
			cnt++;
		}
	}
	return cnt;
}

unsigned int selectSub_ijk(Tensor_Coordinates *m_ijk, crsMatrix * Ls, unsigned int from, unsigned int to, Tensor_Coordinates *res_ijk)
{
	unsigned long long int step1 = Ls->N;
	unsigned long long int step2 = step1 * step1;
	unsigned int cnt = 0;
	for (int i = from; i < to; i++)
	{
		unsigned int j = m_ijk->coord1[i];
		if ((Ls->RowIndex[j + 1] - Ls->RowIndex[j]) > 0)
		{
			res_ijk->coord1[cnt] = m_ijk->coord1[i];
			res_ijk->coord2[cnt] = m_ijk->coord2[i];
			res_ijk->coord3[cnt] = m_ijk->coord3[i];
			res_ijk->data[cnt] = m_ijk->data[i];
			res_ijk->hash[cnt] = m_ijk->coord1[i] + m_ijk->coord2[i] * step2 + step1 * m_ijk->coord3[i];
			cnt++;
		}
	}
	sort_matrix(res_ijk);
	return cnt;
}

void multM_ijk(Tensor_Coordinates *m_ijk, crsMatrix * Ls, ulli from, ulli to, bool conj, Tensor_Coordinates *res_ijk)
{
	for (ulli i = from; i < to; i++)
	{
		unsigned int j = m_ijk->coord1[i];
		dcomplex v1 = Ls->Value[Ls->RowIndex[j]];
		if (conj) v1.im = -v1.im;
		dcomplex v2 = m_ijk->data[i];
		res_ijk->coord1[i] = m_ijk->coord1[i];
		res_ijk->coord2[i] = m_ijk->coord2[i];
		res_ijk->coord3[i] = m_ijk->coord3[i];
		res_ijk->data[i].re = v1.re * v2.re - v1.im * v2.im;
		res_ijk->data[i].im = v1.re * v2.im + v1.im * v2.re;
		res_ijk->hash[i] = m_ijk->hash[i];
	}
}

void reduseM_ijk(Tensor_Coordinates *m_ijk, ulli N_Mat, ulli from, ulli to)
{
	ulli step1 = N_Mat;
	ulli cnt = from;
	for (ulli i = from + 1; i < to; i++)
	{
		m_ijk->coord1[i] = 0;
		if ((m_ijk->coord2[cnt] != m_ijk->coord2[i]) ||
			(m_ijk->coord3[cnt] != m_ijk->coord3[i]))
		{
			cnt++;
			m_ijk->coord1[i] = 0;
			m_ijk->coord2[cnt] = m_ijk->coord2[i];
			m_ijk->coord3[cnt] = m_ijk->coord3[i];
			m_ijk->data[cnt] = m_ijk->data[i];
			m_ijk->hash[cnt] = m_ijk->hash[i];
		}
		else
		{
			m_ijk->data[cnt].re += m_ijk->data[i].re;
			m_ijk->data[cnt].im += m_ijk->data[i].im;
		}
	}
	m_ijk->k = cnt + 1;
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

ulli cntTmpRs(unsigned int N_mat, unsigned int * indexF, unsigned int * indexZ)
{
	ulli cnt = 0;
	for (int i = 0; i < N_mat; i++)
	{
		cnt += (indexF[i + 1] - indexF[i])*(indexZ[i + 1] - indexZ[i]);
	}
	return cnt;
}

ulli multTmpRsSTD(ulli N_mat, unsigned int * indexF, unsigned int * indexZ, Tensor_Coordinates *z_ijk, Tensor_Coordinates *f_ijk, vector<map<int, dcomplex> > & mat)
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
			//res_ijk->coord1[cnt + j] = 0;
			//res_ijk->coord2[cnt + j] = z_ijk->coord3[j2];
			//res_ijk->coord3[cnt + j] = f_ijk->coord3[j1];
			//res_ijk->hash[cnt + j] = res_ijk->coord2[cnt + j] + res_ijk->coord3[cnt + j] * step1;
			dcomplex v1 = z_ijk->data[j2];
			dcomplex v2 = f_ijk->data[j1];


			v = mat[f_ijk->coord3[j1]][z_ijk->coord3[j2]];
			v.re += v1.re * v2.re - v1.im * v2.im;
			v.im += v1.re * v2.im + v1.im * v2.re;
			mat[f_ijk->coord3[j1]][z_ijk->coord3[j2]] = v;
		}
		cnt += c1 * c2;

	}
	//	sort_matrix(res_ijk);
	return cnt;
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
            int fi = l_mat->RowIndex[ind_j+1];
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

ulli multTmpRs(ulli N_mat, unsigned int * indexF, unsigned int * indexZ, Tensor_Coordinates *z_ijk, Tensor_Coordinates *f_ijk, Tensor_Coordinates *res_ijk)
{
	ulli step1 = N_mat;

	unsigned int c1, c2, j1, j2;
	ulli cnt = 0;
	for (unsigned int i = 0; i < N_mat; i++)
	{
		c1 = (indexF[i + 1] - indexF[i]);
		c2 = (indexZ[i + 1] - indexZ[i]);
		for (unsigned int j = 0; j < c1 * c2; j++)
		{
			j1 = indexF[i] + j % c1;
			j2 = indexZ[i] + j / c1;
			res_ijk->coord1[cnt + j] = 0;
			res_ijk->coord2[cnt + j] = z_ijk->coord3[j2];
			res_ijk->coord3[cnt + j] = f_ijk->coord3[j1];
			res_ijk->hash[cnt + j] = res_ijk->coord2[cnt + j] + res_ijk->coord3[cnt + j] * step1;
			dcomplex v1 = z_ijk->data[j2];
			dcomplex v2 = f_ijk->data[j1];
			res_ijk->data[cnt + j].re = v1.re * v2.re - v1.im * v2.im;
			res_ijk->data[cnt + j].im = v1.re * v2.im + v1.im * v2.re;
		}
		cnt += c1 * c2;

	}
	sort_matrix(res_ijk);
	return cnt;
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

	ulli uN = m-> l_mat->N;
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
