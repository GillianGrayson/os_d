#include "CalcRs.h"
#include "sortTensor.h"
#include <omp.h>
#include <math.h>
#include "coef_coord.h"
#include "stdToCrs.h"

#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#endif

//unsigned int countSub_ijk(Tensor_Coordinates *m_ijk, crsMatrix * Ls, unsigned int from, unsigned int to)
//{
//  unsigned int cnt = 0;
//  for (unsigned int i = from; i < to; i++)
//  {
//    unsigned int j = m_ijk->coord1[i];
//    if ((Ls->RowIndex[j + 1] - Ls->RowIndex[j]) > 0)
//    {
//      cnt++;
//    }
//  }
//  return cnt;
//}

//unsigned int selectSub_ijk(Tensor_Coordinates *m_ijk, crsMatrix * Ls, 
//  unsigned int from, unsigned int to, Tensor_Coordinates *res_ijk)
//{
//  unsigned long long int step1 = Ls->N;
//  unsigned long long int step2 = step1 * step1;
//  unsigned int cnt = 0;
//  for (int i = from; i < to; i++)
//  {
//    unsigned int j = m_ijk->coord1[i];
//    if ((Ls->RowIndex[j + 1] - Ls->RowIndex[j]) > 0)
//    {
//      res_ijk->coord1[cnt] = m_ijk->coord1[i];
//      res_ijk->coord2[cnt] = m_ijk->coord2[i];
//      res_ijk->coord3[cnt] = m_ijk->coord3[i];
//      res_ijk->data[cnt] = m_ijk->data[i];
//      res_ijk->hash[cnt] = m_ijk->coord1[i] + m_ijk->coord2[i] * step2 + step1 * m_ijk->coord3[i];
//      cnt++;
//    }
//  }
//  sort_matrix(res_ijk);
//  return cnt;
//}

void multM_ijk(Tensor_Coordinates *m_ijk, crsMatrix * Ls,
  ulli from, ulli to, bool conj, Tensor_Coordinates *res_ijk)
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
  m_ijk->k = cnt+1;
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

ulli multTmpRsSTD(ulli N_mat, unsigned int * indexF, unsigned int * indexZ,
	Tensor_Coordinates *z_ijk, Tensor_Coordinates *f_ijk, vector<map<int, dcomplex> > & mat)
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

ulli multTmpRs(ulli N_mat, unsigned int * indexF, unsigned int * indexZ,
  Tensor_Coordinates *z_ijk, Tensor_Coordinates *f_ijk, Tensor_Coordinates *res_ijk)
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
      Z_ijk->data[cnt]   = f_ijk->data[j1];
      Z_ijk->hash[cnt]   = f_ijk->hash[j1];

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
        Z_ijk->data[cnt].im =  d_ijk->data[j2].re;
        Z_ijk->hash[cnt]   = d_ijk->hash[j2];
        j2++;
      }
      else
      {
        Z_ijk->coord1[cnt] = d_ijk->coord1[j2];
        Z_ijk->coord2[cnt] = d_ijk->coord2[j2];
        Z_ijk->coord3[cnt] = d_ijk->coord3[j2];
        Z_ijk->data[cnt]   = f_ijk->data[j1];
        Z_ijk->data[cnt].re += -d_ijk->data[j2].im;
        Z_ijk->data[cnt].im += d_ijk->data[j2].re;
        Z_ijk->hash[cnt] = d_ijk->hash[j2];
        j1++;
        j2++;
      }
    }
    cnt++;
  }
  while ((j1 < f_ijk->k) )
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
    Z_ijk->data[cnt].im =  d_ijk->data[j2].re;
    Z_ijk->hash[cnt] = d_ijk->hash[j2];
    j2++;
    cnt++;
  }


  return cnt;
}


//crsMatrix* calcSubRs(Model *m, crsMatrix ** f_d_mat, crsMatrix ** f_d_H_mat, ulli start, ulli finish);


void calcRs(Model *m)
{
  dcomplex sum;
  sum.re = 1.0;
  sum.im = 0.0;

  int N_mat = m->N_mat;

  crsMatrix * Rs_tmp;
  crsMatrix ** RsTh;
  crsMatrix * Rs;// = new crsMatrix(m->N_mat, 1);

  
  Tensor_Coordinates * subF_ijk;
  Tensor_Coordinates * subD_ijk;
  Tensor_Coordinates * subZ_ijk;
  Tensor_Coordinates * l_F_ijk;
  Tensor_Coordinates * l_Z_ijk;
//  Tensor_Coordinates * tmpRs1;
//  Tensor_Coordinates * tmpRs2;

  ulli cnt;

  unsigned int *indexZ = new unsigned int[N_mat + 1];
  unsigned int *indexF = new unsigned int[N_mat + 1];

  //unsigned int cntSubFijk = countSub_ijk(m->f_ijk, m->l_mat, 0, m->f_ijk->k);
  //unsigned int cntSubDijk = countSub_ijk(m->d_ijk, m->l_mat, 0, m->d_ijk->k);
  //printf("cntSubFijk %u\n", cntSubFijk);
  //printf("cntSubDijk %u\n", cntSubDijk);
  //subF_ijk = create_matrix(cntSubFijk);
  //subD_ijk = create_matrix(cntSubDijk);
  //selectSub_ijk(m->f_ijk, m->l_mat, 0, m->f_ijk->k, subF_ijk);
  //selectSub_ijk(m->d_ijk, m->l_mat, 0, m->d_ijk->k, subD_ijk);
  
  cnt = dijk_coord_sym(m->l_mat, m->N + 1);
  printf("cntSubDijk %u\n", cnt);
  subD_ijk = new Tensor_Coordinates;
  dijk_coord_ch(subD_ijk, m->l_mat, cnt+1, m->N + 1);
  
  fprintf(m->memlog, "dijk_coord_ch(subD_ijk, m->l_mat, cnt+1, m->N + 1); %lf\n", number_of_allocs);
  fflush(m->memlog);


  subF_ijk = m->f_ijk;
  sort_matrix(subF_ijk);
  sort_matrix(subD_ijk);
  
  calcZ_ijk(subF_ijk, subD_ijk, subZ_ijk);

  l_F_ijk = create_matrix(subF_ijk->k);
  l_Z_ijk = create_matrix(subZ_ijk->k);

  multM_ijk(subZ_ijk, m->l_mat, 0, subZ_ijk->k, false, l_Z_ijk);
  multM_ijk(subF_ijk, m->l_mat, 0, subF_ijk->k, true, l_F_ijk);

  fprintf(m->memlog, "multM_ijk %lf\n", number_of_allocs);
  fflush(m->memlog);

  reduseM_ijk(l_Z_ijk, N_mat, 0, l_Z_ijk->k);
  reduseM_ijk(l_F_ijk, N_mat, 0, l_F_ijk->k);
  
  ismap = true;
  vector<map<int, dcomplex> > mat(N_mat);
  
  createIndex(l_Z_ijk, N_mat, indexZ);
  createIndex(l_F_ijk, N_mat, indexF);
  //ulli cntRsTmp1 = cntTmpRs(N_mat, indexF, indexZ);
  //tmpRs1 = create_matrix(cntRsTmp1);
  //multTmpRs(N_mat, indexF, indexZ, l_Z_ijk, l_F_ijk, tmpRs1);
  //reduseM_ijk(tmpRs1, N_mat, 0, tmpRs1->k);

  multTmpRsSTD(N_mat, indexF, indexZ, l_Z_ijk, l_F_ijk, mat);

  fprintf(m->memlog, "multTmpRsSTD %lf\n", number_of_allocs);
  fflush(m->memlog);


  for (ulli i = 0; i < subZ_ijk->k; i++)
  {
    subZ_ijk->data[i].im = -subZ_ijk->data[i].im;
  }

  l_Z_ijk->k = subZ_ijk->k;
  l_F_ijk->k = subF_ijk->k;
  multM_ijk(subZ_ijk, m->l_mat, 0, subZ_ijk->k, true, l_Z_ijk);
  multM_ijk(subF_ijk, m->l_mat, 0, subF_ijk->k, false, l_F_ijk);
  reduseM_ijk(l_Z_ijk, N_mat, 0, l_Z_ijk->k);
  reduseM_ijk(l_F_ijk, N_mat, 0, l_F_ijk->k);
  createIndex(l_Z_ijk, N_mat, indexZ);
  createIndex(l_F_ijk, N_mat, indexF);
  //ulli cntRsTmp2 = cntTmpRs(N_mat, indexF, indexZ);
  //tmpRs2 = create_matrix(cntRsTmp2);
  //multTmpRs(N_mat, indexF, indexZ, l_Z_ijk, l_F_ijk, tmpRs2);
  //reduseM_ijk(tmpRs2, N_mat, 0, tmpRs2->k);

  multTmpRsSTD(N_mat, indexF, indexZ, l_Z_ijk, l_F_ijk, mat);

  fprintf(m->memlog, "multTmpRsSTD %lf\n", number_of_allocs);
  fflush(m->memlog);

  ismap_print = true;

  crsMatrix * Rs_std = stdToCrs(mat, N_mat);

//  free_matrix(subF_ijk);
  free_matrix(subD_ijk);
  free_matrix(subZ_ijk);
  free_matrix(l_F_ijk );
  free_matrix(l_Z_ijk );
//  delete subF_ijk;
  delete subD_ijk;
  delete subZ_ijk;
  delete l_F_ijk;
  delete l_Z_ijk;


 // crsMatrix *RsCo = new crsMatrix(N_mat, tmpRs2->k);
 //
 // for (ulli i = 0; i < tmpRs1->k; i++)
 // {
 //   RsCo->Value[i].re = tmpRs1->data[i].re + tmpRs2->data[i].re;
 //   RsCo->Value[i].im = tmpRs1->data[i].im + tmpRs2->data[i].im;
 //   RsCo->Col[i] = tmpRs1->coord2[i];
 // }
 //
 // for (int i = 0; i < N_mat; i++)
 // {
 //   RsCo->RowIndex[i] = N_mat+1;
 // }
 // RsCo->RowIndex[N_mat] = tmpRs1->k;
 //
 // for (ulli i = 0; i < tmpRs1->k; i++)
 // {
 //   unsigned int j = tmpRs1->coord3[i];
 //   if (RsCo->RowIndex[j] == N_mat + 1)
 //   {
 //     RsCo->RowIndex[tmpRs1->coord3[i]] = i;
 //   }
 // }
 //
 // for (int i = N_mat - 1; i >= 0; i--)
 // {
 //   if (RsCo->RowIndex[i] == N_mat + 1)
 //   {
 //     RsCo->RowIndex[i] = RsCo->RowIndex[i + 1];
 //   }
 // }

 // printMatrixVal(RsCo);
//  printMatrixVal(Rs_std);
  Rs = Rs_std;



//  free_matrix(tmpRs1);
//  free_matrix(tmpRs2);
//  delete tmpRs1;
//  delete tmpRs2;
  delete[] indexZ;
  delete[] indexF;

/////  Rs = new crsMatrix(m->N_mat, 1);
/////  Rs->Col[0] = 1;
/////  Rs->RowIndex[0] = 1;
/////  for(int i = 1; i <= N_mat; i++)
/////  {
/////    Rs->RowIndex[i] = 1;
/////  }
/////
/////  int nThread = 1;
///////  omp_set_num_threads(3);
/////#pragma omp parallel
/////  {
/////    #pragma omp single
/////    nThread = omp_get_num_threads();
/////  }
/////  printf("nThread %d \n", nThread);
/////
/////  int SubTask = 1;
/////  int CntTask = nThread * SubTask;
/////
/////  RsTh = new crsMatrix *[CntTask];
/////
/////  crsMatrix ** f_mat = m->f_mat;
/////  crsMatrix ** f_H_mat = m->f_H_mat;
/////  crsMatrix ** d_mat = m->d_mat;
/////
/////  for(int i = 0; i < N_mat; i++)
/////  {
/////    toOneBase(*(f_mat[i]));
/////    toOneBase(*(d_mat[i]));
/////    toOneBase(*(f_H_mat[i]));
/////  }
/////
/////  dcomplex sum_i;
/////  sum_i.re = 0.0;
/////  sum_i.im = 1.0;
/////
/////  crsMatrix ** f_d_mat = new crsMatrix *[N_mat];
/////  crsMatrix ** f_d_H_mat = new crsMatrix *[N_mat];
/////
/////#pragma omp for schedule(dynamic)
/////  for(int i = 0; i < N_mat; i++)
/////  {
/////    f_d_mat[i] = new crsMatrix;
/////    SparseMKLAddOne(*(f_mat[i]), sum_i, *(d_mat[i]), *(f_d_mat[i]), true);
/////    f_d_H_mat[i] = new crsMatrix(*(f_d_mat[i]));
/////    for(int conj_i = 0; conj_i < f_d_H_mat[i]->NZ; conj_i++)
/////        f_d_H_mat[i]->Value[conj_i].im = -f_d_H_mat[i]->Value[conj_i].im;
/////  }
/////  
/////  int size = N_mat / CntTask;
///////  int *start = new int[CntTask];
///////  int *finish = new int[CntTask];
/////
///////  start[0] = 0;
///////  finish[0] = size;
///////  if((N_mat % CntTask) != 0) finish[0]++;
/////  printf("N_mat %d \n", N_mat);
///////  printf("%d %d \n", start[0], finish[0]);
///////  for(int i = 1; i < CntTask; i++)
///////  {
///////    start[i] = finish[i - 1];
///////    finish[i] = start[i] + size;
///////    if((N_mat % CntTask) >  i)finish[i]++;
///////    printf("%d %d \n", start[i], finish[i]);
///////  }
/////
/////
/////  int id = 0;
/////  int portion = 10;
/////
/////  int maxNZ = 0;
/////  int workThread = MIN(1, nThread);
/////
/////  #pragma omp parallel
/////  {
/////    int tid = omp_get_thread_num();
/////    int localTskID;
/////    if (tid < workThread)
/////    {
/////      while (true)
/////      {
/////#pragma omp critical
/////      {
/////        localTskID = id;
/////        id++;
/////      }
/////        int start = localTskID * portion;
/////        int finish = (localTskID + 1) * portion;
/////        if (start > N_mat) break;
/////        if (finish > N_mat) {
/////          finish = N_mat;
/////          printf("%d %d \n", start, finish);
/////        }
/////        if (finish - start > 0)
/////        {
/////          RsTh[tid] = calcSubRs(m, f_d_mat, f_d_H_mat, start, finish);
/////        }
/////        else
/////        {
/////          RsTh[tid] = NULL;
/////        }
/////#pragma omp critical
/////        {
/////          if (RsTh[tid] != NULL)
/////          {
/////            Rs_tmp = new crsMatrix;
/////            SparseMKLAddOne(*Rs, sum, *RsTh[tid], *Rs_tmp);
/////            if (maxNZ < RsTh[tid]->NZ) maxNZ = RsTh[tid]->NZ;
/////            delete RsTh[tid];
/////            delete Rs;
/////            Rs = Rs_tmp;
/////          }
/////        }
/////      }
/////    }
/////  }
/////  printf("*************Rs maxNx %d\n", maxNZ);
/////  //delete [] start;
/////  //delete [] finish;
/////  
///////  for(int i = 0; i < CntTask; i++)
///////  {
///////      Rs_tmp = new crsMatrix;
///////      SparseMKLAddOne(*Rs, sum, *RsTh[i], *Rs_tmp);
///////      delete RsTh[i];
///////      delete Rs;
///////      Rs = Rs_tmp;
///////  }
/////  delete[] RsTh;
/////
/////  
/////  for(int i = 0; i < N_mat; i++)
/////  {
/////    delete f_d_mat[i];
/////    delete f_d_H_mat[i];
/////  }
/////  delete[] f_d_mat;
/////  delete[] f_d_H_mat;
/////  //Rs = calcSubRs(m, 0, N_mat);
/////  toZeroBase(*Rs);
/////  for(int i = 0; i < N_mat; i++)
/////  {
/////    toZeroBase(*(f_mat[i]));
/////    toZeroBase(*(d_mat[i]));
/////    toZeroBase(*(f_H_mat[i]));
/////  }
///////  delete FkT;
///////  delete FiT;
//  printMatrixVal(Rs);

  double mm = -m->conf.g / 4.0;
  for (unsigned int i = 0; i < Rs->NZ; i++)
  {
    Rs->Value[i].re *= mm;
    Rs->Value[i].im *= mm;
  }
  m->Rs = Rs;
}


//////crsMatrix* calcSubRs(Model *m, 
//////                     crsMatrix ** f_d_mat, crsMatrix ** f_d_H_mat, 
//////                     int start, int finish)
//////{
//////  int N_mat = m->N_mat;
//////  dcomplex sum_i;
//////  sum_i.re = 0.0;
//////  sum_i.im = 1.0;
//////
//////  dcomplex sum;
//////  sum.re = 1.0;
//////  sum.im = 0.0;
//////
//////  crsMatrix * tmp;
//////  crsMatrix * MatSDi, * MatSDk;
//////  crsMatrix * M1, *M2, *Rsum;
//////
//////  crsMatrix * Rs_tmp;
//////  crsMatrix * Rs = new crsMatrix(m->N_mat, 1);
//////  
//////  Rs->Col[0] = 1;
//////  Rs->RowIndex[0] = 1;
//////  for(int i = 1; i <= N_mat; i++)
//////  {
//////    Rs->RowIndex[i] = 1;
//////  }
//////
////////  crsMatrix ** f_mat = m->f_mat;
////////  crsMatrix ** f_H_mat = m->f_H_mat;
////////  crsMatrix ** d_mat = m->d_mat;
//////  crsMatrix *  As    = m->a_mat;
//////
//////  M1 = new crsMatrix;
//////  M2 = new crsMatrix;
////////  MatSDi = new crsMatrix;
////////  MatSDk = new crsMatrix;
//////  Rsum = new crsMatrix;
//////      
////////  for(int i = 0; i < N_mat; i++)
////////  {
////////    toOneBase(*(f_mat[i]));
////////    toOneBase(*(d_mat[i]));
////////    toOneBase(*(f_H_mat[i]));
////////  }
//////      
//////  for(int i = start; i < finish; i++)
//////  {
////////    SparseMKLAddOne(*(f_mat[i]), sum_i, *(d_mat[i]), *MatSDi, true);
//////    MatSDi = f_d_mat[i];
//////
//////    int st = As->RowIndex[i];
//////    int fn = As->RowIndex[i + 1];
//////
//////    for(int j = st; j < fn; j++)
//////    {
//////      int k = As->Col[j];
//////      //SparseMKLAddOne(*(f_mat[k]), sum_i, *(d_mat[k]), *MatSDk, true);
//////      MatSDk = f_d_H_mat[k];
//////
//////      //for(int conj_i = 0; conj_i < MatSDk->NZ; conj_i++)
//////      //  MatSDk->Value[conj_i].im = -MatSDk->Value[conj_i].im;
//////
//////      SparseMKLMultOne(*(f_H_mat[k]), *MatSDi, *M1, true);
//////      SparseMKLMultOne(*(f_H_mat[i]), *MatSDk, *M2, true);
//////
//////      SparseMKLAddOne(*M1, sum, *M2, *Rsum, true);
//////
//////      dcomplex val1 = As->Value[j];
//////      for(int ii = 0; ii < Rsum->NZ; ii ++)
//////      {
//////        dcomplex val2 = Rsum->Value[ii];
//////        Rsum->Value[ii].re = val1.re *  val2.re - val1.im *  val2.im;
//////        Rsum->Value[ii].im = val1.re *  val2.im + val1.im *  val2.re;
//////      }
//////
//////      Rs_tmp = new crsMatrix;
//////      SparseMKLAddOne(*Rs, sum, *Rsum, *Rs_tmp);
//////      delete Rs;
//////      Rs = Rs_tmp;
//////    }
//////  }
//////
////////  toZeroBase(*Rs);
////////  for(int i = 0; i < N_mat; i++)
////////  {
////////    toZeroBase(*(f_mat[i]));
////////    toZeroBase(*(d_mat[i]));
////////    toZeroBase(*(f_H_mat[i]));
////////  }
////////  delete FkT;
////////  delete FiT;
//////
//////  delete Rsum;
//////  
////////  delete MatSDi;
////////  delete MatSDk;
//////  
//////  delete M1;
//////  delete M2;
//////
//////  //double mm = -m->conf.g / 4.0 / m->N;
//////  //for(int i = 0; i < Rs->NZ; i++)
//////  //{
//////  //  Rs->Value[i].re *= mm;
//////  //  Rs->Value[i].im *= mm;
//////  //}
//////  return Rs;
//////}