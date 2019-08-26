#include "CalcQs.h"
#include "sortTensor.h"
#include <string.h>
#include "coef_coord.h"

//ulli countSelect(Tensor_Coordinates * f_ijk, crsMatrix *hMat, ulli from, ulli to)
//{
//  unsigned int cnt = 0;
//  for (unsigned int i = from; i < to; i++)
//  {
//    int j = f_ijk->coord1[i];
//    cnt += hMat->RowIndex[j + 1] - hMat->RowIndex[j];
//  }
//  return cnt;
//}
//
//void dataSelect(int N_mat, Tensor_Coordinates * f_ijk, crsMatrix *hMat, 
//                unsigned int from, unsigned int to, Tensor_Coordinates * res)
//{
//  unsigned int k = 0;
//  unsigned int cnt = 0;
//  for (unsigned int i = from; i < to; i++)
//  {
//    int j = f_ijk->coord1[i];
//    cnt = hMat->RowIndex[j + 1] - hMat->RowIndex[j];
//    if (cnt > 0)
//    {
//      dcomplex val1 = f_ijk->data[i];
//      val1.im = val1.im;
//      dcomplex val2 = hMat->Value[hMat->RowIndex[j]];
//      res->coord1[k] = f_ijk->coord1[i];
//      res->coord2[k] = f_ijk->coord2[i];
//      res->coord3[k] = f_ijk->coord3[i];
//      res->data[k].re = val1.re * val2.re - val1.im * val2.im;
//      res->data[k].im = val1.re * val2.im + val1.re * val2.im;
//
//      printf("%lf %lf %d %d %d\n", res->data[k].re, res->data[k].im, res->coord1[k], res->coord2[k], res->coord3[k]);
//
//      res->hash[k] = 0;
//      unsigned long long int step = N_mat;
//      res->hash[k] += res->coord2[k] + step * res->coord3[k];
//      k++;
//    }
//  }
//}

void calc_CooQs(int N_mat, Model *m, Tensor_Coordinates * f_ijk, crsMatrix *hMat, crsMatrix *&res)
{
  Tensor_Coordinates * select = new Tensor_Coordinates;
  int *hash = new int[N_mat];
  int *rowi = new int[N_mat+1];
  dcomplex *hash_calc = new dcomplex[N_mat];
  ulli cnt;
  int *hash_col = new int[N_mat];
  int col_cnt;

//  cnt = countSelect(f_ijk, hMat, 0, f_ijk->k);
//
//  select->k = cnt;
//  select->coord1 = new unsigned int[cnt];
//  select->coord2 = new unsigned int[cnt];
//  select->coord3 = new unsigned int[cnt];
//  select->data = new dcomplex[cnt];
//  select->hash = new unsigned long long int[cnt];
//  printf("\n count %u \n\n", cnt);
//
//  dataSelect(N_mat, f_ijk, hMat, 0, f_ijk->k, select);
  cnt = fijk_coord_sym(hMat, m->N + 1);
  fijk_coord_ch(select, hMat, cnt+1, m->N + 1);

  fprintf(m->memlog, "fijk_coord_ch(select, hMat, cnt+1, m->N + 1); %lf\n", number_of_allocs);
  fflush(m->memlog);

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
  for (int i = N_mat-1; i >0; i--)
  {
    if (rowi[i] == -1)
    {
      rowi[i] = rowi[i + 1];
    }
  }
  rowi[N_mat] = select->k;

  res = new crsMatrix(N_mat, 1);
  
  memset(hash, 0, sizeof(int) * N_mat);
  int countNotZero = 0;
  for (int j = 0; j < N_mat; j++)
  {
    col_cnt = 0;
    res->RowIndex[j] = countNotZero;
    int start, finish;
    start = rowi[j];
    finish = rowi[j + 1];
    for (int k = start; k < finish; k++)
    {
      if (hash[select->coord2[k]] == 0)
      {
        hash[select->coord2[k]] = 1;
        hash_col[col_cnt] = select->coord2[k];
        col_cnt++;
        countNotZero++;
      }
    }
    for (int k = 0; k < col_cnt; k++)
    {
      hash[hash_col[k]] = 0;
    }
  }
  res->RowIndex[N_mat] = countNotZero;
  res->setNZ(countNotZero);
  
  memset(hash, 0, sizeof(int) * N_mat);
  countNotZero = 0;
  for (int j = 0; j < N_mat; j++)
  {
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
      hash[res->Col[k]] = 0;
    }
  }

  crsMatrix * R = new crsMatrix(*res);
  Transpose(*res, *R, false);
  Transpose(*R, *res, false);

//  printMatrixVal(res);

  delete R;
  delete[] hash;
  delete[] hash_calc;
  delete[] hash_col;
  delete[] rowi;

  free_matrix(select);
  delete select;
}


void calcQs(Model * m)
{
  int N_mat = m->N_mat;
  crsMatrix * h = m->h;
  
//  crsMatrix ** f_H_mat = m->f_H_mat;
//  crsMatrix * Qs;
//  crsMatrix * resSum, *tmp;
//
//  
//  int rSize = f_H_mat[0]->N;
//  int *hash = new int [rSize];
//  dcomplex *hash_calc = new dcomplex [rSize];
  
  crsMatrix * QsCo;
  calc_CooQs(N_mat, m, m->f_ijk, h, QsCo);

//  unsigned int cnt = fijk_coord_sym(h, m->N + 1);
//  printf("fijk_coord_sym(h, m->N + 1); %u\n", cnt);

//  Qs = new crsMatrix(rSize, 1);
//  
//  int countNotZero = 0;
//  for(int j = 0;j < rSize; j++)
//  {
//    Qs->RowIndex[j] = countNotZero;
//    memset(hash, 0, sizeof(int) * rSize);
//    for(int i = 0; i < h->NZ; i++)
//    {
//      int ind = h->Col[i];
//      int start, finish;
//      start  = f_H_mat[ind]->RowIndex[j];
//      finish = f_H_mat[ind]->RowIndex[j+1];
//      for(int k = start; k < finish; k++)
//      {
//        if(hash[f_H_mat[ind]->Col[k]] == 0)
//        {
//          hash[f_H_mat[ind]->Col[k]] = 1;
//          countNotZero++;
//        }
//      }
//    }
//  }
//  Qs->RowIndex[rSize] = countNotZero;
//  Qs->setNZ(countNotZero);
//
//  countNotZero = 0;
//  for(int j = 0;j < rSize; j++)
//  {
//    int start, finish;
//    memset(hash_calc, 0, sizeof(dcomplex) * rSize);
//    memset(hash, 0, sizeof(int) * rSize);
//    for(int i = 0; i < h->NZ; i++)
//    {
//      int ind = h->Col[i];
//      start  = f_H_mat[ind]->RowIndex[j];
//      finish = f_H_mat[ind]->RowIndex[j+1];
//      for(int k = start; k < finish; k++)
//      {
//        int column = f_H_mat[ind]->Col[k];
//        dcomplex v1 = h->Value[i];
//        dcomplex v2 = f_H_mat[ind]->Value[k];
//        hash_calc[column].re += v1.re * v2.re - v1.im * v2.im;
//        if(hash[column] == 0)
//        {
//          hash[column] = 1;
//          Qs->Col[countNotZero] = column;
//          countNotZero++;
//        }
//      }
//    }
//    start  = Qs->RowIndex[j];
//    finish = Qs->RowIndex[j+1];
//    for(int k = start; k < finish; k++)
//    {
//      Qs->Value[k] = hash_calc[Qs->Col[k]];
//    }
//  }
//
//  resSum = new crsMatrix(*Qs);
//  Transpose(*Qs, *resSum, false);
//  Transpose(*resSum, *Qs, false);
//
////  Qs->Col[0] = 0;
////  Qs->RowIndex[0] = 0;
////  for(int i = 1; i <= N_mat; i++)
////  {
////    Qs->RowIndex[i] = 0;
////  }
////
////  resSum = new crsMatrix;
////  for(int i = 0; i < h->NZ; i++)
////  {
////    int ind = h->Col[i];
//////    if((h[i].re != 0.0) || (h[i].im != 0.0))
////    {
////      
////      SparseMKLAddT(*Qs, h->Value[i], *(f_mat[ind]), *resSum);
////      tmp = Qs;
////      Qs = resSum;
////      resSum = tmp;
////    }
////  }
//  
//  delete resSum;
//  delete []hash;
//  delete []hash_calc;

//  printMatrixVal(QsCo);
  m->Qs = QsCo;
}


void calcQEs(Model * m)
{
  int N_mat = m->N_mat;
  crsMatrix * he = m->he;
  crsMatrix * QEsCo;
  calc_CooQs(N_mat, m, m->f_ijk, he, QEsCo);

//  crsMatrix ** f_mat = m->f_mat;
//  crsMatrix ** f_H_mat = m->f_H_mat;
//  crsMatrix * QEs;
//  crsMatrix * resSum;
//
//  int rSize = f_H_mat[0]->N;
//  int *hash = new int [rSize];
//  dcomplex *hash_calc = new dcomplex [rSize];
//  
//  QEs = new crsMatrix(rSize, 1);
//  
//  int countNotZero = 0;
//  for(int j = 0;j < rSize; j++)
//  {
//    QEs->RowIndex[j] = countNotZero;
//    memset(hash, 0, sizeof(int) * rSize);
//    for(int i = 0; i < he->NZ; i++)
//    {
//      int ind = he->Col[i];
//      int start, finish;
//      start  = f_H_mat[ind]->RowIndex[j];
//      finish = f_H_mat[ind]->RowIndex[j+1];
//      for(int k = start; k < finish; k++)
//      {
//        if(hash[f_H_mat[ind]->Col[k]] == 0)
//        {
//          hash[f_H_mat[ind]->Col[k]] = 1;
//          countNotZero++;
//        }
//      }
//    }
//  }
//  QEs->RowIndex[rSize] = countNotZero;
//  QEs->setNZ(countNotZero);
//
//  countNotZero = 0;
//  for(int j = 0;j < rSize; j++)
//  {
//    int start, finish;
//    memset(hash_calc, 0, sizeof(dcomplex) * rSize);
//    memset(hash, 0, sizeof(int) * rSize);
//    for(int i = 0; i < he->NZ; i++)
//    {
//      int ind = he->Col[i];
//      start  = f_H_mat[ind]->RowIndex[j];
//      finish = f_H_mat[ind]->RowIndex[j+1];
//      for(int k = start; k < finish; k++)
//      {
//        int column = f_H_mat[ind]->Col[k];
//        dcomplex v1 = he->Value[i];
//        dcomplex v2 = f_H_mat[ind]->Value[k];
//        hash_calc[column].re += v1.re * v2.re - v1.im * v2.im;
//        if(hash[column] == 0)
//        {
//          hash[column] = 1;
//          QEs->Col[countNotZero] = column;
//          countNotZero++;
//        }
//      }
//    }
//    start  = QEs->RowIndex[j];
//    finish = QEs->RowIndex[j+1];
//    for(int k = start; k < finish; k++)
//    {
//      QEs->Value[k] = hash_calc[QEs->Col[k]];
//    }
//  }
//
//  resSum = new crsMatrix(*QEs);
//  Transpose(*QEs, *resSum, false);
//  Transpose(*resSum, *QEs, false);
//
//
//  delete resSum;
//  delete []hash;
//  delete []hash_calc;
// 
//
////  printMatrixVal(QEs);

  m->QEs = QEsCo;
}