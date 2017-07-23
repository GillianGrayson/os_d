#include "coef_coord.h"
#include "stdio.h"
#include "math.h"
#include "string.h"
#include "sortTensor.h"
#include "Matrix.h"

#define IND(i, j, k) ((ulli)i)*((ulli)N * N - 1)*((ulli)N * N - 1) + ((ulli)j)*((ulli)N * N - 1) + ((ulli)k)

#define IndS(i, j) ((int)(((N - 1) + (N - i))/2.0 * (i) + (j - i - 1)) * 2)
#define IndJ(i, j) ((int)(((N - 1) + (N - i))/2.0 * (i) + (j - i - 1)) * 2 + 1)
#define IndD(l)    (N * (N-1) + l - 1)



void fijk_coord_ch(Tensor_Coordinates * f_ijk, crsMatrix *sel, ulli NZ, int N)
{
  ulli ind;
  ulli size = NZ+1;
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
      f_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

      f_ijk->data[f_ijk->k].re = -(i) / sqrt((double)(i)* (i + 1));
      ind = f_ijk->coord1[f_ijk->k] = IndD(i);
      f_ijk->coord2[f_ijk->k] = IndS(i, j);
      f_ijk->coord3[f_ijk->k] = IndJ(i, j);
      f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
      f_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

      f_ijk->data[f_ijk->k].re = -(i) / sqrt((double)(i)* (i + 1));
      ind = f_ijk->coord1[f_ijk->k] = IndJ(i, j);
      f_ijk->coord2[f_ijk->k] = IndD(i);
      f_ijk->coord3[f_ijk->k] = IndS(i, j);
      f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
      f_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

      f_ijk->data[f_ijk->k].re = (i) / sqrt((double)(i)* (i + 1));
      ind = f_ijk->coord1[f_ijk->k] = IndS(i, j);
      f_ijk->coord2[f_ijk->k] = IndD(i);
      f_ijk->coord3[f_ijk->k] = IndJ(i, j);
      f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
      f_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

      f_ijk->data[f_ijk->k].re = (i) / sqrt((double)(i)* (i + 1));
      ind = f_ijk->coord1[f_ijk->k] = IndJ(i, j);
      f_ijk->coord2[f_ijk->k] = IndS(i, j);
      f_ijk->coord3[f_ijk->k] = IndD(i);
      f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
      f_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

      f_ijk->data[f_ijk->k].re = -(i) / sqrt((double)(i)* (i + 1));
      ind = f_ijk->coord1[f_ijk->k] = IndS(i, j);
      f_ijk->coord2[f_ijk->k] = IndJ(i, j);
      f_ijk->coord3[f_ijk->k] = IndD(i);
      f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
      f_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];
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
        f_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

        f_ijk->data[f_ijk->k].re = 1.0 / sqrt((double)(m)* (m + 1));
        ind = f_ijk->coord1[f_ijk->k] = IndD(m);
        f_ijk->coord2[f_ijk->k] = IndS(i, j);
        f_ijk->coord3[f_ijk->k] = IndJ(i, j);
        f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
        f_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

        f_ijk->data[f_ijk->k].re = 1.0 / sqrt((double)(m)* (m + 1));
        ind = f_ijk->coord1[f_ijk->k] = IndJ(i, j);
        f_ijk->coord2[f_ijk->k] = IndD(m);
        f_ijk->coord3[f_ijk->k] = IndS(i, j);
        f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
        f_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

        f_ijk->data[f_ijk->k].re = -1.0 / sqrt((double)(m)* (m + 1));
        ind = f_ijk->coord1[f_ijk->k] = IndS(i, j);
        f_ijk->coord2[f_ijk->k] = IndD(m);
        f_ijk->coord3[f_ijk->k] = IndJ(i, j);
        f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
        f_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

        f_ijk->data[f_ijk->k].re = -1.0 / sqrt((double)(m)* (m + 1));
        ind = f_ijk->coord1[f_ijk->k] = IndJ(i, j);
        f_ijk->coord2[f_ijk->k] = IndS(i, j);
        f_ijk->coord3[f_ijk->k] = IndD(m);
        f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
        f_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

        f_ijk->data[f_ijk->k].re = 1.0 / sqrt((double)(m)* (m + 1));
        ind = f_ijk->coord1[f_ijk->k] = IndS(i, j);
        f_ijk->coord2[f_ijk->k] = IndJ(i, j);
        f_ijk->coord3[f_ijk->k] = IndD(m);
        f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
        f_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];
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
      f_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

      f_ijk->data[f_ijk->k].re = (1 + j) / sqrt((double)(j)* (j + 1));
      ind = f_ijk->coord1[f_ijk->k] = IndD(j);
      f_ijk->coord2[f_ijk->k] = IndS(i, j);
      f_ijk->coord3[f_ijk->k] = IndJ(i, j);
      f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
      f_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

      f_ijk->data[f_ijk->k].re = (1 + j) / sqrt((double)(j)* (j + 1));
      ind = f_ijk->coord1[f_ijk->k] = IndJ(i, j);
      f_ijk->coord2[f_ijk->k] = IndD(j);
      f_ijk->coord3[f_ijk->k] = IndS(i, j);
      f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
      f_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

      f_ijk->data[f_ijk->k].re = -(1 + j) / sqrt((double)(j)* (j + 1));
      ind = f_ijk->coord1[f_ijk->k] = IndS(i, j);
      f_ijk->coord2[f_ijk->k] = IndD(j);
      f_ijk->coord3[f_ijk->k] = IndJ(i, j);
      f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
      f_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

      f_ijk->data[f_ijk->k].re = -(1 + j) / sqrt((double)(j)* (j + 1));
      ind = f_ijk->coord1[f_ijk->k] = IndJ(i, j);
      f_ijk->coord2[f_ijk->k] = IndS(i, j);
      f_ijk->coord3[f_ijk->k] = IndD(j);
      f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
      f_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

      f_ijk->data[f_ijk->k].re = (1 + j) / sqrt((double)(j)* (j + 1));
      ind = f_ijk->coord1[f_ijk->k] = IndS(i, j);
      f_ijk->coord2[f_ijk->k] = IndJ(i, j);
      f_ijk->coord3[f_ijk->k] = IndD(j);
      f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
      f_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];
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
            f_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

            f_ijk->data[f_ijk->k].re = 1.0 / sqrt(2.0); //-i*Jik*[Sij,Sjk]  //and symmetric Skj, Jki
            ind = f_ijk->coord1[f_ijk->k] = IndJ(i, k);
            f_ijk->coord2[f_ijk->k] = IndS(i, j);
            f_ijk->coord3[f_ijk->k] = IndS(j, k);
            f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
            f_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

            f_ijk->data[f_ijk->k].re = -1.0 / sqrt(2.0); //-i*Sjk*[Sij,Jik]  //and symmetric Skj, Jki
            ind = f_ijk->coord1[f_ijk->k] = IndS(j, k);
            f_ijk->coord2[f_ijk->k] = IndS(i, j);
            f_ijk->coord3[f_ijk->k] = IndJ(i, k);
            f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
            f_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

            f_ijk->data[f_ijk->k].re = -1.0 / sqrt(2.0); //-i*Sik*[Sij,Jjk]  //and symmetric Ski, Jkj
            ind = f_ijk->coord1[f_ijk->k] = IndS(i, k);
            f_ijk->coord2[f_ijk->k] = IndS(i, j);
            f_ijk->coord3[f_ijk->k] = IndJ(j, k);
            f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
            f_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

            f_ijk->data[f_ijk->k].re = 1.0 / sqrt(2.0); //-i*Sjk*[Jij,Sik]  //and symmetric Ski, Skj
            ind = f_ijk->coord1[f_ijk->k] = IndS(j, k);
            f_ijk->coord2[f_ijk->k] = IndJ(i, j);
            f_ijk->coord3[f_ijk->k] = IndS(i, k);
            f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
            f_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

            f_ijk->data[f_ijk->k].re = -1.0 / sqrt(2.0); //-i*Sik*[Jij,Sjk]  //and symmetric Skj, Ski
            ind = f_ijk->coord1[f_ijk->k] = IndS(i, k);
            f_ijk->coord2[f_ijk->k] = IndJ(i, j);
            f_ijk->coord3[f_ijk->k] = IndS(j, k);
            f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
            f_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

            f_ijk->data[f_ijk->k].re = 1.0 / sqrt(2.0); //-i*Jjk*[Jij,Jik]  //and symmetric Jkj, Jki
            ind = f_ijk->coord1[f_ijk->k] = IndJ(j, k);
            f_ijk->coord2[f_ijk->k] = IndJ(i, j);
            f_ijk->coord3[f_ijk->k] = IndJ(i, k);
            f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
            f_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

            f_ijk->data[f_ijk->k].re = -1.0 / sqrt(2.0); //-i*Jik*[Jij,Jjk]  //and symmetric Jki, Jkj
            ind = f_ijk->coord1[f_ijk->k] = IndJ(i, k);
            f_ijk->coord2[f_ijk->k] = IndJ(i, j);
            f_ijk->coord3[f_ijk->k] = IndJ(j, k);
            f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
            f_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];
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
            f_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

            f_ijk->data[f_ijk->k].re = 1.0 / sqrt(2.0);
            ind = f_ijk->coord1[f_ijk->k] = IndJ(i, k);
            f_ijk->coord2[f_ijk->k] = IndS(i, j);
            f_ijk->coord3[f_ijk->k] = IndS(k, j);
            f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
            f_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

            f_ijk->data[f_ijk->k].re = -1.0 / sqrt(2.0);
            ind = f_ijk->coord1[f_ijk->k] = IndS(k, j);
            f_ijk->coord2[f_ijk->k] = IndS(i, j);
            f_ijk->coord3[f_ijk->k] = IndJ(i, k);
            f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
            f_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

            f_ijk->data[f_ijk->k].re = 1.0 / sqrt(2.0);
            ind = f_ijk->coord1[f_ijk->k] = IndS(i, k);
            f_ijk->coord2[f_ijk->k] = IndS(i, j);
            f_ijk->coord3[f_ijk->k] = IndJ(k, j);
            f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
            f_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

            f_ijk->data[f_ijk->k].re = 1.0 / sqrt(2.0);
            ind = f_ijk->coord1[f_ijk->k] = IndS(k, j);
            f_ijk->coord2[f_ijk->k] = IndJ(i, j);
            f_ijk->coord3[f_ijk->k] = IndS(i, k);
            f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
            f_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

            f_ijk->data[f_ijk->k].re = -1.0 / sqrt(2.0);
            ind = f_ijk->coord1[f_ijk->k] = IndS(i, k);
            f_ijk->coord2[f_ijk->k] = IndJ(i, j);
            f_ijk->coord3[f_ijk->k] = IndS(k, j);
            f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
            f_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

            f_ijk->data[f_ijk->k].re = -1.0 / sqrt(2.0);
            ind = f_ijk->coord1[f_ijk->k] = IndJ(k, j);
            f_ijk->coord2[f_ijk->k] = IndJ(i, j);
            f_ijk->coord3[f_ijk->k] = IndJ(i, k);
            f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
            f_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

            f_ijk->data[f_ijk->k].re = 1.0 / sqrt(2.0);
            ind = f_ijk->coord1[f_ijk->k] = IndJ(i, k);
            f_ijk->coord2[f_ijk->k] = IndJ(i, j);
            f_ijk->coord3[f_ijk->k] = IndJ(k, j);
            f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
            f_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];
          }
          if (k < i)
          {
            f_ijk->data[f_ijk->k].re = -1.0 / sqrt(2.0);
            ind = f_ijk->coord1[f_ijk->k] = IndJ(k, j);
            f_ijk->coord2[f_ijk->k] = IndS(i, j);
            f_ijk->coord3[f_ijk->k] = IndS(k, i);
            f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
            f_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

            f_ijk->data[f_ijk->k].re = -1.0 / sqrt(2.0);
            ind = f_ijk->coord1[f_ijk->k] = IndJ(k, i);
            f_ijk->coord2[f_ijk->k] = IndS(i, j);
            f_ijk->coord3[f_ijk->k] = IndS(k, j);
            f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
            f_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

            f_ijk->data[f_ijk->k].re = 1.0 / sqrt(2.0);
            ind = f_ijk->coord1[f_ijk->k] = IndS(k, j);
            f_ijk->coord2[f_ijk->k] = IndS(i, j);
            f_ijk->coord3[f_ijk->k] = IndJ(k, i);
            f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
            f_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

            f_ijk->data[f_ijk->k].re = 1.0 / sqrt(2.0);
            ind = f_ijk->coord1[f_ijk->k] = IndS(k, i);
            f_ijk->coord2[f_ijk->k] = IndS(i, j);
            f_ijk->coord3[f_ijk->k] = IndJ(k, j);
            f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
            f_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

            f_ijk->data[f_ijk->k].re = 1.0 / sqrt(2.0);
            ind = f_ijk->coord1[f_ijk->k] = IndS(k, j);
            f_ijk->coord2[f_ijk->k] = IndJ(i, j);
            f_ijk->coord3[f_ijk->k] = IndS(k, i);
            f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
            f_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

            f_ijk->data[f_ijk->k].re = -1.0 / sqrt(2.0);
            ind = f_ijk->coord1[f_ijk->k] = IndS(k, i);
            f_ijk->coord2[f_ijk->k] = IndJ(i, j);
            f_ijk->coord3[f_ijk->k] = IndS(k, j);
            f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
            f_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

            f_ijk->data[f_ijk->k].re = 1.0 / sqrt(2.0);
            ind = f_ijk->coord1[f_ijk->k] = IndJ(k, j);
            f_ijk->coord2[f_ijk->k] = IndJ(i, j);
            f_ijk->coord3[f_ijk->k] = IndJ(k, i);
            f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
            f_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

            f_ijk->data[f_ijk->k].re = -1.0 / sqrt(2.0);
            ind = f_ijk->coord1[f_ijk->k] = IndJ(k, i);
            f_ijk->coord2[f_ijk->k] = IndJ(i, j);
            f_ijk->coord3[f_ijk->k] = IndJ(k, j);
            f_ijk->hash[f_ijk->k] = f_ijk->coord1[f_ijk->k] + f_ijk->coord2[f_ijk->k] * step2 + step1 * f_ijk->coord3[f_ijk->k];
            f_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];
          }
        }
      }
    }
  }
}

void dijk_coord_ch(Tensor_Coordinates * d_ijk, crsMatrix *sel, ulli NZ, int N)
{
  ulli ind;
  ulli size = NZ+1;
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
      d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

      d_ijk->data[d_ijk->k].re = -(i) / sqrt((double)(i)* (i + 1));
      ind = d_ijk->coord1[d_ijk->k] = IndS(i, j);
      d_ijk->coord2[d_ijk->k] = IndD(i);
      d_ijk->coord3[d_ijk->k] = IndS(i, j);
      d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
      d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

      d_ijk->data[d_ijk->k].re = -(i) / sqrt((double)(i)* (i + 1));
      ind = d_ijk->coord1[d_ijk->k] = IndD(i);
      d_ijk->coord2[d_ijk->k] = IndJ(i, j);
      d_ijk->coord3[d_ijk->k] = IndJ(i, j);
      d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
      d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

      d_ijk->data[d_ijk->k].re = -(i) / sqrt((double)(i)* (i + 1));
      ind = d_ijk->coord1[d_ijk->k] = IndJ(i, j);
      d_ijk->coord2[d_ijk->k] = IndD(i);
      d_ijk->coord3[d_ijk->k] = IndJ(i, j);
      d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
      d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

      d_ijk->data[d_ijk->k].re = -(i) / sqrt((double)(i)* (i + 1));
      ind = d_ijk->coord1[d_ijk->k] = IndS(i, j);
      d_ijk->coord2[d_ijk->k] = IndS(i, j);
      d_ijk->coord3[d_ijk->k] = IndD(i);
      d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
      d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

      d_ijk->data[d_ijk->k].re = -(i) / sqrt((double)(i)* (i + 1));
      ind = d_ijk->coord1[d_ijk->k] = IndJ(i, j);
      d_ijk->coord2[d_ijk->k] = IndJ(i, j);
      d_ijk->coord3[d_ijk->k] = IndD(i);
      d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
      d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];
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
        d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

        d_ijk->data[d_ijk->k].re = 1.0 / sqrt((double)(m)* (m + 1));
        ind = d_ijk->coord1[d_ijk->k] = IndS(i, j);
        d_ijk->coord2[d_ijk->k] = IndD(m);
        d_ijk->coord3[d_ijk->k] = IndS(i, j);
        d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
        d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

        d_ijk->data[d_ijk->k].re = 1.0 / sqrt((double)(m)* (m + 1));
        ind = d_ijk->coord1[d_ijk->k] = IndD(m);
        d_ijk->coord2[d_ijk->k] = IndJ(i, j);
        d_ijk->coord3[d_ijk->k] = IndJ(i, j);
        d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
        d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

        d_ijk->data[d_ijk->k].re = 1.0 / sqrt((double)(m)* (m + 1));
        ind = d_ijk->coord1[d_ijk->k] = IndJ(i, j);
        d_ijk->coord2[d_ijk->k] = IndD(m);
        d_ijk->coord3[d_ijk->k] = IndJ(i, j);
        d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
        d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

        d_ijk->data[d_ijk->k].re = 1.0 / sqrt((double)(m)* (m + 1));
        ind = d_ijk->coord1[d_ijk->k] = IndS(i, j);
        d_ijk->coord2[d_ijk->k] = IndS(i, j);
        d_ijk->coord3[d_ijk->k] = IndD(m);
        d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
        d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

        d_ijk->data[d_ijk->k].re = 1.0 / sqrt((double)(m)* (m + 1));
        ind = d_ijk->coord1[d_ijk->k] = IndJ(i, j);
        d_ijk->coord2[d_ijk->k] = IndJ(i, j);
        d_ijk->coord3[d_ijk->k] = IndD(m);
        d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
        d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];
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
    d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

    d_ijk->data[d_ijk->k].re = (1 - j) / sqrt((double)(j)* (j + 1));
    ind = d_ijk->coord1[d_ijk->k] = IndS(i, j);
    d_ijk->coord2[d_ijk->k] = IndD(j);
    d_ijk->coord3[d_ijk->k] = IndS(i, j);
    d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
    d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

    d_ijk->data[d_ijk->k].re = (1 - j) / sqrt((double)(j)* (j + 1));
    ind = d_ijk->coord1[d_ijk->k] = IndD(j);
    d_ijk->coord2[d_ijk->k] = IndJ(i, j);
    d_ijk->coord3[d_ijk->k] = IndJ(i, j);
    d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
    d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

    d_ijk->data[d_ijk->k].re = (1 - j) / sqrt((double)(j)* (j + 1));
    ind = d_ijk->coord1[d_ijk->k] = IndJ(i, j);
    d_ijk->coord2[d_ijk->k] = IndD(j);
    d_ijk->coord3[d_ijk->k] = IndJ(i, j);
    d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
    d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

    d_ijk->data[d_ijk->k].re = (1 - j) / sqrt((double)(j)* (j + 1));
    ind = d_ijk->coord1[d_ijk->k] = IndS(i, j);
    d_ijk->coord2[d_ijk->k] = IndS(i, j);
    d_ijk->coord3[d_ijk->k] = IndD(j);
    d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
    d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

    d_ijk->data[d_ijk->k].re = (1 - j) / sqrt((double)(j)* (j + 1));
    ind = d_ijk->coord1[d_ijk->k] = IndJ(i, j);
    d_ijk->coord2[d_ijk->k] = IndJ(i, j);
    d_ijk->coord3[d_ijk->k] = IndD(j);
    d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
    d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];
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
      d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

      d_ijk->data[d_ijk->k].re = (1 - j) / sqrt((double)(j)* (j + 1));
      ind = d_ijk->coord1[d_ijk->k] = IndS(i, j);
      d_ijk->coord2[d_ijk->k] = IndD(j);
      d_ijk->coord3[d_ijk->k] = IndS(i, j);
      d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
      d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

      d_ijk->data[d_ijk->k].re = (1 - j) / sqrt((double)(j)* (j + 1));
      ind = d_ijk->coord1[d_ijk->k] = IndD(j);
      d_ijk->coord2[d_ijk->k] = IndJ(i, j);
      d_ijk->coord3[d_ijk->k] = IndJ(i, j);
      d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
      d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

      d_ijk->data[d_ijk->k].re = (1 - j) / sqrt((double)(j)* (j + 1));
      ind = d_ijk->coord1[d_ijk->k] = IndJ(i, j);
      d_ijk->coord2[d_ijk->k] = IndD(j);
      d_ijk->coord3[d_ijk->k] = IndJ(i, j);
      d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
      d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

      d_ijk->data[d_ijk->k].re = (1 - j) / sqrt((double)(j)* (j + 1));
      ind = d_ijk->coord1[d_ijk->k] = IndS(i, j);
      d_ijk->coord2[d_ijk->k] = IndS(i, j);
      d_ijk->coord3[d_ijk->k] = IndD(j);
      d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
      d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

      d_ijk->data[d_ijk->k].re = (1 - j) / sqrt((double)(j)* (j + 1));
      ind = d_ijk->coord1[d_ijk->k] = IndJ(i, j);
      d_ijk->coord2[d_ijk->k] = IndJ(i, j);
      d_ijk->coord3[d_ijk->k] = IndD(j);
      d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
      d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];
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
        d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

        d_ijk->data[d_ijk->k].re = 2.0 / sqrt((double)(z)* (z + 1));
        ind = d_ijk->coord1[d_ijk->k] = IndS(i, j);
        d_ijk->coord2[d_ijk->k] = IndD(z);
        d_ijk->coord3[d_ijk->k] = IndS(i, j);
        d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
        d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

        d_ijk->data[d_ijk->k].re = 2.0 / sqrt((double)(z)* (z + 1));
        ind = d_ijk->coord1[d_ijk->k] = IndD(z);
        d_ijk->coord2[d_ijk->k] = IndJ(i, j);
        d_ijk->coord3[d_ijk->k] = IndJ(i, j);
        d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
        d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

        d_ijk->data[d_ijk->k].re = 2.0 / sqrt((double)(z)* (z + 1));
        ind = d_ijk->coord1[d_ijk->k] = IndJ(i, j);
        d_ijk->coord2[d_ijk->k] = IndD(z);
        d_ijk->coord3[d_ijk->k] = IndJ(i, j);
        d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
        d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

        d_ijk->data[d_ijk->k].re = 2.0 / sqrt((double)(z)* (z + 1));
        ind = d_ijk->coord1[d_ijk->k] = IndS(i, j);
        d_ijk->coord2[d_ijk->k] = IndS(i, j);
        d_ijk->coord3[d_ijk->k] = IndD(z);
        d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
        d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

        d_ijk->data[d_ijk->k].re = 2.0 / sqrt((double)(z)* (z + 1));
        ind = d_ijk->coord1[d_ijk->k] = IndJ(i, j);
        d_ijk->coord2[d_ijk->k] = IndJ(i, j);
        d_ijk->coord3[d_ijk->k] = IndD(z);
        d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
        d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];
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
            d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

            d_ijk->data[d_ijk->k].re = 1.0 / sqrt(2.0);  //Sik*{Sij,Sjk}  //and symmetric Skj, Ski
            ind = d_ijk->coord1[d_ijk->k] = IndS(i, k);
            d_ijk->coord2[d_ijk->k] = IndS(i, j);
            d_ijk->coord3[d_ijk->k] = IndS(j, k);
            d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
            d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

            d_ijk->data[d_ijk->k].re = 1.0 / sqrt(2.0);  //Jjk*{Sij,Jik}  //and symmetric Jkj, Jki
            ind = d_ijk->coord1[d_ijk->k] = IndJ(j, k);
            d_ijk->coord2[d_ijk->k] = IndS(i, j);
            d_ijk->coord3[d_ijk->k] = IndJ(i, k);
            d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
            d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

            d_ijk->data[d_ijk->k].re = 1.0 / sqrt(2.0);  //Jik*{Sij,Jjk}  //and symmetric Jki, Jkj
            ind = d_ijk->coord1[d_ijk->k] = IndJ(i, k);
            d_ijk->coord2[d_ijk->k] = IndS(i, j);
            d_ijk->coord3[d_ijk->k] = IndJ(j, k);
            d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
            d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

            d_ijk->data[d_ijk->k].re = -1.0 / sqrt(2.0);  //Jjk*{Jij,Sik}  //and symmetric Ski, Jkj
            ind = d_ijk->coord1[d_ijk->k] = IndJ(j, k);
            d_ijk->coord2[d_ijk->k] = IndJ(i, j);
            d_ijk->coord3[d_ijk->k] = IndS(i, k);
            d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
            d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

            d_ijk->data[d_ijk->k].re = 1.0 / sqrt(2.0);  //Jik*{Jij,Sjk}  //and symmetric Skj, Jki
            ind = d_ijk->coord1[d_ijk->k] = IndJ(i, k);
            d_ijk->coord2[d_ijk->k] = IndJ(i, j);
            d_ijk->coord3[d_ijk->k] = IndS(j, k);
            d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
            d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

            d_ijk->data[d_ijk->k].re = 1.0 / sqrt(2.0);  //Sjk*{Jij,Jik}  //and symmetric Skj, Jki
            ind = d_ijk->coord1[d_ijk->k] = IndS(j, k);
            d_ijk->coord2[d_ijk->k] = IndJ(i, j);
            d_ijk->coord3[d_ijk->k] = IndJ(i, k);
            d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
            d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

            d_ijk->data[d_ijk->k].re = -1.0 / sqrt(2.0);  //Sik*{Jij,Jjk}  //and symmetric Ski, Jkj
            ind = d_ijk->coord1[d_ijk->k] = IndS(i, k);
            d_ijk->coord2[d_ijk->k] = IndJ(i, j);
            d_ijk->coord3[d_ijk->k] = IndJ(j, k);
            d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
            d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];
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
            d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

            d_ijk->data[d_ijk->k].re = 1.0 / sqrt(2.0);
            ind = d_ijk->coord1[d_ijk->k] = IndS(i, k);
            d_ijk->coord2[d_ijk->k] = IndS(i, j);
            d_ijk->coord3[d_ijk->k] = IndS(k, j);
            d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
            d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

            d_ijk->data[d_ijk->k].re = -1.0 / sqrt(2.0);
            ind = d_ijk->coord1[d_ijk->k] = IndJ(k, j);
            d_ijk->coord2[d_ijk->k] = IndS(i, j);
            d_ijk->coord3[d_ijk->k] = IndJ(i, k);
            d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
            d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

            d_ijk->data[d_ijk->k].re = -1.0 / sqrt(2.0);
            ind = d_ijk->coord1[d_ijk->k] = IndJ(i, k);
            d_ijk->coord2[d_ijk->k] = IndS(i, j);
            d_ijk->coord3[d_ijk->k] = IndJ(k, j);
            d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
            d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

            d_ijk->data[d_ijk->k].re = 1.0 / sqrt(2.0);
            ind = d_ijk->coord1[d_ijk->k] = IndJ(k, j);
            d_ijk->coord2[d_ijk->k] = IndJ(i, j);
            d_ijk->coord3[d_ijk->k] = IndS(i, k);
            d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
            d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

            d_ijk->data[d_ijk->k].re = 1.0 / sqrt(2.0);
            ind = d_ijk->coord1[d_ijk->k] = IndJ(i, k);
            d_ijk->coord2[d_ijk->k] = IndJ(i, j);
            d_ijk->coord3[d_ijk->k] = IndS(k, j);
            d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
            d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

            d_ijk->data[d_ijk->k].re = 1.0 / sqrt(2.0);
            ind = d_ijk->coord1[d_ijk->k] = IndS(k, j);
            d_ijk->coord2[d_ijk->k] = IndJ(i, j);
            d_ijk->coord3[d_ijk->k] = IndJ(i, k);
            d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
            d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

            d_ijk->data[d_ijk->k].re = 1.0 / sqrt(2.0);
            ind = d_ijk->coord1[d_ijk->k] = IndS(i, k);
            d_ijk->coord2[d_ijk->k] = IndJ(i, j);
            d_ijk->coord3[d_ijk->k] = IndJ(k, j);
            d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
            d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];
          }
          if (k < i)
          {
            d_ijk->data[d_ijk->k].re = 1.0 / sqrt(2.0);
            ind = d_ijk->coord1[d_ijk->k] = IndS(k, j);
            d_ijk->coord2[d_ijk->k] = IndS(i, j);
            d_ijk->coord3[d_ijk->k] = IndS(k, i);
            d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
            d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

            d_ijk->data[d_ijk->k].re = 1.0 / sqrt(2.0);
            ind = d_ijk->coord1[d_ijk->k] = IndS(k, i);
            d_ijk->coord2[d_ijk->k] = IndS(i, j);
            d_ijk->coord3[d_ijk->k] = IndS(k, j);
            d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
            d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

            d_ijk->data[d_ijk->k].re = 1.0 / sqrt(2.0);
            ind = d_ijk->coord1[d_ijk->k] = IndJ(k, j);
            d_ijk->coord2[d_ijk->k] = IndS(i, j);
            d_ijk->coord3[d_ijk->k] = IndJ(k, i);
            d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
            d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

            d_ijk->data[d_ijk->k].re = 1.0 / sqrt(2.0);
            ind = d_ijk->coord1[d_ijk->k] = IndJ(k, i);
            d_ijk->coord2[d_ijk->k] = IndS(i, j);
            d_ijk->coord3[d_ijk->k] = IndJ(k, j);
            d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
            d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

            d_ijk->data[d_ijk->k].re = 1.0 / sqrt(2.0);
            ind = d_ijk->coord1[d_ijk->k] = IndJ(k, j);
            d_ijk->coord2[d_ijk->k] = IndJ(i, j);
            d_ijk->coord3[d_ijk->k] = IndS(k, i);
            d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
            d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

            d_ijk->data[d_ijk->k].re = -1.0 / sqrt(2.0);
            ind = d_ijk->coord1[d_ijk->k] = IndJ(k, i);
            d_ijk->coord2[d_ijk->k] = IndJ(i, j);
            d_ijk->coord3[d_ijk->k] = IndS(k, j);
            d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
            d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

            d_ijk->data[d_ijk->k].re = -1.0 / sqrt(2.0);
            ind = d_ijk->coord1[d_ijk->k] = IndS(k, j);
            d_ijk->coord2[d_ijk->k] = IndJ(i, j);
            d_ijk->coord3[d_ijk->k] = IndJ(k, i);
            d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
            d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

            d_ijk->data[d_ijk->k].re = 1.0 / sqrt(2.0);
            ind = d_ijk->coord1[d_ijk->k] = IndS(k, i);
            d_ijk->coord2[d_ijk->k] = IndJ(i, j);
            d_ijk->coord3[d_ijk->k] = IndJ(k, j);
            d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
            d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];
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
    d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

    d_ijk->data[d_ijk->k].re = 2.0 / sqrt((double)(j)* (j + 1));
    ind = d_ijk->coord1[d_ijk->k] = IndD(i);
    d_ijk->coord2[d_ijk->k] = IndD(i);
    d_ijk->coord3[d_ijk->k] = IndD(j);
    d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
    d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

    d_ijk->data[d_ijk->k].re = 2.0 / sqrt((double)(j)* (j + 1));
    ind = d_ijk->coord1[d_ijk->k] = IndD(i);
    d_ijk->coord2[d_ijk->k] = IndD(j);
    d_ijk->coord3[d_ijk->k] = IndD(i);
    d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
    d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];
  }

  for (int i = 2; i < N; i++)
  {
    d_ijk->data[d_ijk->k].re = 2.0 * (1 - i) / sqrt((double)(i)* (i + 1));
    ind = d_ijk->coord1[d_ijk->k] = IndD(i);
    d_ijk->coord2[d_ijk->k] = IndD(i);
    d_ijk->coord3[d_ijk->k] = IndD(i);
    d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
    d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

    for (int j = i + 1; j < N; j++)
    {
      d_ijk->data[d_ijk->k].re = 2.0 / sqrt((double)(j)* (j + 1));
      ind = d_ijk->coord1[d_ijk->k] = IndD(j);
      d_ijk->coord2[d_ijk->k] = IndD(i);
      d_ijk->coord3[d_ijk->k] = IndD(i);
      d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
      d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

      d_ijk->data[d_ijk->k].re = 2.0 / sqrt((double)(j)* (j + 1));
      ind = d_ijk->coord1[d_ijk->k] = IndD(i);
      d_ijk->coord2[d_ijk->k] = IndD(i);
      d_ijk->coord3[d_ijk->k] = IndD(j);
      d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
      d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];

      d_ijk->data[d_ijk->k].re = 2.0 / sqrt((double)(j)* (j + 1));
      ind = d_ijk->coord1[d_ijk->k] = IndD(i);
      d_ijk->coord2[d_ijk->k] = IndD(j);
      d_ijk->coord3[d_ijk->k] = IndD(i);
      d_ijk->hash[d_ijk->k] = d_ijk->coord1[d_ijk->k] + d_ijk->coord2[d_ijk->k] * step2 + step1 * d_ijk->coord3[d_ijk->k];
      d_ijk->k+= sel->RowIndex[ind + 1] - sel->RowIndex[ind];
    }
  }
}
