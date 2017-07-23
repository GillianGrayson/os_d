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


