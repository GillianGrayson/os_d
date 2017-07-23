#include "Matrix_op.h"
#include <stdio.h>
#include "mkl.h"

// ��������� 2 ���������� ������� � ������� CRS (3 �������, ���������� � ����)
// ���������� C = A * B, C - � ������� CRS (3 �������, ���������� � ����)
//                       ������ ��� C � ������ ��������� �� ����������
// ���������� ������� ���������� ��������: 0 - ��, 1 - �� ��������� ������� (N)
// ���������� ����� ������
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

  // �������� ��������� ��� ������ ������� MKL
  // ��������������� ������� A � B � �������
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

  // ������������ �������, ����������� C = op(A) * B
  char trans;
  
  trans = 'N'; // ������� � ���, op(A) = A - �� ����� ��������������� A

// ������ ��������, �������� �� ��, ��� ����� ���������� ������
// request = 0: ������ ��� �������������� ������� �.�. �������� �������
// ���� �� �� �����, ������� ������ ���������� ��� �������� ����������,
// ����������:
// 1) �������� ������ ��� ������� �������� ����� ic: "���-�� �����+1" ���������;
// 2) ������� ������� � ���������� request = 1 - � ������� ic ����� �������� 
//                                                         ��������� �������
// 3) �������� ������ ��� �������� c � jc 
//    (���-�� ��������� = ic[���-�� �����]-1)
// 4) ������� ������� � ���������� request = 2
  int request;

// ��� ���� ������������� ������: ���� ����������� ���������, ����� �� 
// ������������� ������� A, B � C. � ��� ��������������, ��� ��� �������
// �����������, �������������, �������� ������� "No-No-Yes", �������
// ������������� ������ ��������, ����� ����� ����� �� 1 �� 7 ������������
  int sort = 8;

// ���������� ��������� ���������.
// ������������ ������ ���� request = 0
  int nzmax = -1;

// ��������� ����������
  int info;
  request = 1;
  
  if(!resize)
  {
  // ������� ������ ��� ������� � ������� C
    C.RowIndex = new int[n + 1];
  // ��������� ���������� ��������� ��������� � ������� C
    C.Value = 0;
    C.Col = 0;
  }

  mkl_zcsradd(&trans, &request, &sort, &n, &n,  
              (MKL_Complex16 *)A.Value, A.Col, A.RowIndex,
              (MKL_Complex16 *) (&beta),
              (MKL_Complex16 *)B.Value, B.Col, B.RowIndex, 
              (MKL_Complex16 *)C.Value, C.Col, C.RowIndex, 
              &nzmax, &info);
  int nzc = C.RowIndex[n] - 1;
  if(!resize)
  {
    C.Value = new dcomplex[nzc];
    C.Col = new int[nzc];
  } 
  else
  {
    C.resize(N, nzc);
  }
// ��������� C = A * B
  request = 2;
  mkl_zcsradd(&trans, &request, &sort, &n, &n,  
              (MKL_Complex16 *)A.Value, A.Col, A.RowIndex,
              (MKL_Complex16 *) (&beta),
              (MKL_Complex16 *)B.Value, B.Col, B.RowIndex, 
              (MKL_Complex16 *)C.Value, C.Col, C.RowIndex, 
              &nzmax, &info);
  C.N = n;
  C.NZ = nzc;

  // �������� � ����������� ���� ������� A, B � �
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

// ��������� 2 ���������� ������� � ������� CRS (3 �������, ���������� � ����)
// ���������� C = A * B, C - � ������� CRS (3 �������, ���������� � ����)
//                       ������ ��� C � ������ ��������� �� ����������
// ���������� ������� ���������� ��������: 0 - ��, 1 - �� ��������� ������� (N)
// ���������� ����� ������
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

  // �������� ��������� ��� ������ ������� MKL
  // ��������������� ������� A � B � �������
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

  // ������������ �������, ����������� C = op(A) * B
  char trans;
  
  trans = 'T'; // ������� � ���, op(A) = A - �� ����� ��������������� A

  int request;

  int sort = 8;

  int nzmax = -1;

// ��������� ����������
  int info;
  request = 1;
  
  if(!resize)
  {
  // ������� ������ ��� ������� � ������� C
    C.RowIndex = new int[n + 1];
  // ��������� ���������� ��������� ��������� � ������� C
    C.Value = 0;
    C.Col = 0;
  }

  mkl_zcsradd(&trans, &request, &sort, &n, &n,  
              (MKL_Complex16 *)A.Value, A.Col, A.RowIndex,
              (MKL_Complex16 *) (&beta),
              (MKL_Complex16 *)B.Value, B.Col, B.RowIndex, 
              (MKL_Complex16 *)C.Value, C.Col, C.RowIndex, 
              &nzmax, &info);
  int nzc = C.RowIndex[n] - 1;
  if(!resize)
  {
    C.Value = new dcomplex[nzc];
    C.Col = new int[nzc];
  } 
  else
  {
    C.resize(N, nzc);
  }
// ��������� C = A * B
  request = 2;
  mkl_zcsradd(&trans, &request, &sort, &n, &n,  
              (MKL_Complex16 *)A.Value, A.Col, A.RowIndex,
              (MKL_Complex16 *) (&beta),
              (MKL_Complex16 *)B.Value, B.Col, B.RowIndex, 
              (MKL_Complex16 *)C.Value, C.Col, C.RowIndex, 
              &nzmax, &info);
  C.N = n;
  C.NZ = nzc;

  // �������� � ����������� ���� ������� A, B � �
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


// ��������� 2 ���������� ������� � ������� CRS (3 �������, ���������� � ����)
// ���������� C = A * B, C - � ������� CRS (3 �������, ���������� � ����)
//                       ������ ��� C � ������ ��������� �� ����������
// ���������� ������� ���������� ��������: 0 - ��, 1 - �� ��������� ������� (N)
// ���������� ����� ������
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

  // �������� ��������� ��� ������ ������� MKL
  // ��������������� ������� A � B � �������
  int i, j;

  // ������������ �������, ����������� C = op(A) * B
  char trans = 'N'; // ������� � ���, op(A) = A - �� ����� ��������������� A

// ������ ��������, �������� �� ��, ��� ����� ���������� ������
// request = 0: ������ ��� �������������� ������� �.�. �������� �������
// ���� �� �� �����, ������� ������ ���������� ��� �������� ����������,
// ����������:
// 1) �������� ������ ��� ������� �������� ����� ic: "���-�� �����+1" ���������;
// 2) ������� ������� � ���������� request = 1 - � ������� ic ����� �������� 
//                                                         ��������� �������
// 3) �������� ������ ��� �������� c � jc 
//    (���-�� ��������� = ic[���-�� �����]-1)
// 4) ������� ������� � ���������� request = 2
  int request;

// ��� ���� ������������� ������: ���� ����������� ���������, ����� �� 
// ������������� ������� A, B � C. � ��� ��������������, ��� ��� �������
// �����������, �������������, �������� ������� "No-No-Yes", �������
// ������������� ������ ��������, ����� ����� ����� �� 1 �� 7 ������������
  int sort = 8;

// ���������� ��������� ���������.
// ������������ ������ ���� request = 0
  int nzmax = -1;

// ��������� ����������
  int info;
  request = 1;
  
  if(!resize)
  {
  // ������� ������ ��� ������� � ������� C
    C.RowIndex = new int[n + 1];
  // ��������� ���������� ��������� ��������� � ������� C
    C.Value = 0;
    C.Col = 0;
  }

  mkl_zcsradd(&trans, &request, &sort, &n, &n,  
              (MKL_Complex16 *)A.Value, A.Col, A.RowIndex,
              (MKL_Complex16 *) (&beta),
              (MKL_Complex16 *)B.Value, B.Col, B.RowIndex, 
              (MKL_Complex16 *)C.Value, C.Col, C.RowIndex, 
              &nzmax, &info);
  int nzc = C.RowIndex[n] - 1;
  if(!resize)
  {
    C.Value = new dcomplex[nzc];
    C.Col = new int[nzc];
  } 
  else
  {
    C.resize(N, nzc);
  }
// ��������� C = A * B
  request = 2;
  mkl_zcsradd(&trans, &request, &sort, &n, &n,  
              (MKL_Complex16 *)A.Value, A.Col, A.RowIndex,
              (MKL_Complex16 *) (&beta),
              (MKL_Complex16 *)B.Value, B.Col, B.RowIndex, 
              (MKL_Complex16 *)C.Value, C.Col, C.RowIndex, 
              &nzmax, &info);
  C.N = n;
  C.NZ = nzc;

  // �������� � ����������� ���� ������� A, B � �

  return 0;
}
