#ifndef __MATRIX__
#define __MATRIX__

#include <stdio.h>
// Complex (double is the base datatype for real and imaginary parts)
typedef struct dcomplex
{
  double re;
  double im;
} dcomplex;

typedef unsigned long long int ulli;

struct Tensor_Coordinates
{
  dcomplex * data;
  unsigned int * coord1;
  unsigned int * coord2;
  unsigned int * coord3;
  ulli * hash;
  ulli k;
};

struct Tensor_Coordinates_1
{
  dcomplex * data;
  unsigned int coord1;
  unsigned int * coord2;
  unsigned int * coord3;
  unsigned int N;
};

void free_matrix(Tensor_Coordinates * matrix);
Tensor_Coordinates * create_matrix(int NZ);

struct crsMatrix
{
protected:
  int sizeMem;

public:
  int N;  // Размер матрицы (N x N)
  int NZ; // Кол-во ненулевых элементов

  // Массив значений (размер NZ)
  dcomplex* Value;
  // Массив номеров столбцов (размер NZ)
  int* Col;
  // Массив индексов строк (размер N + 1)
  int* RowIndex;

  crsMatrix()
  {
    N = 0;
    sizeMem = NZ = 0;
    Value = NULL;
    Col = NULL;
    RowIndex = NULL;
  }

  void resize(int _N, int _NZ=0)
  {
    if(N < _N)
    {
      if(RowIndex != NULL)
      {
        delete [] RowIndex;
      }
      N = _N;
      RowIndex = new int[N + 1];
    }
    N = _N;
      
    NZ = _NZ;
    if(_NZ > sizeMem)
    {
      sizeMem = NZ;
      if(Col != NULL)
      {
        delete[] Col;
      }
      Col = new int[NZ];

      if(Value != NULL)
      {
        delete[] Value;
      }
      Value = new dcomplex[NZ];
    }
  }

  crsMatrix(int _N, int _NZ)
  {
    N = _N;
    sizeMem = NZ = _NZ;
    
    Value = new dcomplex[NZ];
    Col = new int[NZ];
    RowIndex = new int[N + 1];
    
    for(int i = 0; i < NZ; i++)
    {
      Value[i].re = 0.0;
      Value[i].im = 0.0;
    }
  }

  void setNZ(int _NZ)
  {
    sizeMem = NZ = _NZ;
    
    if (Value != NULL)
    {
      delete[] Value;
    }
    
    if (Col != NULL)
    {
      delete[] Col;
    }

    Value = new dcomplex[NZ];
    Col = new int[NZ];
    for(int i = 0; i < NZ; i++)
    {
      Value[i].re = 0.0;
      Value[i].im = 0.0;
    }
  }

  crsMatrix(const crsMatrix &mat)
  {
    N = mat.N;
    sizeMem = NZ = mat.NZ;
    
    Value = new dcomplex[NZ];
    Col = new int[NZ];
    RowIndex = new int[N + 1];

    for(int i = 0; i < N + 1; i++)
    {
      RowIndex[i] = mat.RowIndex[i];
    }

    for(int i = 0; i < NZ; i++)
    {
      Value[i].re = mat.Value[i].re;
      Value[i].im = mat.Value[i].im;
      Col[i] = mat.Col[i];
    }

  }

  ~crsMatrix()
  {
    delete[] Value;
    delete[] Col;
    delete[] RowIndex;
  }
};

#endif