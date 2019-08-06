#include "Matrix.h"

void free_matrix(Tensor_Coordinates * matrix)
{
  delete[](matrix->data);
  delete[](matrix->coord1);
  delete[](matrix->coord2);
  delete[](matrix->coord3);
  delete[](matrix->hash);
//  delete matrix;
  matrix->k = 0;
}

Tensor_Coordinates * create_matrix(int NZ)
{
  Tensor_Coordinates * matrix = new Tensor_Coordinates;
  matrix->data = new dcomplex[NZ];
  matrix->coord1 = new unsigned int [NZ];
  matrix->coord2 = new unsigned int [NZ];
  matrix->coord3 = new unsigned int [NZ];
  matrix->hash = new unsigned long long int[NZ];
  matrix->k = NZ;
  return matrix;
}
