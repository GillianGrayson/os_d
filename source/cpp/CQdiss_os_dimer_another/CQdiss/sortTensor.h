#ifndef __SORT_TENSOR__
#define __SORT_TENSOR__

#include "Matrix.h"

ulli max_bit(ulli * mas, unsigned int n);

inline void swap(Tensor_Coordinates  * mas, unsigned int i, unsigned int j);

void msd_sort(Tensor_Coordinates * mas, unsigned int from, unsigned int to, ulli bit, int threads_level);

void sort_matrix(Tensor_Coordinates * matrix);

#endif