#include "sortTensor.h"

ulli max_bit(ulli * mas, unsigned int n)
{
  ulli mask = 0;
  //#pragma omp parallel for  reduction(|:mask)
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

