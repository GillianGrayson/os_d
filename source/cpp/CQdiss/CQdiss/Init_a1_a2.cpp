#include "Init_a1_a2.h"
#include "initH.h"
#include <math.h>

crsMatrix * createA1mat(int N)
{
  crsMatrix * mat;
  mat  = new crsMatrix(N + 1, N + 1);
  mat->RowIndex[0] = 0;

  for(int i = 0; i < mat->N; i++)
  {
    mat->Col[i] = i;
    mat->Value[i].re = 2.0 * i - N;
    mat->RowIndex[i + 1] = mat->RowIndex[i] + 1;
  }

  return mat;
}

crsMatrix * createA2mat(int N)
{
  int i, j;
  crsMatrix * mat;
  mat  = new crsMatrix(N + 1, N * 2);
  mat->RowIndex[0] = 0;
  mat->Col[0] = 1;
  mat->RowIndex[1] = 1;
  for(i = 1; i < N; i++)
    mat->RowIndex[i + 1] = mat->RowIndex[i] + 2;
  mat->RowIndex[i + 1] = mat->RowIndex[i] + 1;
  double val;
  val = N;
  mat->Value[0].im = sqrt(val);
  for(i = 1; i < N; i++)
  {
    j = mat->RowIndex[i];
    mat->Col[j] = i-1;
    mat->Col[j+1] = i+1;
    val = (N - i + 1)  * (i + 0);
    mat->Value[j].im = -sqrt(val);
    val = (N - i + 0)  * (i + 1);
    mat->Value[j+1].im = sqrt(val);
  }
  j = mat->RowIndex[i];
  mat->Col[j] = i - 1;
  val = (N - i + 1)  * (i + 0);
  mat->Value[j].im = -sqrt(val);

  return mat;
}

void init_a1_a2_opt(Model * m)
{
	int N = m->N;
	FMatrixs *Fs = m->Fs;

	crsMatrix * A1 = createA1mat(N);
	crsMatrix * A2 = createA2mat(N);

	int N_mat = m->N_mat;

	crsMatrix * a1_mat = new crsMatrix(N_mat, N_mat);
	crsMatrix * a2_mat = new crsMatrix(N_mat, N_mat);

	to_F_basis(A1, a1_mat);
	to_F_basis(A2, a2_mat);
	
	crsMatrix * a1_i_a2_mat = new crsMatrix(N_mat, a2_mat->NZ + a1_mat->NZ);
	int k = 0;
	int k1, k2;
	for (int i = 0; i < N_mat; i++)
	{
		a1_i_a2_mat->RowIndex[i] = k;
		int c1 = a1_mat->RowIndex[i + 1] - a1_mat->RowIndex[i];
		int c2 = a2_mat->RowIndex[i + 1] - a2_mat->RowIndex[i];
		//    a_res->RowIndex[i] = k * a2_mat->NZ;
		if (c1 > 0)
		{
			a1_i_a2_mat->Value[k] = a1_mat->Value[a1_mat->RowIndex[i]];
		}
		if (c2 > 0)
		{
			a1_i_a2_mat->Value[k].re += a2_mat->Value[a2_mat->RowIndex[i]].im;
			a1_i_a2_mat->Value[k].im -= a2_mat->Value[a2_mat->RowIndex[i]].re;
		}
		if (c1 + c2 > 0)
		{
			a1_i_a2_mat->Col[k] = i;
			k++;
		}
	}
	a1_i_a2_mat->RowIndex[N_mat] = k;
	a1_i_a2_mat->NZ = k;
//	printMatrixVal(a2_mat);
//	  printMatrixVal(a1_i_a2_mat);
	m->l_mat = a1_i_a2_mat;
	delete a1_mat;
	delete a2_mat;

	delete A1;
	delete A2;
}

void init_a1_a2(Model * m)
{
  int N = m->N;
  FMatrixs *Fs = m->Fs;

  crsMatrix * A1 = createA1mat(N);
//  printMatrix(A1);

  crsMatrix * A2 = createA2mat(N);
//  printMatrix(A2);

  int N_mat = (N + 1) * (N + 1) - 1;

  crsMatrix * a1_mat = new crsMatrix(N_mat, N_mat);
  crsMatrix * a2_mat = new crsMatrix(N_mat, N_mat);

  int k = 0;

  crsMatrix * res;
  int cnt;
  a1_mat->RowIndex[0] = 0;
  for(int i = 0; i < N_mat; i++)
  {
    res = new crsMatrix;
    SparseMKLMult(*A1, *(Fs->F[i+1]), *res);
    cnt = trace_struct(*res);
    if(cnt > 0)
    {
      a1_mat->Value[k] = trace(*res);
      a1_mat->Col[k] = i;
      k++;
    }
    a1_mat->RowIndex[i + 1] = k;
    delete res;
  }
  a1_mat->NZ = k;

  //printMatrixVal(a1_mat);

  k = 0;
  a2_mat->RowIndex[0] = 0;
  for(int i = 0; i < N_mat; i++)
  {
    res = new crsMatrix;
    SparseMKLMult(*A2, *(Fs->F[i+1]), *res);
    cnt = trace_struct(*res);
    int c = a1_mat->RowIndex[i + 1] - a1_mat->RowIndex[i];
    if((cnt > 0) || (c>0))
    {
//      printf("%d %d\n", cnt, c);
      a2_mat->Value[k] = trace(*res);
      a2_mat->Col[k] = i;
      k++;
    }
    a2_mat->RowIndex[i + 1] = k;
    delete res;
  }
  a2_mat->NZ = k;
//  printMatrixVal(a2_mat);
//  crsMatrix * a_res = new crsMatrix(N_mat, k * k);
  crsMatrix * a1_i_a2_mat = new crsMatrix(N_mat, a2_mat->NZ);
  //a1_i_a2_mat -> RowIndex[0] = 0;
  //a1_i_a2_mat -> RowIndex[1] = a2_mat->NZ;
  k = 0;
  int k1, k2;
  for(int i = 0; i < N_mat; i++)
  {
    a1_i_a2_mat->RowIndex[i] = k;
    int c1 = a1_mat->RowIndex[i + 1] - a1_mat->RowIndex[i];
    int c2 = a2_mat->RowIndex[i + 1] - a2_mat->RowIndex[i];
//    a_res->RowIndex[i] = k * a2_mat->NZ;
    if(c1 > 0)
    {
      a1_i_a2_mat->Value[k] = a1_mat->Value[a1_mat->RowIndex[i]];
    }
    if(c2 > 0)
    {
      a1_i_a2_mat->Value[k].re += a2_mat->Value[a2_mat->RowIndex[i]].im;
      a1_i_a2_mat->Value[k].im -= a2_mat->Value[a2_mat->RowIndex[i]].re;
    }
    if(c1 + c2 > 0)
    {
      a1_i_a2_mat->Col[k] = i;
      k++;
    }
  }
  a1_i_a2_mat->RowIndex[N_mat] = k;
  a1_i_a2_mat->NZ = k;
//  a_res->RowIndex[N_mat] = k * a2_mat->NZ;
//  
//  for(int i = 0; i < k; i++)
//  {
//    for(int j = 0; j < k; j++)
//    {
//      a_res->Value[i * k + j].re = 0.0;
//      a_res->Value[i * k + j].re += a1_i_a2_mat->Value[i].re * a1_i_a2_mat->Value[j].re;
//      a_res->Value[i * k + j].re -= a1_i_a2_mat->Value[i].im *(- a1_i_a2_mat->Value[j].im);
//      a_res->Value[i * k + j].im = 0.0;
//      a_res->Value[i * k + j].im += a1_i_a2_mat->Value[i].re *(- a1_i_a2_mat->Value[j].im);
//      a_res->Value[i * k + j].im += a1_i_a2_mat->Value[i].im * a1_i_a2_mat->Value[j].re;
//      a_res->Col[i * k + j] = a1_i_a2_mat->Col[j];
//    }
//  }

//  printMatrixVal(a1_i_a2_mat);
//  printMatrixVal(a_res);

//  delete a1_i_a2_mat;

  m->l_mat = a1_i_a2_mat;
//  m->a_mat = a_res;


//  vector<pair<int , dcomplex> > * a_std;
//  a_std = new vector<pair<int , dcomplex> >[N_mat];
//
//  //костыль. нужно переделать на матричные операции
//  for(int i = 0; i < N_mat; i++)
//  {
//    for(int j = 0; j < N_mat; j++)
//    {
//      dcomplex v1, v2, v3, v4, res1,res2, res;
//      int c1, c2, c3, c4;
//      c1 = a1_mat->RowIndex[i + 1] - a1_mat->RowIndex[i];
//      c2 = a1_mat->RowIndex[j + 1] - a1_mat->RowIndex[j];
//      c3 = a2_mat->RowIndex[i + 1] - a2_mat->RowIndex[i];
//      c4 = a2_mat->RowIndex[j + 1] - a2_mat->RowIndex[j];
//
//      if((c1 + c3) * (c2 + c4) > 0)
//      {
//
//        v1.re = a1_mat->Value[a1_mat->RowIndex[i]].re * c1;
//        v3.re = a2_mat->Value[a2_mat->RowIndex[i]].re * c3;
//        v2.re = a1_mat->Value[a1_mat->RowIndex[j]].re * c2;
//        v4.re = a2_mat->Value[a2_mat->RowIndex[j]].re * c4;
//
//        v1.im = a1_mat->Value[a1_mat->RowIndex[i]].im * c1;
//        v3.im = a2_mat->Value[a2_mat->RowIndex[i]].im * c3;
//        v2.im = a1_mat->Value[a1_mat->RowIndex[j]].im * c2;
//        v4.im = a2_mat->Value[a2_mat->RowIndex[j]].im * c4;
//        
//        res1.re = v1.re + v3.im;
//        res1.im = v1.im - v3.re;
//        
//        res2.re = v2.re + v4.im;
//        res2.im = v2.im - v4.re;
//
//        res2.im = -res2.im;
//        
//        res.re = res1.re * res2.re - res1.im * res2.im;
//        res.im = res1.re * res2.im + res1.im * res2.re;
//        a_std[i].push_back(make_pair(j, res));
//
//      }
//    }
//  }
//
//
  delete a1_mat;
  delete a2_mat;
   
  delete A1;
  delete A2;
//
//  m->a_mat = stdToCrs(a_std, N_mat);

//  delete [] a_std;
//  printMatrixVal(m->a_mat);
}