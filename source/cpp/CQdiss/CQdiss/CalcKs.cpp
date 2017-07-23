#include "CalcKs.h"
#include "coef_coord.h"

void calcKs(Model *m)
{
  int N = m->N;
  int N_mat = m->N_mat;
  crsMatrix *Ks_tmp;
  crsMatrix *FsT;
  dcomplex  *Ks = m->Ks;
//  crsMatrix *As = m->a_mat;
//  crsMatrix **Fs = m->f_mat;
  crsMatrix *AsT;

  crsMatrix * l_mat = m->l_mat;


  ulli cnt = fijk_coord_sym(m->l_mat, m->N + 1);
  printf("cntSubFijk %u\n", cnt);
  m->f_ijk = new Tensor_Coordinates;
  fijk_coord_ch(m->f_ijk, m->l_mat, cnt+1, m->N + 1);
  Tensor_Coordinates * f_ijk = m->f_ijk;

  fprintf(m->memlog, "fijk_coord_ch(m->f_ijk, m->l_mat, cnt+1, m->N + 1); %lf\n", number_of_allocs);
  fflush(m->memlog);

  for (int i = 0; i < N_mat; i++)
  {
    Ks[i].re = 0.0;
    Ks[i].im = 0.0;
  }

  for (ulli i = 0; i < f_ijk->k; i++)
  {
    int ii = f_ijk->coord1[i];
    int ji = f_ijk->coord2[i];
    int ki = f_ijk->coord3[i];
    int c1 = l_mat->RowIndex[ii + 1] - l_mat->RowIndex[ii];
    int c2 = l_mat->RowIndex[ji + 1] - l_mat->RowIndex[ji];
    if (c1 * c2 > 0)
    {
      dcomplex v1 = l_mat->Value[l_mat->RowIndex[ii]];
      dcomplex v2 = l_mat->Value[l_mat->RowIndex[ji]];
      v2.im = -v2.im;
      dcomplex v3 = f_ijk->data[i];
      dcomplex tmp;
      tmp.re = v1.re * v2.re - v1.im * v2.im;
      tmp.im = v1.re * v2.im + v1.im * v2.re;

      Ks[ki].re += tmp.re * v3.re - tmp.im * v3.im;
      Ks[ki].im += tmp.re * v3.im + tmp.im * v3.re;
    }
  }

/////  printVectorVal(Ks, N_mat);
/////
/////  for(int i = 0; i < N_mat; i++)
/////  {
/////    Ks[i].re = 0.0;
/////    Ks[i].im = 0.0;
/////  }
///////  printMatrixVal(As);
/////  for(int i = 0; i < N_mat; i++)
/////  {
/////    AsT = As;
/////    //AsT = new crsMatrix(*(As));
/////    //Transpose(*(As), *AsT);
/////    FsT = Fs[i];
/////    FsT = new crsMatrix(*(Fs[i]));
/////    Transpose(*(Fs[i]), *FsT, false);
///////    printMatrixVal(FsT);
/////    for(int j = 0; j < N_mat; j++)
/////    {
/////      int ii,jj, m1, m2;
/////      ii = As->RowIndex[i];
/////      m1 = As->RowIndex[i + 1];
/////      jj = FsT->RowIndex[j];
/////      m2 = FsT->RowIndex[j + 1];
/////
/////      while((ii < m1)&&(jj < m2))
/////      {
/////        if(AsT->Col[ii] < FsT->Col[jj])
/////        {  
/////          ii++;
/////        } 
/////        else
/////        {
/////          if(AsT->Col[ii] > FsT->Col[jj])
/////          {
/////            jj++;
/////          } 
/////          else
/////          {
/////            dcomplex as, fs;
/////            as = AsT->Value[ii];
/////            fs = FsT->Value[jj];
/////            Ks[j].re += as.re * fs.re - as.im * fs.im;
/////            Ks[j].im += as.re * fs.im + as.im * fs.re;
/////
/////            ii++;
/////            jj++;
/////          }
/////        }
/////      }
/////
/////    }
/////    delete FsT;
/////    //delete AsT;
/////  }
/////  
/////  printVectorVal(Ks, N_mat);

  dcomplex val;
  for(int i = 0; i < N_mat; i++)
  {
    Ks[i].re *= (-m->conf.g) / (N + 1);
    Ks[i].im *= ( m->conf.g) / (N + 1);
    val = Ks[i];
    Ks[i].re = val.im;
    Ks[i].im = val.re;
  }
}