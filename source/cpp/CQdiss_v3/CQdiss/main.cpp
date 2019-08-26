#include <omp.h>

#include "read_config.h"
#include "Model.h"
#include "InitFs.h"
#include "initH.h"
#include "Init_a1_a2.h"
#include "Init_f_d.h"
#include "TranspFs.h"

#include "CalcQs.h"
#include "CalcKs.h"
#include "CalcRs.h"
#include "CalcGs.h"

#include "linSolv.h"
#include "calcODE.h"
#include "calcODE_real.h"

#include "CalcRho.h"

#include "CalcEigs.h"
#include "CalcTraseRO2.h"

#include "testODE.h"

double number_of_allocs = 0;
/*
unsigned int number_of_allocs_tmp = 0;

void* operator new(std::size_t size) throw(std::bad_alloc){
  number_of_allocs+=size;
  void *p = malloc(size);
  if (!p) throw std::bad_alloc();
  return p;
}

void* operator new  [](std::size_t size) throw(std::bad_alloc) {
  number_of_allocs += size;
  void *p = malloc(size);
  if (!p) throw std::bad_alloc();
  return p;
}

void* operator new  [](std::size_t size, const std::nothrow_t&) throw() {
  number_of_allocs += size;
  return malloc(size);
}
void* operator new   (std::size_t size, const std::nothrow_t&) throw(){
  number_of_allocs += size;
  return malloc(size);
}


void operator delete(void* ptr) throw() { free(ptr); }
void operator delete (void* ptr, const std::nothrow_t&) throw() { free(ptr); }
void operator delete[](void* ptr) throw() { free(ptr); }
void operator delete[](void* ptr, const std::nothrow_t&) throw() { free(ptr); }
//*/

int _main(int argc, char ** argv)
{
  printf("current path (%s)\n\n", argv[0]);
  if (argc >1)
    omp_set_num_threads(atoi(argv[1]));

  double time, all_time;
  ConfigParam param;
  read_config(param, "config.txt");

  //param.g = 1.0;
  //param.E0 = 0.0;
  //param.U = 5.0;
  //param.J = 1.0;
  //param.N = 100;
  //param.N = 3;

  FILE * memlog = fopen("memlog.txt", "w");

  Model * model;
  model = createModel(param.N, param);

  fprintf(memlog, "model = createModel(param.N, param); %zu\n", number_of_allocs);

  all_time = omp_get_wtime();
  time = omp_get_wtime();
  initFs(model->Fs, model->N);
  time = omp_get_wtime() - time;
  printf("initFs: %2.4lf\n", time);
  fflush(stdout);
  
  fprintf(memlog, "initFs(model->Fs, model->N); %zu\n", number_of_allocs);
  number_of_allocs = 0;

  time = omp_get_wtime();
 // init_h_vector(model);
  init_h_vector_opt(model);
  time = omp_get_wtime() - time;
  printf("init_h_vector: %2.4lf\n", time);
  fflush(stdout);

  fprintf(memlog, "init_h_vector(model); %zu\n", number_of_allocs);
  number_of_allocs = 0;

  time = omp_get_wtime();
  init_he_vector_opt(model);
  time = omp_get_wtime() - time;
  printf("init_he_vector: %2.4lf\n", time);
  fflush(stdout);
  
  fprintf(memlog, "init_he_vector(model); %zu\n", number_of_allocs);
  number_of_allocs = 0;

  time = omp_get_wtime();
//  init_a1_a2(model);
  init_a1_a2_opt(model);
  time = omp_get_wtime() - time;
  printf("init_a1_a2: %2.4lf\n", time);
  fflush(stdout);
  
  fprintf(memlog, "init_a1_a2(model); %zu\n", number_of_allocs);
  number_of_allocs = 0;

  //time = omp_get_wtime();
  //init_f_d(model);
  //time = omp_get_wtime() - time;
  //printf("init_f_d (6 for): %2.4lf\n", time);

  time = omp_get_wtime();
  init_f_d_valentin(model);
  time = omp_get_wtime() - time;
  printf("init_f_d (valentin): %2.4lf\n", time);
  fflush(stdout);
  fprintf(memlog, "init_f_d_valentin(model); %zu\n", number_of_allocs);
  number_of_allocs = 0;
  //time = omp_get_wtime();
  //transpFs(model);
  //time = omp_get_wtime() - time;
  //printf("transpF: %2.4lf\n", time);
  //fprintf(memlog, "transpFs(model); %d\n", number_of_allocs);

  time = omp_get_wtime();
  calcQs(model);
  time = omp_get_wtime() - time;
  printf("calcQs: %2.4lf\n", time);
  fflush(stdout);
//  printMatrixVal(model->Qs);
  fprintf(memlog, "calcQs(model); %zu\n", number_of_allocs);
  number_of_allocs = 0;

  time = omp_get_wtime();
  calcQEs(model);
  time = omp_get_wtime() - time;
  printf("calcQEs: %2.4lf\n", time);
  fflush(stdout);
//  printMatrixVal(model->QEs);
  fprintf(memlog, "calcQEs(model); %zu\n", number_of_allocs);
  number_of_allocs = 0;

  time = omp_get_wtime();
  calcKs(model);
  time = omp_get_wtime() - time;
  printf("calcKs: %2.4lf\n", time);
  fflush(stdout);
  fprintf(memlog, "calcKs(model); %zu\n", number_of_allocs);
  number_of_allocs = 0;

  time = omp_get_wtime();
  calcRs(model);
  time = omp_get_wtime() - time;
  printf("calcRs: %2.4lf\n", time);
  fflush(stdout);
  fprintf(memlog, "calcRs(model); %zu\n", number_of_allocs);
  number_of_allocs = 0;

  time = omp_get_wtime();
  calcGs(model);
  time = omp_get_wtime() - time; 
  printf("calcGs: %2.4lf\n", time);
  fflush(stdout);
  fprintf(memlog, "calcGs(model); %zu\n", number_of_allocs);
  number_of_allocs = 0;

//  printMatrixVal(model->Gs);
//  printVectorVal(model->Ks, model->N_mat);
//  printMatrixVal(model->Rs);
//  printMatrix(model->Rs);
//  printMatrixVal(model->Gs);

  if(model->conf.hasDriving == 0)
  {
    time = omp_get_wtime();
    linSolv(model);
    //linSolvReal(model);
    time = omp_get_wtime() - time;
    printf("linSolv: %2.4lf\n", time);
    saveMatrix("Gs.txt", model->Gs);
  }
  else
  {
    time = omp_get_wtime();
    initRhoODE(model);
//    printVectorVal(model->RhoF, model->N_mat);
    time = omp_get_wtime() - time;
    printf("initRhoODE: %2.4lf\n", time);
        
//    dcomplex diff;
//  
//    FILE * f = fopen("diff_itr.txt", "w");
  
    time = omp_get_wtime();
  
//    for(int itr = 0; itr < model->conf.N_T; itr++)
//    {
//      calcODE(model, model->conf.h, 
//              model->conf.NSTEP, itr * model->conf.T);
////      diff = calcDiffIter(model);
////      fprintf(f, "%e %e \n", diff.re, diff.im);
//    }
  
    complex_to_real(model->Gs->Value, model->Gs->NZ);
    complex_to_real(model->QEs->Value, model->QEs->NZ);
    complex_to_real(model->Ks, model->N_mat);
    complex_to_real(model->RhoF, model->N_mat);
    toOneBase(*(model->Gs));
    toOneBase(*(model->QEs));

    for(int itr = 0; itr < model->conf.N_T; itr++)
    {
      calcODE_real(model, model->conf.h, 
              model->conf.NSTEP, itr * model->conf.T);
    }

    toZeroBase(*(model->Gs));
    toZeroBase(*(model->QEs));
    real_to_complex(model->Gs->Value, model->Gs->NZ);
    real_to_complex(model->QEs->Value, model->QEs->NZ);
    real_to_complex(model->Ks, model->N_mat);
    real_to_complex(model->RhoF, model->N_mat);


    time = omp_get_wtime() - time;
  
//    fclose(f);
    printf("calcODE: %2.4lf\n", time);
	fflush(stdout);
  }
  fprintf(memlog, "calcODE_real; %zu\n", number_of_allocs);
  number_of_allocs = 0;

  linSolvCheck(model);
//  printVectorVal(model->RhoF, model->N_mat);

  time = omp_get_wtime();
  calcRho(model);
  time = omp_get_wtime() - time;
  printf("calcRho: %2.4lf\n", time);
  fprintf(memlog, "calcRho(model); %zu\n", number_of_allocs);
  number_of_allocs = 0;

  all_time = omp_get_wtime() - all_time;

  printf("all time: %2.4lf\n", all_time);
//  printMatrixVal(model->Rho);
  saveMatrix("Rho.txt", model->Rho);

  time = omp_get_wtime();
  calcTraseRO2(model);
  time = omp_get_wtime() - time;
  printf("calcTraseRO2: %2.4lf\n", time);
  fprintf(memlog, "calcTraseRO2(model); %zu\n", number_of_allocs);
  number_of_allocs = 0;
//  testODE(model);

  if(model->conf.CalcEig == 1)
  {
    time = omp_get_wtime();
    calcEig(model);
    time = omp_get_wtime() - time;
    printf("calcEig: %2.4lf\n", time);
  }

  saveVectorVal("RhoF.txt", model->RhoF, 1, model->N_mat);
  saveAbsMatrixVal("absRho.txt", model->Rho);
  saveAngleMatrixVal("angleRho.txt", model->Rho);

  freeModel(model);
  
  fclose(memlog);

  return 0;
}
