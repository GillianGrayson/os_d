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


double number_of_allocs = 0.0;

std::map<void *, double> mem_map;

FILE * mem_time;

bool ismap = false;
bool ismap_print = false;
double map_size = 0.0;

bool program = true;
bool isEnd = false;
void* operator new(std::size_t size) throw(std::bad_alloc){
  number_of_allocs += size / 1024.0 / 1024.0;
  void *p = malloc(size);

  if (mem_time != NULL)
  {
    if (ismap)
    {
      map_size += size / 1024.0 / 1024.0;
    }
    if (ismap_print)
    {
      fprintf(mem_time, "%lf \n", map_size);
      fprintf(mem_time, "%lf \n", -map_size);
      fflush(mem_time);
      ismap = false;
      ismap_print = false;
    }
    if (program && !ismap)
    {
      program = false;
      std::pair<void *, double> pr(p, size / 1024.0 / 1024.0);
      program = false;
      mem_map.insert(pr);

      fprintf(mem_time, "%lf \n", mem_map[p]);
      fflush(mem_time);
    }
  }
  program = true;
  if (!p) throw std::bad_alloc();
  return p;
}
void* operator new   (std::size_t size, const std::nothrow_t&) throw(){
  number_of_allocs += size / 1024.0 / 1024.0;
  void *p = malloc(size);

  if (mem_time != NULL)
  {
    if (ismap)
    {
      map_size += size / 1024.0 / 1024.0;
    }
    if (ismap_print)
    {
      fprintf(mem_time, "%lf \n", map_size);
      fprintf(mem_time, "%lf \n", -map_size);
      fflush(mem_time);
      ismap = false;
      ismap_print = false;
    }
    if (program && !ismap)
    {
      program = false;
      std::pair<void *, double> pr(p, size / 1024.0 / 1024.0);
      program = false;
      mem_map.insert(pr);

      fprintf(mem_time, "%lf \n", mem_map[p]);
      fflush(mem_time);
    }
  }
  program = true;

  return p;
}

void* operator new  [](std::size_t size) throw(std::bad_alloc) {
  number_of_allocs += size / 1024.0 / 1024.0;
  void *p = malloc(size);
  
  if (mem_time != NULL)
  {
    if (ismap)
    {
      map_size += size / 1024.0 / 1024.0;
    }
    if (ismap_print)
    {
      fprintf(mem_time, "%lf \n", map_size);
      fprintf(mem_time, "%lf \n", -map_size);
      fflush(mem_time);
      ismap = false;
      ismap_print = false;
    }
    if (program && !ismap)
    {
      program = false;
      std::pair<void *, double> pr(p, size / 1024.0 / 1024.0);
      program = false;
      mem_map.insert(pr);

      fprintf(mem_time, "%lf \n", mem_map[p]);
      fflush(mem_time);
    }
  }
  program = true;

  if (!p) throw std::bad_alloc();
  return p;
}

void* operator new  [](std::size_t size, const std::nothrow_t&) throw() {
  number_of_allocs += size / 1024.0 / 1024.0;
  void *p = malloc(size);
  
  if (mem_time != NULL)
  {
    if (ismap)
    {
      map_size += size / 1024.0 / 1024.0;
    }
    if (ismap_print)
    {
      fprintf(mem_time, "%lf \n", map_size);
      fprintf(mem_time, "%lf \n", -map_size);
      fflush(mem_time);
      ismap = false;
    }
    if (program && !ismap)
    {
      program = false;
      std::pair<void *, double> pr(p, size / 1024.0 / 1024.0);
      program = false;
      mem_map.insert(pr);

      fprintf(mem_time, "%lf \n", mem_map[p]);
      fflush(mem_time);
    }
  }
  program = true;

  return p;
}


void operator delete(void* ptr) throw() {  
  if (!isEnd) if (mem_map.count(ptr) > 0){
    fprintf(mem_time, "%lf \n", -mem_map[ptr]);
    fflush(mem_time);
  }
free(ptr);
}
void operator delete (void* ptr, const std::nothrow_t&) throw() {  
  if (!isEnd)if (mem_map.count(ptr) > 0) {
    fprintf(mem_time, "%lf \n", -mem_map[ptr]);
    fflush(mem_time);
  }
free(ptr);
}
void operator delete[](void* ptr) throw() { 
  if (!isEnd)if (mem_map.count(ptr) > 0) {
    fprintf(mem_time, "%lf \n", -mem_map[ptr]);
    fflush(mem_time);
  }
free(ptr);
}
void operator delete[](void* ptr, const std::nothrow_t&) throw() { 
  if (!isEnd)if (mem_map.count(ptr) > 0) {
    fprintf(mem_time, "%lf \n", -mem_map[ptr]);
    fflush(mem_time);
  }
free(ptr);
}
//*/

int main(int argc, char ** argv)
{
  mem_time = fopen("mem_time.txt", "w");
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

  //FILE * memlog = fopen("memlog.txt", "w");

  Model * model;
  fprintf(mem_time, "model = createModel(param.N, param); \n");

  model = createModel(param.N, param);
  FILE * memlog = model->memlog;


  fprintf(memlog, "model = createModel(param.N, param); %lf\n", number_of_allocs);
  fflush(memlog);
  
  all_time = omp_get_wtime();
//  time = omp_get_wtime();
//  initFs(model->Fs, model->N);
//  time = omp_get_wtime() - time;
//  printf("initFs: %2.4lf\n", time);
//  fflush(stdout);
//  fprintf(memlog, "initFs(model->Fs, model->N); %lf\n", number_of_allocs);

  number_of_allocs = 0;

  fprintf(mem_time, "init_h_vector(model); \n");
  time = omp_get_wtime();
  // init_h_vector(model);
  init_h_vector_opt(model);
  time = omp_get_wtime() - time;
  printf("init_h_vector: %2.4lf\n", time);
  fflush(stdout);

  fprintf(memlog, "init_h_vector(model); %lf\n", number_of_allocs);
  fflush(memlog);
  number_of_allocs = 0;
  
  fprintf(mem_time, "init_he_vector(model); \n");
  time = omp_get_wtime();
  init_he_vector_opt(model);
  time = omp_get_wtime() - time;
  printf("init_he_vector: %2.4lf\n", time);
  fflush(stdout);

  fprintf(memlog, "init_he_vector(model); %lf\n", number_of_allocs);
  fflush(memlog);
  number_of_allocs = 0;
  
  fprintf(mem_time, "init_a1_a2(model); \n");
  time = omp_get_wtime();
  //  init_a1_a2(model);
  init_a1_a2_opt(model);
  time = omp_get_wtime() - time;
  printf("init_a1_a2: %2.4lf\n", time);
  fflush(stdout);

  fprintf(memlog, "init_a1_a2(model); %lf\n", number_of_allocs);
  fflush(memlog);
  number_of_allocs = 0;
  
  //time = omp_get_wtime();
  //init_f_d(model);
  //time = omp_get_wtime() - time;
  //printf("init_f_d (6 for): %2.4lf\n", time);

  fprintf(mem_time, "init_f_d_valentin(model); \n");
  time = omp_get_wtime();
  init_f_d_valentin(model);
  time = omp_get_wtime() - time;
  printf("init_f_d (valentin): %2.4lf\n", time);
  fflush(stdout);
  fprintf(memlog, "init_f_d_valentin(model); %lf\n", number_of_allocs);
  fflush(memlog);
  number_of_allocs = 0;
  //time = omp_get_wtime();
  //transpFs(model);
  //time = omp_get_wtime() - time;
  //printf("transpF: %2.4lf\n", time);
  //fprintf(memlog, "transpFs(model); %d\n", number_of_allocs);

  fprintf(mem_time, "calcQs(model); \n");
  time = omp_get_wtime();
  calcQs(model);
  time = omp_get_wtime() - time;
  printf("calcQs: %2.4lf\n", time);
  fflush(stdout);
  //  printMatrixVal(model->Qs);
  fprintf(memlog, "calcQs(model); %lf\n", number_of_allocs);
  fflush(memlog);
  number_of_allocs = 0;
  
  
  fprintf(mem_time, "calcQEs(model); \n");
  time = omp_get_wtime();
  calcQEs(model);
  time = omp_get_wtime() - time;
  printf("calcQEs: %2.4lf\n", time);
  fflush(stdout);
  //  printMatrixVal(model->QEs);
  fprintf(memlog, "calcQEs(model); %lf\n", number_of_allocs);
  fflush(memlog);
  number_of_allocs = 0;
  
  
  fprintf(mem_time, "calcKs(model); \n");
  time = omp_get_wtime();
  calcKs(model);
  time = omp_get_wtime() - time;
  printf("calcKs: %2.4lf\n", time);
  fflush(stdout);
  fprintf(memlog, "calcKs(model); %lf\n", number_of_allocs);
  fflush(memlog);
  number_of_allocs = 0;

  fprintf(mem_time, "calcRs(model); \n");
  time = omp_get_wtime();
  calcRs(model);
  time = omp_get_wtime() - time;
  printf("calcRs: %2.4lf\n", time);
  fflush(stdout);
  fprintf(memlog, "calcRs(model); %lf\n", number_of_allocs);
  fflush(memlog);
  number_of_allocs = 0;

  fprintf(mem_time, "calcGs(model); \n");
  time = omp_get_wtime();
  calcGs(model);
  time = omp_get_wtime() - time;
  printf("calcGs: %2.4lf\n", time);
  fflush(stdout);
  fprintf(memlog, "calcGs(model); %lf\n", number_of_allocs);
  fflush(memlog);
  number_of_allocs = 0;

  //saveMatrix_coor("Gs_p.txt", model->Gs);

  //  printMatrixVal(model->Gs);
  //  printVectorVal(model->Ks, model->N_mat);
  //  printMatrixVal(model->Rs);
  //  printMatrix(model->Rs);
  //  printMatrixVal(model->Gs);

  fprintf(mem_time, "calcODE(model); \n");

  if (model->conf.hasDriving == 0)
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

    for (int itr = 0; itr < model->conf.N_T; itr++)
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
  fprintf(memlog, "calcODE_real; %lf\n", number_of_allocs);
  number_of_allocs = 0;

  linSolvCheck(model);
  //  printVectorVal(model->RhoF, model->N_mat);

  fprintf(mem_time, "calcRho_fill(model); \n");
  time = omp_get_wtime();
//  calcRho(model);
//  printMatrixVal(model->Rho);
  calcRho_fill(model);
  //printMatrixVal(model->Rho);
  time = omp_get_wtime() - time;
  printf("calcRho: %2.4lf\n", time);
  fprintf(memlog, "calcRho(model); %lf\n", number_of_allocs);
  fflush(memlog);
  number_of_allocs = 0;

  all_time = omp_get_wtime() - all_time;

  printf("all time: %2.4lf\n", all_time);
  dcomplex tr = trace(*(model->Rho));
  printf("Trace: %lf\n", tr.re);
  saveMatrix("Rho.txt", model->Rho);

  fprintf(mem_time, "calcTraseRO2(model); \n");
  time = omp_get_wtime();
  calcTraseRO2(model);
  time = omp_get_wtime() - time;
  printf("calcTraseRO2: %2.4lf\n", time);
  fprintf(memlog, "calcTraseRO2(model); %lf\n", number_of_allocs);
  fflush(memlog);
  number_of_allocs = 0;
  //  testODE(model);

  
  fprintf(mem_time, "other \n");
  if (model->conf.CalcEig == 1)
  {
    time = omp_get_wtime();
    calcEig(model);
    time = omp_get_wtime() - time;
    printf("calcEig: %2.4lf\n", time);
  }


  saveVectorVal("RhoF.txt", model->RhoF, 1, model->N_mat);
  saveMatrixVal("rho_f.txt", model->Rho);
  saveAbsMatrixVal("absRho.txt", model->Rho);
  saveAngleMatrixVal("angleRho.txt", model->Rho);

  freeModel(model);
  isEnd = true;
//  fclose(memlog);
  fclose(mem_time);
  return 0;
}
