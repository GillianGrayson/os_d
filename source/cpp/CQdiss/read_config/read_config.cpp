#include <stdio.h>
#include "read_config.h"
#include "read_config_util.h"


void output_setting(ConfigParam &param)
{
  printf("############# parameter task #############\n");
  printf("N, %d\n", param.N);

  printf("J , %lf \n",param.J );
  printf("E0, %lf \n",param.E0);
  printf("U , %lf \n",param.U );
  printf("g , %lf \n",param.g );
  printf("CalcEig, %d\n", param.CalcEig);

  printf("hasDriving, %d \n", param.hasDriving);
  if(param.hasDriving == 1)
  {
     printf("T, %lf\n", param.T);
     printf("A0, %lf\n", param.A0);
     printf("N T, %d\n", param.N_T);
     printf("h, %lf\n", param.h);
  }
  
  printf("############# parameter task #############\n");
}

void read_config(ConfigParam &param, char * file_name)
{
  FILE * f_conf = fopen(file_name, "r");
  if(f_conf == NULL)
  {
    printf("file %s not found \n", file_name);
    printf("use default setting \n");
  } 
  else
  {
    read_file_conf(param, f_conf);
  }
  output_setting(param);
}
