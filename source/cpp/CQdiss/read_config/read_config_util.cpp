#include "read_config.h"
#include <stdio.h>
#include <string.h>

void set_param(ConfigParam &param,char *s, double value)
{
  if(strcmp(s, "N,")== 0)
  {
    param.N = (int)value;
  }
  if(strcmp(s, "N_T,")== 0)
  {
    param.N_T = (int)value;
  }
  if(strcmp(s, "NSTEP,")== 0)
  {
    param.NSTEP = (int)value;
    param.h  = param.T / param.NSTEP;    
  }
  
  if(strcmp(s, "E0,")== 0)
  {
    param.E0 = value;
  }

  if(strcmp(s, "J,")== 0)
  {
    param.J = value;
  }
  
  if(strcmp(s, "U,")== 0)
  {
    param.U = value;
  }
  
  if(strcmp(s, "g,")== 0)
  {
    param.g = value / param.N;
  }
  if(strcmp(s, "A0,")== 0)
  {
    param.A0 = value;
  }
  if(strcmp(s, "w,")== 0)
  {
    param.w = value;
    param.T  = 2.0 * PI / param.w;
    param.h  = param.T / param.NSTEP;
    
  }
  if(strcmp(s, "hasDriving,") == 0)
  {
    param.hasDriving = (int) value;
  }
  if (strcmp(s, "driving_type,") == 0)
  {
	  param.driving_type = (int)value;
  }
  if(strcmp(s, "CalcEig,") == 0)
  {
    param.CalcEig = (int) value;
  }
}

void read_file_conf(ConfigParam &param, FILE *f_conf)
{
  char s[80];
  double val=1;
  int cnt;
  while(!feof(f_conf))
  {
    cnt = fscanf(f_conf, "%s", s);
    if(cnt > 0)
    {
      fscanf(f_conf, "%lf", &val);
      //printf("*%s*%lf*\n", s, val);
      set_param(param, s, val);
    }
  }
  
  fclose(f_conf);
}