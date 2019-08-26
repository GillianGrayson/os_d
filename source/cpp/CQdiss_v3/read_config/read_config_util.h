#ifndef __READ_CONFIG_UTIL__
#define __READ_CONFIG_UTIL__

#include "Config.h"

void set_param(ConfigParam &param,char *s, double value);
void read_file_conf(ConfigParam &param, FILE *f_conf);

#endif
