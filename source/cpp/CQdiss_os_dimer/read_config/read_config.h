#ifndef __READ_CONFIG__
#define __READ_CONFIG__

#include "Config.h"

void set_param(RunParam& rp, ConfigParam& cp, string str, string val);
void init_params(RunParam& rp, ConfigParam& cp, char* file_name);
void output_params(RunParam& rp, ConfigParam& param);

#endif
