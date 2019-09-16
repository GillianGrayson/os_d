#pragma once
#include "config.h"
#include "data.h"

void dump_adr_single(AllData * ad, int tr_id, bool append);

void dump_adr_avg(AllData * ad, bool append);

void dump_adr_avg_mult(AllData * ad, bool append, int begin_traj_id, int end_traj_id);

void dump_cd(AllData * ad);
