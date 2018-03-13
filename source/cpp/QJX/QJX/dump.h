#pragma once
#include "config.h"
#include "data.h"

void dump_adr_single(AllData * ad, int tr_id, bool append);

void dump_adr_avg(AllData * ad, bool append);

void dump_std(AllData * ad);

void dump_lpn(AllData * ad);

void dump_cd(AllData * ad);

void dump_evo_std(AllData * ad);

void dump_evo_lpn(AllData * ad);