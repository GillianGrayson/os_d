#pragma once
#include "config.h"
#include "data.h"

void dump_adr_single(AllData * ad, int tr_id, bool append);

void dump_adr_avg(AllData * ad, bool append);

void update_evo_std(AllData * ad, int dump_id);

void update_evo_lpn(AllData * ad, int dump_id);

void dump_std(AllData * ad);

void dump_lpn(AllData * ad);

void dump_cd(AllData * ad);

void dump_evo_std(AllData * ad);

void dump_evo_lpn(AllData * ad);