#pragma once
#include "config.h"
#include "data.h"

Split * init_split_structure(RunParam * rp, ConfigParam * cp, MainData * md);
void init_split_branches(Split * branch, int branch_id, RunParam * rp, ConfigParam * cp, MainData * md);