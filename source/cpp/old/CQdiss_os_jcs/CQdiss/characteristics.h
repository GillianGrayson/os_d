#pragma once
#include "data.h"
#include "Model.h"

void characteristics_std(Model *m, RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd, int dump_id);
void characteristics_deep(Model *m, RunParam &rp, ConfigParam &cp, MainData &md, PropData &pd, int dump_id);
