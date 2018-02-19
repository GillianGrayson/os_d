#pragma once
#include "config.h"
#include "data.h"
#include "data_qj.h"

void dump_adr_single(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int tr_id, bool append);

void dump_adr_avg(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, bool append);

void update_evo_std(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int dump_id);

void update_evo_lpn(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd, int dump_id);

void dump_std(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd);

void dump_lpn(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd);

void dump_cd(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd);

void dump_evo_std(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd);

void dump_evo_lpn(RunParam * rp, ConfigParam * cp, MainData * md, QJData * qjd);