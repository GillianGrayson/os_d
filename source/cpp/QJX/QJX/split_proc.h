#pragma once
#include "config.h"
#include "data.h"

Split * init_split_structure_deep(RunParam * rp, ConfigParam * cp, MainData * md);

Split * init_split_structure(RunParam * rp, ConfigParam * cp, MainData * md);

void init_split_branches(Split * branch, int branch_id, RunParam * rp, ConfigParam * cp, MainData * md);

void copy_branch_not_member(Split * src, Split * dst);

void copy_struct_not_member(Split * src, Split * dst);

void delete_branch(Split * branch);

void delete_split_struct(Split * head);

void delete_branch_not_member(Split * branch);

void delete_split_struct_not_member(Split * head);