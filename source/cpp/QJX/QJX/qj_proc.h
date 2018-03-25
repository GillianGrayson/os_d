#pragma once
#include "config.h"
#include "data.h"

void init_split_branches(Split * branch, int branch_id, AllData * ad);

void copy_branch_not_member(Split * src, Split * dst);

void copy_struct_not_member(Split * src, Split * dst);

void delete_branch(Split * branch);

void delete_split_struct(Split * head);

void delete_branch_not_member(Split * branch);

void delete_split_struct_not_member(Split * head);

void prop_step(MKL_Complex16 * phi, MKL_Complex16 * matrix, MKL_Complex16 * res, int sys_size);

void one_period_branch(AllData * ad, Split * head, int tr_id, Split * branch);

void one_sub_period_deep(AllData * ad, int tr_id, int part_id, int thread_id);