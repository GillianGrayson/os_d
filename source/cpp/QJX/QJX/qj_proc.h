#pragma once
#include "config.h"
#include "data.h"

void init_splits_deep(AllData * ad);
void init_splits(AllData * ad);

void free_splits_deep(AllData * ad);
void free_splits(AllData * ad);

Split * init_split_structure_cd(AllData * ad);

Split * init_split_structure(AllData * ad);

void init_split_branches(Split * branch, int branch_id, AllData * ad);

void copy_branch_not_member(Split * src, Split * dst);

void copy_struct_not_member(Split * src, Split * dst);

void delete_branch(Split * branch);

void delete_split_struct(Split * head);

void delete_branch_not_member(Split * branch);

void delete_split_struct_not_member(Split * head);