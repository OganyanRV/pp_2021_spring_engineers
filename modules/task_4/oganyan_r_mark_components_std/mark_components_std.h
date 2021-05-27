// Copyright 2021 Oganyan Robert


#ifndef MODULES_TASK_4_OGANYAN_R_MARK_COMPONENTS_STD_MARK_COMPONENTS_STD_H_
#define MODULES_TASK_4_OGANYAN_R_MARK_COMPONENTS_STD_MARK_COMPONENTS_STD_H_


#include <utility>
#include <cstdint>
#include <stdexcept>
#include <vector>
#include <queue>
#include "../../modules/task_2/oganyan_r_mark_components_omp/Disjoint_Set_Union.h"




std::pair<std::vector<int>, int> MarkComponentsSeq_3(std::vector<int> *img,
                                                   int height, int width);

std::pair<std::vector<int>, int> MarkComponentsParStd(std::vector<int> *img,
                                                   int height, int width, int num_proc);

void bfs_3(std::vector<int> *img, std::pair<int, int> start,
         int *number, int width, int height);


#endif  //  MODULES_TASK_4_OGANYAN_R_MARK_COMPONENTS_STD_MARK_COMPONENTS_STD_H_
