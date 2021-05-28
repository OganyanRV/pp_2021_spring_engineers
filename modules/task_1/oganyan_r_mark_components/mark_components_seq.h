//
// Created by ogrob on 5/27/21.
//

#ifndef MODULES_TASK_1_OGANYAN_R_MARK_COMPONENTS_MARK_COMPONENTS_SEQ_H_
#define MODULES_TASK_1_OGANYAN_R_MARK_COMPONENTS_MARK_COMPONENTS_SEQ_H_

#include <vector>
#include <cstdint>
#include <queue>
#include <stdexcept>
#include <utility>
#include "../../modules/task_2/oganyan_r_mark_components_omp/Disjoint_Set_Union.h"


std::pair<std::vector<int>, int> MarkComponentsSeq(std::vector<int> *img,
                                                   int height, int width);

#endif  //  MODULES_TASK_1_OGANYAN_R_MARK_COMPONENTS_MARK_COMPONENTS_SEQ_H_

