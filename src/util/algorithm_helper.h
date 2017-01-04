//
// Created by cheyulin on 12/24/16.
//

#ifndef CODES_YCHE_ALGORITHM_HELPER_H
#define CODES_YCHE_ALGORITHM_HELPER_H

#include <map>
#include "parallel_utils/detail/dataflow_scheduler.h"
#include "parallel_utils/detail/reduce_scheduler.h"

namespace yche {
    template<typename Algorithm>
    void ExecuteAlgorithm(int thread_num, unique_ptr<Algorithm> &algorithm_ptr, map<int, int> &index_name_map) {
        cout << "Reduce Enabled" << endl;
        DataFlower<Algorithm> parallelizer(thread_num, std::move(algorithm_ptr));
        parallelizer.ParallelExecute();
        algorithm_ptr = std::move(parallelizer.algorithm_ptr_);

#ifndef NOT_COUT_COMMUNITY_RESULT
        auto communities_ptr_vec = std::move(algorithm_ptr->overlap_community_vec_);
        cout << "comm_size:" << communities_ptr_vec->size() << endl;
        for (auto &&community_ptr:*communities_ptr_vec) {
            for (auto member_id:*community_ptr) {
                cout << index_name_map[member_id] << "(" << member_id << "),";
            }
            cout << endl;
        }
#endif
    }
}

#endif //CODES_YCHE_ALGORITHM_HELPER_H
