//
// Created by cheyulin on 8/8/16.
//

#ifndef CODES_YCHE_FINE_GRAINED_MERGE_SCHEDULER_H
#define CODES_YCHE_FINE_GRAINED_MERGE_SCHEDULER_H

#include <semaphore.h>
#include <pthread.h>

#include <memory>
#include <vector>
#include <iostream>

#include <boost/range.hpp>

#include "deprecated/parallel_utils/parallel_configuration.h"
#include "thread_pool_breakable.h"

using namespace std;
using namespace yche;

namespace yche {
    template<typename Container2D, typename Predicate, typename TrueCallback, typename FalseCallback>
    class BSPer {
    public:
        BSPer(int thread_count, vector<Container2D> &reduce_data_ptr_vector,
              Predicate pair_computation_func, TrueCallback success_action_func, FalseCallback fail_action_func)
                : thread_count_(thread_count), thread_pool_(thread_count),
                  reduce_data_ptr_vector_(reduce_data_ptr_vector), predicate_(pair_computation_func),
                  true_callback_(success_action_func), false_callback_(fail_action_func) {
        }

        Container2D Execute() {
            while (reduce_data_ptr_vector_.size() > 1) {
                right_reduce_data_index_ = reduce_data_ptr_vector_.size() - 1;
                ReduceComputation();
                reduce_data_ptr_vector_.erase(reduce_data_ptr_vector_.end() - 1);
            }
            return reduce_data_ptr_vector_[0];
        }

        virtual  ~BSPer() {
        }

    private:
        ThreadPoolBreakable thread_pool_;
        int thread_count_;

        vector<Container2D> reduce_data_ptr_vector_;
        unsigned long right_reduce_data_index_;

        Predicate predicate_;
        TrueCallback true_callback_;
        FalseCallback false_callback_;

        void ReduceComputation();
    };

    //Parallel the inner for loop
    template<typename ReduceDataType, typename ComputationFuncType, typename ActionFuncType, typename FailActionFuncType>
    void BSPer<ReduceDataType, ComputationFuncType,
            ActionFuncType, FailActionFuncType>::ReduceComputation() {
        int round_num = 0;
        //Here make use of thread pool
        for (auto &left_ele:*reduce_data_ptr_vector_[right_reduce_data_index_]) {
            bool is_terminate_in_advance = false;
            for (auto &right_ele:*reduce_data_ptr_vector_[0]) {
                thread_pool_.AddTask([&right_ele, &left_ele, this]() {
                    if (this->predicate_(left_ele, right_ele)) {
                        return BreakWithCallBackRetType(true, []() {
                            this->true_callback_(left_ele, right_ele);
                        });
                    } else
                        return BreakWithCallBackRetType();
                });
            }
            thread_pool_.WaitForBreakOrTerminate(is_terminate_in_advance);
            if (!is_terminate_in_advance) {
                false_callback_(left_ele, reduce_data_ptr_vector_[0]);
            }
        }
    }
}

#endif //CODES_YCHE_FINE_GRAINED_MERGE_SCHEDULER_H
