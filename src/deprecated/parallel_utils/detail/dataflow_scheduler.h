//
// Created by cheyulin on 4/25/16.
//

#ifndef CODES_YCHE_PARALLELIZER_H
#define CODES_YCHE_PARALLELIZER_H

#include <thread>

#include "deprecated/parallel_utils/parallel_configuration.h"
#include "reduce_scheduler.h"

namespace yche {
    using namespace std;

    template<typename AlgorithmType>
    class DataFlower {
    private:
        using BasicDataType = typename AlgorithmType::BasicDataType;
        using MergeDataType = typename AlgorithmType::MergeDataType;
        using ReduceDataType = typename AlgorithmType::ReduceDataType;

        vector<thread> thread_list_;
        unique_ptr<vector<unique_ptr<BasicDataType>>> global_computation_task_vec_ptr_;
        vector<pair<size_t, size_t>> local_computation_range_index_vec_;
        vector<vector<unique_ptr<ReduceDataType>>> reduce_task_vectors_;

        int thread_count_;
        int idle_count_;
        pthread_t *thread_handles;
        pthread_mutex_t counter_mutex_lock_;

        vector<bool> is_rec_mail_empty_;
        bool is_end_of_local_computation;

        void RingCommTaskRequest(int thread_id);

        void InitTasks();

        void DoLeftMerging();

    public:
        unique_ptr<AlgorithmType> algorithm_ptr_;

        void ParallelExecute();

        DataFlower(int thread_count, unique_ptr<AlgorithmType> algorithm_ptr) : thread_count_(thread_count) {
            algorithm_ptr_ = std::move(algorithm_ptr);
            thread_handles = new pthread_t[thread_count_];

            pthread_mutex_init(&counter_mutex_lock_, NULL);

            is_rec_mail_empty_.resize(thread_count_, true);
            is_end_of_local_computation = false;
            local_computation_range_index_vec_.resize(thread_count);
            reduce_task_vectors_.resize(thread_count_);
            idle_count_ = 0;
        }

        virtual ~DataFlower() {
            pthread_mutex_destroy(&counter_mutex_lock_);
            delete[]thread_handles;
        }
    };

    template<typename AlgorithmType>
    void DataFlower<AlgorithmType>::ParallelExecute() {
        using namespace std::chrono;
        auto start = high_resolution_clock::now();
        InitTasks();
        auto end = high_resolution_clock::now();
        auto elapsed = duration_cast<milliseconds>(end - start).count();
        cout << "Task Init Cost :" << elapsed << '\n' << "Thread_count:" << thread_count_ << endl;

        start = high_resolution_clock::now();
        for (auto id = 0; id < thread_count_; id++) { thread_list_.emplace_back(RingCommTaskRequest, id); }
        for (auto &thread_ref: thread_list_) { thread_ref.join(); }
        end = high_resolution_clock::now();
        elapsed = duration_cast<milliseconds>(end - start).count();
        cout << "Parallel Execution Cost :" << elapsed << '\n' << "Thread_count:" << thread_count_ << endl;

        start = high_resolution_clock::now();
        DoLeftMerging();
        end = high_resolution_clock::now();
        elapsed = duration_cast<milliseconds>(end - start).count();
        cout << "Left Merging Cost: " << elapsed << endl;
    }

    template<typename AlgorithmType>
    void DataFlower<AlgorithmType>::RingCommTaskRequest(int thread_id) {
        auto dst_index = (thread_id + 1) % thread_count_;
        auto src_index = (thread_id - 1 + thread_count_) % thread_count_;
        auto &range_pair = local_computation_range_index_vec_[thread_id];
        auto &local_reduce_queue = reduce_task_vectors_[thread_id];

        while (!is_end_of_local_computation) {
            auto task_size = range_pair.second - range_pair.first + 1;
            if (task_size == 0) {
                if (idle_count_ == thread_count_ - 1) {
                    is_end_of_local_computation = true;
                    break;
                } else {
                    pthread_mutex_lock(&counter_mutex_lock_);
                    //Deal With All Finish and Enter into Idle State
                    if (idle_count_ == thread_count_ - 1)
                        is_end_of_local_computation = true;
                    else
                        idle_count_++;
                    pthread_mutex_unlock(&counter_mutex_lock_);

                    is_rec_mail_empty_[dst_index] = false;
                    while (!is_rec_mail_empty_[dst_index]) {
                        if (local_reduce_queue.size() > 1) {
                            local_reduce_queue.front() = algorithm_ptr_->ReduceComputation(local_reduce_queue.front(),
                                                                                           local_reduce_queue.back());
                            local_reduce_queue.erase(local_reduce_queue.end() - 1);
                        }
                        if (is_end_of_local_computation) {
                            cout << "Finish Local Computation" << thread_id << endl;
                            break;
                        }
                    }
                    pthread_mutex_lock(&counter_mutex_lock_);
                    idle_count_--;
                    pthread_mutex_unlock(&counter_mutex_lock_);
                }
            } else {
                if (task_size > 1) {
                    //Check Flag
                    auto &nei_range_pair = local_computation_range_index_vec_[src_index];
                    if (is_rec_mail_empty_[thread_id] == false) {
                        //update neighbor computation range pair
                        nei_range_pair.second = range_pair.second;
                        nei_range_pair.first = nei_range_pair.second - task_size / 2 + 1;
                        range_pair.second = nei_range_pair.first - 1;

                        is_rec_mail_empty_[thread_id] = true;
                    }
                }
                //Do Local Computation
                auto result = algorithm_ptr_->LocalComputation((*global_computation_task_vec_ptr_)[range_pair.first]);
                range_pair.first++;

                //Directly Put into Reduce Task Vector
                local_reduce_queue.push_back(algorithm_ptr_->WrapMergeDataToReduceData(result));
            }
        }
    }

    template<typename AlgorithmType>
    void DataFlower<AlgorithmType>::DoLeftMerging() {
        vector<ReduceDataType> reduce_data_ptr_vec;
        //Do Left Merging, Current Impl Do not care about the branch cost since it is only called once
        for (auto i = 0; i < thread_count_; i++) {
            auto &local_reduce_queue = reduce_task_vectors_[i];
            while (local_reduce_queue.size() > 0) {
                auto &reduce_data_ptr = std::move(local_reduce_queue.back());
                local_reduce_queue.erase(local_reduce_queue.end() - 1);
                reduce_data_ptr_vec.push_back(std::move(reduce_data_ptr));
            }
        }

        cout << "Before Reducer" << endl;
        cout << "Reduce Data Size:" << reduce_data_ptr_vec.size() << endl;
        Reducer<decltype(reduce_data_ptr_vec), ReduceDataType, decltype(algorithm_ptr_->CmpReduceData), decltype(algorithm_ptr_->ReduceComputation)> reducer(
                thread_count_, reduce_data_ptr_vec, algorithm_ptr_->CmpReduceData, algorithm_ptr_->ReduceComputation);
        algorithm_ptr_->overlap_community_vec_ = std::move(reducer.ParallelExecute());
    }
}


#endif //CODES_YCHE_PARALLELIZER_H
