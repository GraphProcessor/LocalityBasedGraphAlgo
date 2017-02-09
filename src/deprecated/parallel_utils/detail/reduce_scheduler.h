//
// Created by cheyulin on 5/7/16.
//

#ifndef CODES_YCHE_REDUCER_H
#define CODES_YCHE_REDUCER_H

#include "semaphore.h"

#include <thread>

#include <chrono>
#include <iostream>
#include <memory>
#include <vector>
#include <algorithm>

#include "deprecated/parallel_utils/parallel_configuration.h"

namespace yche {
    using namespace std;

    struct TaskIndices {
        size_t res_idx_;
        size_t beg_idx_;
        size_t end_idx_;

        TaskIndices(size_t beg_idx, size_t end_idx) : res_idx_(beg_idx), beg_idx_(beg_idx + 1), end_idx_(end_idx) {}

        void Reschedule(size_t begin_idx, size_t end_idx) {
            beg_idx_ = begin_idx;
            end_idx_ = end_idx;
        }

        size_t GetSize() {
            return end_idx_ - beg_idx_ + 1;
        }
    };

    template<typename Container, typename Data, typename Predicate, typename Combiner>
    class Reducer {
    private:
        size_t thread_count_;
        size_t data_count_;
        size_t idle_count_;
        pthread_t *thread_handles_;

        vector<Data> first_phase_reduce_data_pool_vec_;
        vector<Data> second_phase_global_reduce_data_vector;
        vector<TaskIndices<size_t>> reduce_data_indices_vec_;
        vector<bool> is_rec_mail_empty_;

        vector<sem_t> sem_mail_boxes_;
        vector<pthread_mutex_t> check_indices_mutex_lock_vector_;
        pthread_mutex_t counter_mutex_lock_;
        pthread_mutex_t task_taking_mutex_;
        pthread_cond_t task_taking_cond_var_;
        pthread_barrier_t timestamp_barrier_;

        vector<bool> is_busy_working_;
        bool is_reduce_task_only_one_;
        bool is_end_of_loop_;
        bool is_end_of_reduce_;

        Predicate data_cmp_function_;
        Combiner reduce_compute_function_;

        void RingCommTaskStealAndRequestThreadFunction(size_t thread_id);

        void InitDataPerThread(Container &data_collection);

    public:

        unique_ptr<Data> ParallelExecute();

        Reducer(size_t thread_count, Container &reduce_data_collection,
                Predicate data_cmp_function,
                Combiner compute_function)
                : thread_count_(thread_count), data_cmp_function_(data_cmp_function),
                  reduce_compute_function_(compute_function) {
            thread_handles_ = new pthread_t[thread_count_];
            data_count_ = 0;

            sem_mail_boxes_.resize(thread_count_);
            for (auto i = 0; i < thread_count_; i++) { sem_init(&sem_mail_boxes_[i], 0, 0); }
            check_indices_mutex_lock_vector_.resize(thread_count_);
            for (auto i = 0; i < thread_count_; i++) { pthread_mutex_init(&check_indices_mutex_lock_vector_[i], NULL); }

            pthread_barrier_init(&timestamp_barrier_, NULL, thread_count_);
            pthread_mutex_init(&task_taking_mutex_, NULL);
            pthread_mutex_init(&counter_mutex_lock_, NULL);
            pthread_cond_init(&task_taking_cond_var_, NULL);


            is_end_of_loop_ = false;
            is_end_of_reduce_ = false;
            is_reduce_task_only_one_ = false;
            is_rec_mail_empty_.resize(thread_count_, true);
            is_busy_working_.resize(thread_count_, false);
            idle_count_ = 0;
            //InitLocalData
            InitDataPerThread(reduce_data_collection);
        }

        virtual ~Reducer() {
            for (auto i = 0; i < thread_count_; i++) {
                sem_destroy(&sem_mail_boxes_[i]);
            }
            for (auto i = 0; i < thread_count_; i++) {
                pthread_mutex_destroy(&check_indices_mutex_lock_vector_[i]);
            }
            pthread_mutex_destroy(&task_taking_mutex_);
            pthread_mutex_destroy(&counter_mutex_lock_);
            pthread_cond_destroy(&task_taking_cond_var_);
            pthread_barrier_destroy(&timestamp_barrier_);
            delete[]thread_handles_;
        }
    };

    template<typename DataCollectionType, typename DataType, typename DataCmpFunctionType, typename ComputationFunctionType>
    void Reducer<DataCollectionType, DataType, DataCmpFunctionType, ComputationFunctionType>::InitDataPerThread(
            DataCollectionType &data_collection) {
        data_count_ = data_collection.size();
        if (data_count_ == 1) {
            is_reduce_task_only_one_ = true;
            first_phase_reduce_data_pool_vec_.push_back(std::move(*data_collection.begin()));
            return;
        }
        first_phase_reduce_data_pool_vec_.resize(data_collection.size());
        int idx = 0;
        for (auto &data:data_collection) {
            first_phase_reduce_data_pool_vec_[idx] = std::move(data);
            idx++;
        }
        auto task_per_thread = data_collection.size() / thread_count_;
        cout << task_per_thread << " tasks per thread" << endl;

        reduce_data_indices_vec_.reserve(thread_count_);
        for (auto i = 0; i < thread_count_ - 1; ++i) {
            reduce_data_indices_vec_.emplace(i * task_per_thread, (i + 1) * task_per_thread - 1);
        }
        reduce_data_indices_vec_[thread_count_ - 1].InitTaskIndices((thread_count_ - 1) * task_per_thread,
                                                                    data_collection.size() - 1);
        //Sort From Greatest to Least
        for (auto i = 0; i < thread_count_; i++) {
            sort(first_phase_reduce_data_pool_vec_.begin() + reduce_data_indices_vec_[i].beg_idx_,
                 first_phase_reduce_data_pool_vec_.begin() + reduce_data_indices_vec_[i].end_idx_,
                 data_cmp_function_);
        }
        cout << "Reduce Task Init Finished" << endl;
#ifdef DEBUG
        for (auto i = 0; i < thread_count_; ++i) {
            cout << reduce_data_indices_vec_[i].res_idx_ << ","
                 << reduce_data_indices_vec_[i].beg_idx_ << ","
                 << reduce_data_indices_vec_[i].end_idx_ << endl;
        }
#endif
    }

    template<typename DataCollectionType, typename DataType, typename DataCmpFunctionType, typename ComputationFunctionType>
    void
    Reducer<DataCollectionType, DataType, DataCmpFunctionType, ComputationFunctionType>::RingCommTaskStealAndRequestThreadFunction(
            size_t thread_id) {
        size_t thread_idx = thread_id;
        auto dst_idx = (thread_idx + 1) % thread_count_;
        auto src_idx = (thread_idx - 1 + thread_count_) % thread_count_;
        auto &local_reduce_data_indices = reduce_data_indices_vec_[thread_idx];

        using namespace std::chrono;
        time_point<high_resolution_clock> start_time_point, end_time_point;
        if (thread_idx == 0) {
            start_time_point = high_resolution_clock::now();
        }
        if (thread_count_ < data_count_) {
            while (!is_end_of_loop_) {
                auto reduce_data_size = local_reduce_data_indices.GetReduceTaskSize();
                if (reduce_data_size == 0) {
                    if (idle_count_ == thread_count_ - 1) {
                        is_end_of_loop_ = true;
                        idle_count_ = 0;
                        for (auto i = 0; i < thread_count_; i++)
                            if (is_rec_mail_empty_[i] == false)
                                sem_post(&sem_mail_boxes_[i]);
                    } else {
                        pthread_mutex_lock(&counter_mutex_lock_);
                        if (idle_count_ == thread_count_ - 1) {
                            is_end_of_loop_ = true;
                            idle_count_ = 0;
                            for (auto i = 0; i < thread_count_; i++)
                                if (is_rec_mail_empty_[i] == false)
                                    sem_post(&sem_mail_boxes_[i]);
                        } else
                            idle_count_++;
                        pthread_mutex_unlock(&counter_mutex_lock_);

                        bool is_going_to_request = false;

#ifdef STEAL_ENABLE
                        pthread_mutex_lock(&check_indices_mutex_lock_vector_[dst_idx]);
                        auto available_task_num = reduce_data_indices_vec_[dst_idx].GetReduceTaskSize();
                        if (available_task_num == 0 ||
                            (available_task_num == 1 && is_busy_working_[dst_idx] == true)) {
                            is_going_to_request = true;
                        } else {
                            if (is_busy_working_[dst_idx]) {
                                //Index Should Consider Dst Current Computation
                                auto current_end_idx = reduce_data_indices_vec_[dst_idx].end_idx_-1;
                                auto current_start_idx = current_end_idx - (reduce_data_size + 1) / 2 + 1;
                                pthread_mutex_lock(&check_indices_mutex_lock_vector_[thread_idx]);
                                reduce_data_indices_vec_[thread_idx].RescheduleTaskIndices(current_start_idx,
                                                                                             current_end_idx);
                                pthread_mutex_unlock(&check_indices_mutex_lock_vector_[thread_idx]);
                                reduce_data_indices_vec_[dst_idx].end_idx_=current_start_idx;
                            }
                            else {
                                //Not Required to Consider Dst Current Computation, since it not enter critical section
                                auto current_end_idx = reduce_data_indices_vec_[dst_idx].end_idx_;
                                auto current_start_idx = current_end_idx - (reduce_data_size + 1) / 2 + 1;
                                pthread_mutex_lock(&check_indices_mutex_lock_vector_[thread_idx]);
                                reduce_data_indices_vec_[thread_idx].RescheduleTaskIndices(current_start_idx,
                                                                                             current_end_idx);
                                pthread_mutex_unlock(&check_indices_mutex_lock_vector_[thread_idx]);
                                reduce_data_indices_vec_[dst_idx].end_idx_=current_start_idx-1;
                            }
                            //Update Current Thread Index and Dst Index
                        }
                        pthread_mutex_unlock(&check_indices_mutex_lock_vector_[dst_idx]);

                        if (is_going_to_request) {
                            is_rec_mail_empty_[dst_idx] = false;
                            sem_wait(&sem_mail_boxes_[dst_idx]);
                        }
#else
                        is_rec_mail_empty_[dst_idx] = false;
                        sem_wait(&sem_mail_boxes_[dst_idx]);
#endif
                        if (is_end_of_loop_) {
                            break;
                        }
                        pthread_mutex_lock(&counter_mutex_lock_);
                        idle_count_--;
                        pthread_mutex_unlock(&counter_mutex_lock_);
                    }
                } else {
                    if (reduce_data_size > 1) {
                        //Check Flag and Assign Tasks To Left Neighbor
                        if (is_rec_mail_empty_[thread_idx] == false) {
                            auto neighbor_end_idx = reduce_data_indices_vec_[thread_idx].end_computation_idx_;
                            auto neighbor_start_idx = neighbor_end_idx - reduce_data_size / 2 + 1;
                            reduce_data_indices_vec_[src_idx].RescheduleTaskIndices(
                                    neighbor_start_idx, neighbor_end_idx);
                            local_reduce_data_indices.end_computation_idx_ = neighbor_start_idx - 1;
                            is_rec_mail_empty_[thread_idx] = true;
                            sem_post(&sem_mail_boxes_[thread_idx]);
                        }
                    }

                    //To Avoid Enter into Computation when the neighbor is stealing
#ifdef STEAL_ENABLE
                    pthread_mutex_lock(&check_indices_mutex_lock_vector_[thread_idx]);
                    is_busy_working_[thread_idx] = true;
                    //Task Has been Steal
                    if (local_reduce_data_indices.GetReduceTaskSize() == 0)
                    {
                        is_busy_working_[thread_idx]= false;
                        continue;
                    }
                    pthread_mutex_unlock(&check_indices_mutex_lock_vector_[thread_idx]);
#endif
                    //Do reduce computation, use the first max one and the last min one
                    first_phase_reduce_data_pool_vec_[local_reduce_data_indices.result_idx_] = std::move(
                            reduce_compute_function_(
                                    first_phase_reduce_data_pool_vec_[local_reduce_data_indices.result_idx_],
                                    first_phase_reduce_data_pool_vec_[local_reduce_data_indices.end_computation_idx_]));
#ifdef STEAL_ENABLE
                    pthread_mutex_lock(&check_indices_mutex_lock_vector_[thread_idx]);
                    is_busy_working_[thread_idx] = false;
                    local_reduce_data_indices.end_idx_--;
                    pthread_mutex_unlock(&check_indices_mutex_lock_vector_[thread_idx]);
#else
                    local_reduce_data_indices.end_computation_idx_--;
#endif
                }
            }
        }

        //Barrier
        pthread_barrier_wait(&timestamp_barrier_);

        if (thread_idx == 0) {
            end_time_point = high_resolution_clock::now();
            auto duration_time = duration_cast<milliseconds>(end_time_point - start_time_point).count();
            cout << "Before Global Variable Task Acq In First-Phase Reduce:" << duration_time << endl;
        }

#ifndef REDUCE_2ND_PHASE_SEQUENTIAL
        //Reduce Data Size Has become much larger in this phase, Maybe Need Fine-Grained Parallelism
        //Do left things, 1) send data back to global variable 2) use condition variable to synchronize
        unique_ptr<DataType> result_data_ptr = std::move(
                first_phase_reduce_data_pool_vec_[reduce_data_indices_vec_[thread_idx].result_idx_]);
        unique_ptr<DataType> input_data_ptr;
        pthread_barrier_wait(&timestamp_barrier_);
        cout << "Thread Index:" << thread_idx << endl;

        while (!is_end_of_reduce_) {
            pthread_mutex_lock(&task_taking_mutex_);
            //Send result to global vector
            second_phase_global_reduce_data_vector.push_back(std::move(result_data_ptr));
            //Fetch Data : Busy Worker
            if (second_phase_global_reduce_data_vector.size() >= 2) {
                result_data_ptr = std::move(second_phase_global_reduce_data_vector.back());
                second_phase_global_reduce_data_vector.erase(second_phase_global_reduce_data_vector.end() - 1);
                input_data_ptr = std::move(second_phase_global_reduce_data_vector.back());
                second_phase_global_reduce_data_vector.erase(second_phase_global_reduce_data_vector.end() - 1);
                pthread_mutex_unlock(&task_taking_mutex_);

                //Do the computation After release the lock
                if (is_end_of_reduce_)
                    break;
                result_data_ptr = std::move(reduce_compute_function_(result_data_ptr, input_data_ptr));
            }
                //Judge Whether The Whole Computation Finished
            else {
                if (idle_count_ == thread_count_ - 1) {
                    is_end_of_reduce_ = true;
                    pthread_cond_broadcast(&task_taking_cond_var_);
                } else {
                    idle_count_++;
                    //Idle Worker, Go to Cond Var Buffer
                    while (pthread_cond_wait(&task_taking_cond_var_, &task_taking_mutex_) != 0);
                }
                pthread_mutex_unlock(&task_taking_mutex_);
            }
        }
#else
        //Send Back And Ask Corresponding One to Finish
        pthread_mutex_lock(&task_taking_mutex_);
        second_phase_global_reduce_data_vector.resize(thread_count_);
        pthread_mutex_unlock(&task_taking_mutex_);
        second_phase_global_reduce_data_vector[thread_idx] = std::move(
                first_phase_reduce_data_pool_vec_[reduce_data_indices_vec_[thread_idx].res_idx_]);
        pthread_barrier_wait(&timestamp_barrier_);
        cout << "Finished Send Back 2 Global" << endl;
#endif
    }

    template<typename DataCollectionType, typename DataType, typename DataCmpFunctionType, typename ComputationFunctionType>
    void *
    Reducer<DataCollectionType, DataType, DataCmpFunctionType, ComputationFunctionType>::InvokeRingCommThreadFunction(
            void *bundle_input_ptr) {
        auto my_bundle_input_ptr = ((BundleInput *) bundle_input_ptr);
        my_bundle_input_ptr->reducer_ptr_->RingCommTaskStealAndRequestThreadFunction(my_bundle_input_ptr->thread_id_);
        return NULL;
    }

    template<typename DataCollectionType, typename DataType, typename DataCmpFunctionType, typename ComputationFunctionType>
    unique_ptr<DataType>
    Reducer<DataCollectionType, DataType, DataCmpFunctionType, ComputationFunctionType>::ParallelExecute() {
        if (is_reduce_task_only_one_) {
            return std::move(first_phase_reduce_data_pool_vec_[0]);
        } else {
            std::vector<BundleInput *> input_bundle_vec(thread_count_);
            for (auto thread_id = 0; thread_id < thread_count_; thread_id++) {
                input_bundle_vec[thread_id] = new BundleInput();
                input_bundle_vec[thread_id]->reducer_ptr_ = this;
                input_bundle_vec[thread_id]->thread_id_ = thread_id;
                pthread_create(&thread_handles_[thread_id], NULL, this->InvokeRingCommThreadFunction,
                               (void *) input_bundle_vec[thread_id]);
            }

            for (auto thread_id = 0; thread_id < thread_count_; thread_id++) {
                pthread_join(thread_handles_[thread_id], NULL);
            }

            //Delete After All Execution-Flow Join
            for (auto i = 0; i < thread_count_; ++i) {
                delete input_bundle_vec[i];
            }

#ifdef REDUCE_2ND_PHASE_SEQUENTIAL
            cout << "Do Left Reduce In Single Core" << endl;
            //Do Left Reduce
            while (second_phase_global_reduce_data_vector.size() > 1) {
                cout << second_phase_global_reduce_data_vector.size() << endl;
                second_phase_global_reduce_data_vector[0] =
                        std::move(reduce_compute_function_(second_phase_global_reduce_data_vector[0],
                                                           second_phase_global_reduce_data_vector.back()));
                second_phase_global_reduce_data_vector.erase(second_phase_global_reduce_data_vector.end()-1);
            }

#endif
            return std::move(second_phase_global_reduce_data_vector[0]);
        }
    }
}

#endif //CODES_YCHE_REDUCER_H
