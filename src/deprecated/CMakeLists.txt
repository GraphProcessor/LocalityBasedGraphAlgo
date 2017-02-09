set(ParallelDetail parallel_utils/detail/dataflow_scheduler.h parallel_utils/detail/fine_grained_merge_scheduler.h parallel_utils/detail/reduce_scheduler.h)
set(ThreadPoolBase parallel_utils/detail/thread_pool_base.h)
set(ThreadPoolBreak parallel_utils/detail/thread_pool_breakable.h)

add_executable(ThreadPoolTest ${ThreadPoolBase} test/thread_pool_test.cpp)
target_link_libraries(ThreadPoolTest ${CMAKE_THREAD_LIBS_INIT} ${Boost_LIBRARIES})

add_executable(ThreadPoolBreakableTest ${ThreadPoolBase} ${ThreadPoolBreak} ${ParallelDetail} test/thread_pool_breakable_test.cpp)
target_link_libraries(ThreadPoolBreakableTest ${CMAKE_THREAD_LIBS_INIT} ${Boost_LIBRARIES})

