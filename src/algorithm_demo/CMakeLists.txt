add_executable(demo_cis demo_cis.cpp ../algorithm/sequential/cis_sequential_algorithm.cpp ../util/graph_io_helper.cpp)
target_link_libraries(demo_cis ${Boost_LIBRARIES})

add_executable(demo_demon demo_demon.cpp ../algorithm/sequential/demon_sequential_algorithm.cpp ../util/graph_io_helper.cpp)
target_link_libraries(demo_demon ${Boost_LIBRARIES})

add_executable(demo_hk_growth demo_hk_growth.cpp ../algorithm/sequential/hk_grow_sequential_algorithm.cpp ../util/graph_io_helper.cpp)
target_link_libraries(demo_hk_growth ${Boost_LIBRARIES})