add_executable(demo_cis demo_cis.cpp ../algorithm/sequential/cis_sequential_algorithm.cpp ../util/basic_io_helper.h)
target_link_libraries(demo_cis ${Boost_LIBRARIES})

add_executable(demo_demon demo_demon.cpp ../algorithm/sequential/demon_seqential_algorithm.cpp ../algorithm/sequential/demon_sequential_algorithm.h)
target_link_libraries(demo_demon ${Boost_LIBRARIES})

add_executable(demo_hr_growth demo_hr_growth.cpp ../algorithm/sequential/hr_grow_sequential_algorithm.cpp ../algorithm/sequential/hr_grow_sequential_algorithm.h)
target_link_libraries(demo_demon ${Boost_LIBRARIES})