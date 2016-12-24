set(SOURCE_FILES modularity.cpp modularity.h include_header.h input_output_handler.h)
add_executable(ocd_evaluation ${SOURCE_FILES})
target_link_libraries(ocd_evaluation ${Boost_LIBRARIES})
