cmake_minimum_required(VERSION 2.8)
project(Codes_Yche)

#compiler flag
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -O3 -g")

#config some env var with cmake utils
find_package(Boost REQUIRED COMPONENTS regex thread)
find_package(Threads)

#alert some path information
message(${Boost_INCLUDE_DIRS})
message(${Boost_LIBRARIES})
message(${CMAKE_THREAD_LIBS_INIT})

include_directories(${Boost_INCLUDE_DIRS})
include_directories(".")

add_subdirectory(algorithm_demo)
