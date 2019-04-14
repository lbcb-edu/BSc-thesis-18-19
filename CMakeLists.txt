cmake_minimum_required(VERSION 3.2)
project(long_reads_mapper)

set(long_reads_mapper_VERSION_MAJOR 1)
set(long_reads_mapper_VERSION_MINOR 1)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

configure_file(
  "${PROJECT_SOURCE_DIR}/src/long_reads_mapper.hpp.in"
  "${PROJECT_BINARY_DIR}/long_reads_mapper.hpp"
)
add_library(brown_minimizers STATIC src/brown_minimizers.cc)

include_directories("${PROJECT_BINARY_DIR}")

add_subdirectory(vendor/bioparser EXCLUDE_FROM_ALL)
add_subdirectory(vendor/thread_pool EXCLUDE_FROM_ALL)
AUX_SOURCE_DIRECTORY(vendor/ksw2 VARIABLE) 

include_directories(${PROJECT_SOURCE_DIR}/src)

add_executable(long_reads_mapper src/long_reads_mapper.cc)

target_link_libraries(long_reads_mapper bioparser brown_minimizers thread_pool pthread)

