cmake_minimum_required(VERSION 3.2)
project(srmapper)

set(srmapper_VERSION_MAJOR 1)
set(srmapper_VERSION_MINOR 0)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(mapper_build_tests "Build mapper unit tests" OFF)

configure_file(
  "${PROJECT_SOURCE_DIR}/src/srmapper.hpp.in"
  "${PROJECT_BINARY_DIR}/srmapper.hpp"
)

add_library(brown_minimizers STATIC
    src/brown_minimizers.cc)

include_directories("${PROJECT_BINARY_DIR}")

add_subdirectory(vendor/bioparser EXCLUDE_FROM_ALL)
add_subdirectory(vendor/thread_pool EXCLUDE_FROM_ALL)

if (mapper_build_tests)
  include_directories(${PROJECT_SOURCE_DIR}/src)
  add_executable(srmapper_unittest test/srmapper_unittest.cc)
  add_subdirectory(vendor/googletest/googletest EXCLUDE_FROM_ALL)
  target_link_libraries(srmapper_unittest gtest_main)
endif (mapper_build_tests)

add_executable(srmapper src/srmapper.cc src/brown_minimizers.cc)
target_link_libraries(srmapper brown_minimizers bioparser thread_pool pthread)
