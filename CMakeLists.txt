cmake_minimum_required (VERSION 3.2)
project (MapToGraph)

set (sequence_to_graph_mapping_VERSION_MAJOR 1)
set (sequence_to_graph_mapping_VERSION_MINOR 0)

set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

option (sequence_to_graph_mapping_build_tests "Build sequence_to_graph unit tests" OFF)

configure_file (
        "${PROJECT_SOURCE_DIR}/src/sequence_to_graph_mapping.h.in"
        "${PROJECT_BINARY_DIR}/bin/sequence_to_graph_mapping.h"
)

#add_library(white_alignment STATIC src/white_alignment.cxx)
add_library(white_minimizers STATIC src/white_minimizers.cxx)

include_directories("${PROJECT_BINARY_DIR}/bin")

add_subdirectory (vendor/bioparser EXCLUDE_FROM_ALL)
add_subdirectory(vendor/thread_pool EXCLUDE_FROM_ALL)

add_executable (maptograph src/sequence_to_graph_mapping.cxx)
target_link_libraries (maptograph bioparser white_minimizers thread_pool pthread)

#if (sequence_to_graph_mapping_build_tests)
#        include_directories(${PROJECT_SOURCE_DIR}/src)
#
#        add_executable (sequence_to_graph_mapping_test test/sequence_to_graph_mapping_test.cxx)
#
#        add_subdirectory (vendor/googletest/googletest EXCLUDE_FROM_ALL)
#
#        target_link_libraries (sequence_to_graph_mapping_test gtest_main)
#endif (sequence_to_graph_mapping_tests)
