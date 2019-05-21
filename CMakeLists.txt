cmake_minimum_required (VERSION 3.2)
project(OptSeqAlignmentLongGaps)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

#version
set (OptSeqAlignmentLongGaps_VERSION_MAJOR 0)
set (OptSeqAlignmentLongGaps_VERSION_MINOR 9)
set (OptSeqAlignmentLongGaps_VERSION_PATCH 0)

set(CMAKE_CXX_STANDARD 11)

configure_file (
    "${PROJECT_SOURCE_DIR}/src/OptSeqAlignmentLongGapsConfig.h.in"
    "${PROJECT_BINARY_DIR}/OptSeqAlignmentLongGapsConfig.h"
)

include_directories("${PROJECT_BINARY_DIR}")

add_subdirectory(vendor/bioparser EXCLUDE_FROM_ALL)

add_library(OSALG_lib ${PROJECT_SOURCE_DIR}/src/OSALG_lib.cpp)
add_library(OSALG_lib_vector ${PROJECT_SOURCE_DIR}/src/OSALG_lib_vector.cpp)

add_executable(OptSeqAlignmentLongGaps src/OptSeqAlignmentLongGaps.cpp)
target_link_libraries(OptSeqAlignmentLongGaps bioparser OSALG_lib OSALG_lib_vector)
