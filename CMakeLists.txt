cmake_minimum_required(VERSION 3.2)
project(lormap)

set(lormap_VERSION_MAJOR 1)
set(lormap_VERSION_MINOR 0)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(KSW2_FILES  src/ksw2_extd2_sse.c src/ksw2_extf2_sse.c src/ksw2_extz2_sse.c)

configure_file(
  "${PROJECT_SOURCE_DIR}/src/lormap.hpp.in"
  "${PROJECT_BINARY_DIR}/lormap.hpp"
)

add_library(brown_minimizers STATIC src/brown_minimizers.cc)

include_directories("${PROJECT_BINARY_DIR}")

add_subdirectory(vendor/bioparser EXCLUDE_FROM_ALL)
add_subdirectory(vendor/thread_pool EXCLUDE_FROM_ALL)

include_directories(${PROJECT_SOURCE_DIR}/src)

add_executable(lormap src/lormap.cc src/index.cc src/util.cc ${KSW2_FILES})

target_link_libraries(lormap bioparser brown_minimizers thread_pool pthread)

