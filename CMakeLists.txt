cmake_minimum_required(VERSION 2.8)
project(zavrsni)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -march=native")
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_subdirectory(vendor/spoa EXCLUDE_FROM_ALL)
add_subdirectory(vendor/bioparser EXCLUDE_FROM_ALL)

add_library(blue_alignment STATIC src/blue_alignment.cpp)

add_executable(zavrsni src/zavrsni.cpp)

target_link_libraries(zavrsni spoa bioparser blue_alignment)


