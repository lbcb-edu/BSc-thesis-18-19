cmake_minimum_required(VERSION 2.8)
project(zavrsni)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -march=native")
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)


add_subdirectory(vendor/bioparser EXCLUDE_FROM_ALL)

add_executable(containers src/containers.cpp)
add_executable(detection src/detection.cpp)

target_link_libraries(containers bioparser)
target_link_libraries(detection bioparser)

