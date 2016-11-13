cmake_minimum_required(VERSION 2.8)
project(libmsr)

set(CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} -std=c99 -O3 -mavx2")

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
add_library(msr src/mbr_encode.c src/gf.c src/msr_encode.c)

set(TEST_FILES test/test.c)
set(BENCHMARK_FILES test/benchmark.c)

add_executable(test ${TEST_FILES})

add_executable(benchmark ${BENCHMARK_FILES})


target_link_libraries(test msr)
target_link_libraries(benchmark msr)