cmake_minimum_required(VERSION 2.8)
project(libmsr)

set(CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} -std=c99")

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
add_library(msr src/mbr_encode.c src/msr_encode.c src/gf.c)

set(TEST_FILES test/test.c include/gf.h)

add_executable(test ${TEST_FILES})
target_link_libraries(test msr)