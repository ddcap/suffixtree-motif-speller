cmake_minimum_required(VERSION 2.6)
project(suffixtree)

add_executable(suffixtree main.cpp motif.cpp genefamily.cpp suffixtree.cpp)
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -O3 -Wall -pedantic -mpopcnt")

set(CMAKE_BUILD_TYPE Debug)

install(TARGETS suffixtree RUNTIME DESTINATION bin)
