cmake_minimum_required(VERSION 2.6)
project(suffixtree)

add_executable(suffixtree main.cpp suffixtree.cpp genefamily.cpp motif.cpp)
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall -pedantic -mpopcnt")

# set(CMAKE_BUILD_TYPE Debug)

install(TARGETS suffixtree RUNTIME DESTINATION bin)
