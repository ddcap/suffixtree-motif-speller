cmake_minimum_required(VERSION 2.6)
project(motifIterator)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(motifIterator main.cpp motif.cpp genefamily.cpp suffixtree.cpp)
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -O3 -Wall -pedantic -mpopcnt")

#set(CMAKE_BUILD_TYPE Debug)

install(TARGETS motifIterator RUNTIME DESTINATION bin)
