cmake_minimum_required(VERSION 2.6)
project(motifIterator)
#set(CMAKE_PROJECT_VERSION 2)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
add_subdirectory(sparse-map)

add_executable(motifIterator main.cpp motif.cpp genefamily.cpp suffixtree.cpp motifmap.cpp)
target_link_libraries(motifIterator PRIVATE tsl::sparse_map)

SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -O3 -Wall -pedantic -mpopcnt")

#set(CMAKE_BUILD_TYPE Debug)

install(TARGETS motifIterator RUNTIME DESTINATION bin)
