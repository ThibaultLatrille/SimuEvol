cmake_minimum_required(VERSION 3.7)
project(SimuEvol)

set(CMAKE_CXX_STANDARD 14)

# Compilation options
option(COVERAGE_MODE "For coverage mode using g++ " OFF) #OFF by default
option(DEBUG_MODE "Debug mode (with asserts and such) " OFF) #OFF by default
set(ALWAYSON_CXX_FLAGS "--std=c++14 -Wall -Wextra -Wno-unused-parameter -Wfatal-errors ")
if(COVERAGE_MODE)
    set(CMAKE_CXX_FLAGS "-O0 -fprofile-arcs -ftest-coverage ${ALWAYSON_CXX_FLAGS}") # coverage mode
    message("-- INFO: Compiling in coverage mode.\n-- INFO: flags are: " ${CMAKE_CXX_FLAGS})
elseif(DEBUG_MODE)
    set(CMAKE_CXX_FLAGS "-O0 -g ${ALWAYSON_CXX_FLAGS}") # debug mode
    message("-- INFO: Compiling in debug mode.\n-- INFO: flags are: " ${CMAKE_CXX_FLAGS})
else()
    set(CMAKE_CXX_FLAGS "-O3 -DNDEBUG ${ALWAYSON_CXX_FLAGS}") # release mode
    message("-- INFO: Compiling in release mode.\n-- INFO: flags are: " ${CMAKE_CXX_FLAGS})
endif(COVERAGE_MODE)
unset(COVERAGE_MODE CACHE)

include_directories("src")
include_directories("src/lib")
include_directories("utils")

add_executable(SimuDiv src/simu_div.cpp)
add_executable(SimuDfe src/simu_dfe.cpp)
add_executable(SimuGeo src/simu_geo.cpp)
add_executable(SimuFold src/simu_fold.cpp)
add_executable(SimuRelax src/simu_relax.cpp)
add_executable(SimuStab src/simu_stab.cpp)
add_executable(SimuPoly src/simu_poly.cpp)
add_executable(ToyGeo src/toy_geometric.cpp)
add_executable(ToyStab src/toy_stability.cpp)
